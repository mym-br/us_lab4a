/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020, 2025 Marcelo Y. Matuda               *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/

#ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_H
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_H

#include <complex>
#include <cstddef> /* std::size_t */
#include <memory>
#include <vector>

#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>

#include "CoherenceFactor.h"
#include "ExecutionTimeMeasurement.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Matrix.h"
#include "Tensor3.h"
#include "TwoMediumSTAConfiguration.h"
#include "XZ.h"

class dim3;

namespace Lab {

void execTransposeKernel(const dim3& gridDim, const dim3& blockDim,
			float* rawData, float* rawDataT, unsigned int oldSizeX, unsigned int oldSizeY);
void execProcessImageKernel(const dim3& gridDim, const dim3& blockDim,
			float* rawData, unsigned int numGridPoints, float* gridValueRe,
			float* gridValueIm, float* rxApod);
void execProcessImagePCFKernel(const dim3& gridDim, const dim3& blockDim,
			float* rawData, unsigned int numGridPoints, float* gridValueRe,
			float* gridValueIm, float* rxApod, float pcfFactor);



// Two-medium image formation, using analytic signals (each sample is a real-imag vector).
// The final image is a combination of sub-images, using apodization.
//
// Without PCF, this class is slower than VectorialCombinedTwoMediumImagingProcessor,
// using a Core i5-3470. With PCF, this class is faster.
//
// Processing steps:
//   Find row at the interface                               - CPU
//   Calculate delays at the interface                       - CPU
//   Calculate delays in the grid above the interface        - CPU
//   Signal preparation                                      - CPU
//   Apply delays and store the sample                       - CPU
//   Use apodization, [apply PCF] and accumulate the samples - CUDA
//
// The grid must be rectangular.
//
class VectorialCombinedTwoMediumImagingCUDAProcessor {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int txElem;
	};

	VectorialCombinedTwoMediumImagingCUDAProcessor(
			const TwoMediumSTAConfiguration<float>& config,
			std::vector<Matrix<float>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<float>& coherenceFactor,
			float maxFermatBlockSize,
			float peakOffset,
			unsigned int signalStartOffset);
	~VectorialCombinedTwoMediumImagingCUDAProcessor();

	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<float>>& interfacePointList,
		const std::vector<float>& rxApod,
		const Matrix<XZ<float>>& gridXZ,
		Matrix<std::complex<float>>& gridValue);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tMinRowIdxML;
	MeasurementList<double> tMedium1DelayMatrixML;
	MeasurementList<double> tCalculateDelaysML;
	MeasurementList<double> tPrepareDataML;
	MeasurementList<double> tProcessColumnML;
	void execTimeMeasReset() {
		tMinRowIdxML.reset(         EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tMedium1DelayMatrixML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tCalculateDelaysML.reset(   EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tPrepareDataML.reset(       EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tProcessColumnML.reset(     EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	void execTimeMeasShowResults() {
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMinRowIdx:         ", tMinRowIdxML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMedium1DelayMatrix:", tMedium1DelayMatrixML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tCalculateDelays:   ", tCalculateDelaysML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tPrepareData:       ", tPrepareDataML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tProcessColumn:     ", tProcessColumnML);
	}
#endif

private:
	struct CUDAData;

	template<typename TFloat>
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};

	template<typename TFloat> struct CalculateDelays;
	template<typename TFloat> struct PrepareDataWithOneTxElem;
	template<typename TFloat> struct ProcessColumnWithOneTxElem;

	VectorialCombinedTwoMediumImagingCUDAProcessor(const VectorialCombinedTwoMediumImagingCUDAProcessor&) = delete;
	VectorialCombinedTwoMediumImagingCUDAProcessor& operator=(const VectorialCombinedTwoMediumImagingCUDAProcessor&) = delete;
	VectorialCombinedTwoMediumImagingCUDAProcessor(VectorialCombinedTwoMediumImagingCUDAProcessor&&) = delete;
	VectorialCombinedTwoMediumImagingCUDAProcessor& operator=(VectorialCombinedTwoMediumImagingCUDAProcessor&&) = delete;

	const TwoMediumSTAConfiguration<float>& config_;
	std::vector<Matrix<float>>& acqDataList_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<float>& coherenceFactor_;
	float maxFermatBlockSize_;
	const float lambda2_;
	std::size_t signalLength_;
	Tensor3<std::complex<float>, tbb::cache_aligned_allocator<std::complex<float>>> signalTensor_;
	float signalOffset_;
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> minRowIdx_; // for each column
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> firstGridPointIdx_; // for each column
	std::vector<float, tbb::cache_aligned_allocator<float>> xArray_;
	Matrix<float, tbb::cache_aligned_allocator<float>> medium1DelayMatrix_;
	Tensor3<float, tbb::cache_aligned_allocator<float>> delayTensor_;
	unsigned int rawDataN1_;
	unsigned int rawDataN2_;
	std::size_t rawDataSize_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData<float>>> prepareDataTLS_;
	std::unique_ptr<CUDAData> data_;

	unsigned int rxApodSize_;
	unsigned int gridXZN1_;
	unsigned int gridXZN2_;
};

} // namespace Lab

#endif // VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_H
