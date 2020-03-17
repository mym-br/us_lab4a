/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_2_H
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_2_H

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



namespace Lab {

struct VectorialCombinedTwoMediumImagingCUDAProcessor2Data;

// Without PCF, this class is slower than VectorialCombinedTwoMediumImagingProcessor,
// using a Core i5-3470. With PCF, this class is faster.
//
// Uses CUDA in part of the processing.
//
// The grid must be rectangular.
// Requirements:
// - Single precision
//
// numElements must be multiple of 16.
// Tested only with numElements = 32.
class VectorialCombinedTwoMediumImagingCUDAProcessor2 {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int txElem;
	};

	template<typename TFloat>
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};

	VectorialCombinedTwoMediumImagingCUDAProcessor2(
			const TwoMediumSTAConfiguration<float>& config,
			std::vector<Matrix<float>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<float>& coherenceFactor,
			float maxFermatBlockSize,
			float peakOffset,
			unsigned int signalStartOffset);
	~VectorialCombinedTwoMediumImagingCUDAProcessor2();

	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<float>>& interfacePointList,
		const std::vector<float>& rxApod,
		const Matrix<XZ<float>>& gridXZ,
		Matrix<std::complex<float>>& gridValue);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tMinRowIdx;
	MeasurementList<double> tMedium1DelayMatrix;
	MeasurementList<double> tCalculateDelays;
	MeasurementList<double> tPrepareData;
	MeasurementList<double> tProcessColumn;
#endif

private:
	VectorialCombinedTwoMediumImagingCUDAProcessor2(const VectorialCombinedTwoMediumImagingCUDAProcessor2&) = delete;
	VectorialCombinedTwoMediumImagingCUDAProcessor2& operator=(const VectorialCombinedTwoMediumImagingCUDAProcessor2&) = delete;
	VectorialCombinedTwoMediumImagingCUDAProcessor2(VectorialCombinedTwoMediumImagingCUDAProcessor2&&) = delete;
	VectorialCombinedTwoMediumImagingCUDAProcessor2& operator=(VectorialCombinedTwoMediumImagingCUDAProcessor2&&) = delete;

	static std::size_t roundUpToMultipleOfGroupSize(std::size_t n, std::size_t groupSize) {
		std::size_t numGroups = (n + (groupSize - 1)) / groupSize;
		return numGroups * groupSize;
	}

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
	Matrix<float, tbb::cache_aligned_allocator<float>> medium1DelayMatrix_; // (interface_idx, element)
	Tensor3<float, tbb::cache_aligned_allocator<float>> delayTensor_;
	unsigned int rawDataN1_;
	unsigned int rawDataN2_;
	std::size_t rawDataSize_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData<float>>> prepareDataTLS_;
	std::unique_ptr<VectorialCombinedTwoMediumImagingCUDAProcessor2Data> data_;
};

} // namespace Lab

#endif // VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_2_H
