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

// Two-medium image formation, using analytic signals (each sample is a real-imag vector).
// The final image is a combination of sub-images, using apodization.
//
// Processing steps:
//   Find row at the interface                                             - CPU
//   Calculate delays at the interface                                     - CPU
//   Calculate delays in the grid above the interface                      - CUDA
//   Signal preparation                                                    - CPU
//   Apply delays, use apodization, [apply PCF] and accumulate the samples - CUDA
//
// The grid must be rectangular.
//
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

	const TwoMediumSTAConfiguration<float>& config_;
	std::vector<Matrix<float>>& acqDataList_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<float>& coherenceFactor_;
	float maxFermatBlockSize_;
	const float lambda2_;
	unsigned int signalLength_;
	float signalOffset_;
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> minRowIdx_; // for each column
	std::vector<float, tbb::cache_aligned_allocator<float>> xArray_;
	Matrix<float, tbb::cache_aligned_allocator<float>> medium1DelayMatrix_; // (interface_idx, element)
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData<float>>> prepareDataTLS_;
	std::unique_ptr<VectorialCombinedTwoMediumImagingCUDAProcessor2Data> data_;

	unsigned int stepConfigListSize_;
	unsigned int interfacePointListSize_;
	unsigned int rxApodSize_;
	unsigned int gridXZN1_;
	unsigned int gridXZN2_;
};

} // namespace Lab

#endif // VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_CUDA_PROCESSOR_2_H
