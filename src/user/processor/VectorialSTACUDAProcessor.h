/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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

#ifndef VECTORIAL_STA_CUDA_PROCESSOR_H
#define VECTORIAL_STA_CUDA_PROCESSOR_H

#include <memory>
#include <vector>

#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>

#include "ArrayProcessor.h"
#include "CoherenceFactor.h"
#include "ExecutionTimeMeasurement.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "XYZValueFactor.h"



namespace Lab {

struct VectorialSTACUDAProcessorData;

// STA image formation, using analytic signals (each sample is a real-imag vector).
//
// Processing steps:
//   Signal preparation                                                    - CPU
//   Calculate delays                                                      - CUDA
//   Apply delays, use apodization, [apply PCF] and accumulate the samples - CUDA
//
// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
class VectorialSTACUDAProcessor : public ArrayProcessor<XYZValueFactor<float>> {
public:
	struct PrepareDataThreadData {
		Interpolator<float> interpolator;
		HilbertEnvelope<float> envelope;
		std::vector<float, tbb::cache_aligned_allocator<float>> signal;
	};

	VectorialSTACUDAProcessor(
			const STAConfiguration<float>& config,
			STAAcquisition<float>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<float>& coherenceFactor,
			float peakOffset,
			bool calculateEnvelope,
			const std::vector<float>& txApod,
			const std::vector<float>& rxApod);
	virtual ~VectorialSTACUDAProcessor();

	virtual void prepare(unsigned int baseElement);
	virtual void process(Matrix<XYZValueFactor<float>>& gridData);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tAcquisitionML;
	MeasurementList<double> tPrepareDataML;
	MeasurementList<double> tCalculateDelaysML;
	MeasurementList<double> tDelaySumML;
	void execTimeMeasReset(unsigned int n) {
		tAcquisitionML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tPrepareDataML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tCalculateDelaysML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tDelaySumML.reset(       EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	void execTimeMeasShowResults(unsigned int n) {
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tAcquisition:    ", tAcquisitionML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:    ", tPrepareDataML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "tCalculateDelays:", tCalculateDelaysML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "tDelaySum:       ", tDelaySumML);
	}
#endif

private:
	VectorialSTACUDAProcessor(const VectorialSTACUDAProcessor&) = delete;
	VectorialSTACUDAProcessor& operator=(const VectorialSTACUDAProcessor&) = delete;
	VectorialSTACUDAProcessor(VectorialSTACUDAProcessor&&) = delete;
	VectorialSTACUDAProcessor& operator=(VectorialSTACUDAProcessor&&) = delete;

	const STAConfiguration<float>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<float>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<float>& coherenceFactor_;
	typename STAAcquisition<float>::AcquisitionDataType acqData_;
	unsigned int signalLength_;
	float signalOffset_;
	std::vector<float, tbb::cache_aligned_allocator<float>> xArray_;
	bool calculateEnvelope_;//TODO: Remove?
	bool initialized_;
	std::vector<float> txApod_;
	std::vector<float> rxApod_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<VectorialSTACUDAProcessorData> data_;

	unsigned int gridN1_;
	unsigned int gridN2_;
};

} // namespace Lab

#endif // VECTORIAL_STA_CUDA_PROCESSOR_H
