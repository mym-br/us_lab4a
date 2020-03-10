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

#ifndef CYLINDERDETECTIONANDFERMATMETHOD_H_
#define CYLINDERDETECTIONANDFERMATMETHOD_H_

#include <algorithm> /* copy, min, reverse */
#include <cmath>
#include <cstddef> /* std::size_t */
#include <iomanip> /* setprecision */
#include <sstream>
#include <type_traits> /* is_same */
#include <utility> /* make_pair, pair */
#include <vector>

#include "ArcCylindricalWaveProcessor.h"
#include "ArrayProcessor.h"
#include "CCBFPitchCatchProcessor.h"
#include "CoherenceFactor.h"
#include "CrossCorrArcCylindricalWaveProcessor.h"
#include "CrossCorrelationSTAProcessor.h"
#include "CylindricalWaveProcessor.h"
#include "DefaultSTAProcessor.h"
#include "ExecutionTimeMeasurement.h"
#include "FFTWFilter.h"
#include "Fitting.h"
#include "HilbertEnvelope.h"
#include "ImageGrid.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NetworkSTAAcquisition.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SavedSTAAcquisition.h"
#include "SimulatedSTAAcquisition.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "String.h"
#include "TangentCurveGeometryProcessor.h"
#include "Tensor3.h"
#include "Timer.h"
#include "TwoMediumSTAConfiguration.h"
#include "TwoMediumSTAProcessor.h"
#include "Util.h"
#include "VectorialCombinedTwoMediumImagingProcessor.h"
#include "VectorialSTAProcessor.h"
#include "VectorialTwoMediumSTAProcessor.h"
#include "XZComplexValue.h"
#include "XZComplexValueFactor.h"
#include "XZValueFactor.h"
#ifdef USE_OPENCL
# include "VectorialCombinedTwoMediumImagingOCLProcessor.h"
#endif
#ifdef USE_CUDA
# include "VectorialCombinedTwoMediumImagingCUDAProcessor.h"
#endif

#define CYL_DETECT_AND_FERMAT_METHOD_SCF_LOW_PASS_LAMBDAS 2.0

#define CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_GET_ENVELOPE 1
#define CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_USE_SCF 1
//#define CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA 1
//#define CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_PEAK_IMAGES 1
#define CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA 1
//#define CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_ARCS 1

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_GEOMETRY_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_USE_CROSS_CORRELATION 1

#define CYLINDER_DETECTION_AND_FERMAT_METHOD_PITCH_CATCH_USE_CROSS_CORRELATION 1
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_TANGENT_CURVE_USE_CROSS_CORRELATION 1
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_ARCS_USE_CROSS_CORRELATION 1
//#define CYLINDER_DETECTION_AND_FERMAT_METHOD_ARC_USE_ANALYTIC_SIGNAL 1

#define CYLINDER_DETECTION_AND_FERMAT_METHOD_ARC_USE_LOBE_CENTER 1
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_ARC_LOBE_THRESHOLD (0.5)

#define CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION 1
#define CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_USE_CROSS_CORRELATION 1



namespace Lab {

template<typename TFloat>
class CylinderDetectionAndFermatMethod : public Method {
public:
	CylinderDetectionAndFermatMethod(Project& project) : project_(project) { }
	virtual ~CylinderDetectionAndFermatMethod() = default;

	virtual void execute();
private:
	CylinderDetectionAndFermatMethod(const CylinderDetectionAndFermatMethod&) = delete;
	CylinderDetectionAndFermatMethod& operator=(const CylinderDetectionAndFermatMethod&) = delete;
	CylinderDetectionAndFermatMethod(CylinderDetectionAndFermatMethod&&) = delete;
	CylinderDetectionAndFermatMethod& operator=(CylinderDetectionAndFermatMethod&&) = delete;

	void acquireSignals();
	// Detects points using the peaks in STA images.
	void detectPoints();
	// Detects points using the peaks in STA images. The images are formed using sums of analytic signals.
	void detectPointsUsingVectorSum();
	// Detects points using the peaks in STA images. The images are formed using sums of the results of cross-correlations.
	void detectPointsUsingCrossCorrelation();
	// Detects points using geometry. The distances are calculated using one element at a time, as emitter and receiver.
	void detectPointsUsingCCBFPulseEcho();
	// Detects points using geometry. One element emits and others receive.
	void detectPointsUsingCCBFPitchCatch();
	// Detects points using arcs instead of images.
	void detectPointsInArcs();
	// Detects points using only geometry (envelope of arc and ellipses).
	void detectPointsUsingTangentCurveGeometry();
	void fitCircle();
	void execTwoMediumImaging();
	void execCombinedTwoMediumImaging();
	template<typename Proc> void execCombinedTwoMediumImagingCyl();

	void measureSpeed1();
	void measureSpeed1AndDistanceError();
	void measureSpeed1AndDistanceError2();
	void measureDistance();

	Project& project_;
};



template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::acquireSignals()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir      = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto numBaseElemSteps = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition;

	switch (project_.method()) {
	case MethodEnum::cylinder_detection_and_fermat_acquisition_simulated_dp:
	case MethodEnum::cylinder_detection_and_fermat_acquisition_simulated_sp:
		acquisition = std::make_unique<SimulatedSTAAcquisition<TFloat>>(project_, config);
		break;
	case MethodEnum::cylinder_detection_and_fermat_acquisition_network_dp:
	case MethodEnum::cylinder_detection_and_fermat_acquisition_network_sp:
		acquisition = std::make_unique<NetworkSTAAcquisition<TFloat>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;

	Timer t0;

	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		acquisition->prepare(baseElement);
		for (unsigned int txElem = 0; txElem < config.numElements; ++txElem) {
			acquisition->execute(txElem, acqData);
			project_.saveHDF5(acqData, String::Begin() << savedAcqDir << "/base" << baseElement << "_s" << txElem << String::End(), "signal");
		}
	}

	LOG_DEBUG << "ACQ time = " << t0.getTime();
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPoints()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir      = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir        = methodPM->value<std::string>( "output_dir");
	const auto numBaseElemSteps = methodPM->value<unsigned int>("num_base_elem_steps",   3, config.numElementsMux - config.numElements + 1);
	const auto peakOffset       = methodPM->value<TFloat>(      "peak_offset"        , 0.0, 50.0);

	project_.createDirectory(outputDir, false);
	Matrix<XZValueFactor<TFloat>> gridData;

	const TFloat lambda = config.propagationSpeed / config.centerFrequency;
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), lambda, gridData);

	CoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	std::unique_ptr<ArrayProcessor<XZValueFactor<TFloat>>> p =
		std::make_unique<DefaultSTAProcessor<TFloat, XZValueFactor<TFloat>>>(
				config, *acquisition, coherenceFactor, peakOffset);

	std::vector<std::pair<TFloat, TFloat>> peakPositionList(baseElemList.size()); // x, z
	HilbertEnvelope<TFloat> envelope;

	// Find the peaks.
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[PEAK] ---------- baseElement = " << baseElement;

		// Obtain the image.
		p->prepare(baseElement);
		p->process(gridData);

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the raw image... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_raw_image_" << baseElement << String::End(), "image", Util::CopyValueOp());

		LOG_DEBUG << "Saving the factors... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_factor_" << baseElement << String::End(), "factors", Util::CopyFactorOp());
#endif

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_USE_SCF
		// Apply the SCF.
		for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->value *= iter->factor;
			iter->factor = 1;
		}
# ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
# endif
#endif

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_GET_ENVELOPE
		// Get the envelope.
		for (unsigned int i = 0, end = gridData.n1(); i < end; ++i) {
			envelope.calculate(&gridData(i, 0), gridData.n2());
		}
# ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the image (after envelope)... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_env_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
# endif
#endif

		// Divide by r^2 and find the peak.
		TFloat maxAbsValue = -1;
		peakPositionList[baseElemIdx].first = 0;
		peakPositionList[baseElemIdx].second = 0;
		for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
			const TFloat x = iter->x;
			const TFloat z = iter->z;
			const TFloat r = std::sqrt(x * x + z * z);
			iter->value /= r * r;

			const TFloat absVal = std::abs(iter->value);
			if (absVal > maxAbsValue) {
				maxAbsValue = absVal;
				// Convert to the coordinate system centered on the first element.
				peakPositionList[baseElemIdx].first = iter->x + (baseElement + (config.numElements - 1) / TFloat(2)) * config.pitch;
				peakPositionList[baseElemIdx].second = iter->z;
			}
		}
		LOG_DEBUG << "Peak: x = " << peakPositionList[baseElemIdx].first << " z = " << peakPositionList[baseElemIdx].second;

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the modified image... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_env_modif_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
#endif

		if (project_.processingCancellationRequested()) {
			LOG_INFO  << "########## Processing cancelled.";
			return;
		}
	}

	LOG_DEBUG << "Saving the positions of the peaks...";
	project_.saveHDF5(peakPositionList, outputDir + "/point_pos_xz", "pointPositionList");

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, outputDir + "/peak_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, outputDir + "/peak_z", "z", Util::CopyZOp());
#endif
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPointsUsingVectorSum()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir      = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir        = methodPM->value<std::string>( "output_dir");
	const auto upsamplingFactor = methodPM->value<unsigned int>("upsampling_factor"  , 1, 128);
	const auto numBaseElemSteps = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);
	const auto peakOffset       = methodPM->value<TFloat>(      "peak_offset"        , 0.0, 50.0);

	project_.createDirectory(outputDir, false);
	Matrix<XZValueFactor<TFloat>> gridData;

	const TFloat lambda = config.propagationSpeed / config.centerFrequency;
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), lambda, gridData);

	AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');

	const std::vector<TFloat> txApod(config.numElements, 1.0);
	const std::vector<TFloat> rxApod(config.numElements, 1.0);
	std::unique_ptr<ArrayProcessor<XZValueFactor<TFloat>>> p =
		std::make_unique<VectorialSTAProcessor<TFloat, XZValueFactor<TFloat>>>(
					config,
					*acquisition,
					upsamplingFactor,
					coherenceFactor,
					peakOffset,
					true,
					txApod,
					rxApod);

	std::vector<std::pair<TFloat, TFloat>> peakPositionList(baseElemList.size()); // x, z
	std::vector<TFloat> auxPeakX(peakPositionList.size()), auxPeakZ(peakPositionList.size());

	// Find the peaks.
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[PEAK] ---------- baseElement = " << baseElement;

		// Obtain the image.
		p->prepare(baseElement);
		p->process(gridData);

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the raw image... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_raw_image_" << baseElement << String::End(), "image", Util::CopyValueOp());

		LOG_DEBUG << "Saving the factors... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_factor_" << baseElement << String::End(), "factors", Util::CopyFactorOp());
#endif
		if (coherenceFactor.enabled()) {
			// Apply the CF.
			for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
				iter->value *= iter->factor;
				iter->factor = 1;
			}
#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
			LOG_DEBUG << "Saving the image (after CF)... baseElement = " << baseElement;
			project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
#endif
		}

		// Divide by r^2 and find the peak.
		TFloat maxAbsValue = -1;
		peakPositionList[baseElemIdx].first = 0;
		peakPositionList[baseElemIdx].second = 0;
		for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
//			const TFloat x = iter->x;
//			const TFloat z = iter->z;
//			const TFloat r = std::sqrt(x * x + z * z);
//			iter->value /= r * r;

			const TFloat absVal = std::abs(iter->value);
			if (absVal > maxAbsValue) {
				maxAbsValue = absVal;
				// Convert to the coordinate system centered on the first element.
				peakPositionList[baseElemIdx].first = iter->x + (baseElement + (config.numElements - 1) / TFloat(2)) * config.pitch;
				peakPositionList[baseElemIdx].second = iter->z;
			}
		}

		LOG_DEBUG << "Peak: x = " << peakPositionList[baseElemIdx].first << " z = " << peakPositionList[baseElemIdx].second;
#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the modified image... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_env_modif_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
#endif
		for (unsigned int i = 0; i < peakPositionList.size(); ++i) {
			auxPeakX[i] = peakPositionList[i].first;
			auxPeakZ[i] = peakPositionList[i].second;
		}
		project_.showFigure2D(0, "Points", auxPeakX, auxPeakZ, true, true);

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_PEAK_IMAGES
		project_.showFigure3D(1, "Peak detection", &gridData, Project::emptyPointList,
					true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS);
#endif
		if (project_.processingCancellationRequested()) {
			LOG_INFO  << "########## Processing cancelled.";
			return;
		}
	}

	LOG_DEBUG << "Saving the positions of the peaks...";
	project_.saveHDF5(peakPositionList, outputDir + "/point_pos_xz", "pointPositionList");

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, outputDir + "/peak_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, outputDir + "/peak_z", "z", Util::CopyZOp());
#endif
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPointsUsingCrossCorrelation()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"  , 1, 128);
	const auto numBaseElemSteps   = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);

	project_.createDirectory(outputDir, false);
	Matrix<XZValueFactor<TFloat>> gridData;

	const TFloat lambda = config.propagationSpeed / config.centerFrequency;
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), lambda, gridData);

	CoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	std::unique_ptr<ArrayProcessor<XZValueFactor<TFloat>>> p =
		std::make_unique<CrossCorrelationSTAProcessor<TFloat>>(
						config,
						*acquisition,
						upsamplingFactor,
						coherenceFactor,
						refPulse);

	std::vector<std::pair<TFloat, TFloat>> peakPositionList(baseElemList.size()); // x, z
	std::vector<TFloat> auxPeakX(peakPositionList.size()), auxPeakZ(peakPositionList.size());

	// Find the peaks.
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[PEAK] ---------- baseElement = " << baseElement;

		// Obtain the image.
		p->prepare(baseElement);
		p->process(gridData);

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the raw image... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_raw_image_" << baseElement << String::End(), "image", Util::CopyValueOp());

		LOG_DEBUG << "Saving the factors... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_factor_" << baseElement << String::End(), "factors", Util::CopyFactorOp());
#endif

		if (coherenceFactor.enabled()) {
			// Apply the CF.
			for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
				iter->value *= iter->factor;
				iter->factor = 1;
			}
#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
			LOG_DEBUG << "Saving the image (after CF)... baseElement = " << baseElement;
			project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
#endif
		}

		// Find the peak.
		TFloat maxValue = -1;
		peakPositionList[baseElemIdx].first = 0;
		peakPositionList[baseElemIdx].second = 0;
		for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
//			const TFloat x = iter->x;
//			const TFloat z = iter->z;
//			const TFloat r = std::sqrt(x * x + z * z);

			// Divide by r^2.
//			iter->value /= r * r;

			if (iter->value > maxValue) {
				maxValue = iter->value;
				// Convert to the coordinate system centered on the first element.
				peakPositionList[baseElemIdx].first = iter->x + (baseElement + (config.numElements - 1) / TFloat(2)) * config.pitch;
				peakPositionList[baseElemIdx].second = iter->z;
			}
		}

		LOG_DEBUG << "Peak: x = " << peakPositionList[baseElemIdx].first << " z = " << peakPositionList[baseElemIdx].second;
#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
		LOG_DEBUG << "Saving the modified image... baseElement = " << baseElement;
		project_.saveHDF5(gridData, String::Begin() << outputDir << "/peak_cf_env_modif_image_" << baseElement << String::End(), "image", Util::CopyValueOp());
#endif
		for (unsigned int i = 0; i < peakPositionList.size(); ++i) {
			auxPeakX[i] = peakPositionList[i].first;
			auxPeakZ[i] = peakPositionList[i].second;
		}
		project_.showFigure2D(0, "Points", auxPeakX, auxPeakZ, true, true);

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_PEAK_IMAGES
		//project_.showFigure3D(baseElemIdx, "Peak detection", &gridData, &arc,
		project_.showFigure3D(baseElemIdx, "Peak detection", &gridData, Project::emptyPointList,
					true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS);
#endif
		if (project_.processingCancellationRequested()) {
			LOG_INFO  << "########## Processing cancelled.";
			return;
		}
	}

	LOG_DEBUG << "Saving the positions of the peaks...";
	project_.saveHDF5(peakPositionList, outputDir + "/point_pos_xz", "pointPositionList");
#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SAVE_DATA
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, outputDir + "/peak_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, outputDir + "/peak_z", "z", Util::CopyZOp());
#endif
}

// Camacho, J.
// Cruza, J. F.
// Brizuela, J.
// Fritsch, C.
// Automatic dynamic depth focusing for NDT.
// IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control,
// vol. 61, no. 4, pp. 673-684, 2014.
// DOI: 10.1109/TUFFC.2014.2955
template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPointsUsingCCBFPulseEcho()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"  , 2, 128);
	const auto numBaseElemSteps   = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset");
#endif

	project_.createDirectory(outputDir, false);

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	//==============================
	// Distance calculation using high sampling rate.
	//

	Interpolator<TFloat> interp;
	interp.prepare(upsamplingFactor, CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);

	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
	std::vector<TFloat> signal;
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_USE_CROSS_CORRELATION
	std::vector<TFloat> tempSignal;
	FFTWFilter<TFloat> revRefPulseFilter;

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor);
	interp.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
	std::reverse(revRefPulse.begin(), revRefPulse.end());
	revRefPulseFilter.setCoefficients(revRefPulse);

	const unsigned int signalOffset = revRefPulse.size() - 1; // cross-correlation using convolution (revRefPulseFilter)
#else
	HilbertEnvelope<TFloat> envelope;

	const unsigned int signalOffset = (config.samplingFrequency * upsamplingFactor) * peakOffset / config.centerFrequency;
#endif
	LOG_DEBUG << "signalOffset: " << signalOffset;

	std::vector<TFloat> distanceList(config.numElements);
	std::vector<std::pair<TFloat, TFloat>> pointPositionList(baseElemList.size() * (config.numElements - 1)); // x, z
	std::vector<TFloat> auxPointX(pointPositionList.size()), auxPointZ(pointPositionList.size());
	unsigned int deadZoneSamplesUp = ((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed);

	// Acquire an A-scan to obtain its length.
	acquisition->prepare(baseElemList[0]);
	acquisition->execute(0, acqData);
	const std::size_t signalLength = acqData.n2() * upsamplingFactor; // when using cross-correlation, the signal is longer than this
	if (deadZoneSamplesUp >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp (" << deadZoneSamplesUp <<
							") >= signalLength (" << signalLength << ").");
	}

	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[RANGE] ### baseElement: " << baseElement;

		acquisition->prepare(baseElement);

		for (unsigned int groupElemIdx = 0; groupElemIdx < config.numElements; ++groupElemIdx) {

			acquisition->execute(groupElemIdx, acqData);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PULSE_ECHO_USE_CROSS_CORRELATION
			tempSignal.resize(acqData.n2() * upsamplingFactor);
			interp.interpolate(&acqData(groupElemIdx, 0), acqData.n2(), &tempSignal[0]); // rx elem = tx elem

			Util::removeDC(&tempSignal[0], tempSignal.size(), deadZoneSamplesUp);

			// Cross-correlation using convolution.
			revRefPulseFilter.filter(tempSignal, signal);
#else
			signal.resize(acqData.n2() * upsamplingFactor);
			interp.interpolate(&acqData(groupElemIdx, 0), acqData.n2(), &signal[0]); // rx elem = tx elem

			Util::removeDC(&signal[0], signal.size(), deadZoneSamplesUp);

			envelope.calculate(&signal[0], signal.size());
#endif

			TFloat maxValue = 0;
			unsigned int idxMax = 0;
			for (unsigned int i = deadZoneSamplesUp; i < signal.size(); ++i) {
				if (signal[i] > maxValue) {
					maxValue = signal[i];
					idxMax = i;
				}
			}

			const TFloat tPulseEcho = (idxMax > signalOffset) ? ((idxMax - signalOffset) / (config.samplingFrequency * upsamplingFactor)) : 0;
			distanceList[groupElemIdx] = (config.propagationSpeed * tPulseEcho) / 2;

			LOG_DEBUG << "[RANGE] groupElemIdx: " << groupElemIdx << " tPulseEcho: " << tPulseEcho << " distance: " << distanceList[baseElemIdx];
		}

		const unsigned int pointPositionListOffset = baseElemIdx * (config.numElements - 1);
		for (unsigned int groupElemIdx = 0; groupElemIdx < config.numElements - 1; ++groupElemIdx) {
			const TFloat xElem = config.pitch * (baseElement + groupElemIdx);
			const TFloat sinPhi = (distanceList[groupElemIdx] - distanceList[groupElemIdx + 1]) / config.pitch;
			pointPositionList[pointPositionListOffset + groupElemIdx].first = xElem + distanceList[groupElemIdx] * sinPhi;
			const TFloat cosPhi = std::sqrt(1 - sinPhi * sinPhi);
			pointPositionList[pointPositionListOffset + groupElemIdx].second = distanceList[groupElemIdx] * cosPhi;
		}
	}

	for (unsigned int i = 0; i < pointPositionList.size(); ++i) {
		auxPointX[i] = pointPositionList[i].first;
		auxPointZ[i] = pointPositionList[i].second;
	}
	project_.showFigure2D(1, "Points", auxPointX, auxPointZ, true, true);

	LOG_DEBUG << "Saving the positions of the points...";
	project_.saveHDF5(pointPositionList, outputDir + "/point_pos_xz", "pointPositionList");
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPointsUsingCCBFPitchCatch()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"  , 1, 128);
	const auto numBaseElemSteps   = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);
	const auto firstRxElem        = methodPM->value<unsigned int>("first_rx_element"   , 0, config.numElements - 1);
	const auto lastRxElem         = methodPM->value<unsigned int>("last_rx_element"    , firstRxElem + 2, config.numElements - 1);
	const auto txElem             = methodPM->value<unsigned int>("tx_element"         , firstRxElem + 1, lastRxElem - 1);
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PITCH_CATCH_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	constauto peakOffset          = methodPM->value<TFloat>(      "peak_offset");
#endif

	project_.createDirectory(outputDir, false);

	//=====================================================================
	// Initialization.

	const unsigned int numActiveRxElements = lastRxElem - firstRxElem + 1;

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PITCH_CATCH_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;

	// This acquisition is done to obtain the signal length.
	acquisition->prepare(baseElemList[0]);
	acquisition->execute(txElem, acqData);
	const unsigned int origSignalLength = acqData.n2();
	const unsigned int requestedSignalLength = static_cast<unsigned int>(config.acquisitionTime * config.samplingFrequency + 0.5f); // rounds
	const unsigned int signalLength = std::min(origSignalLength, requestedSignalLength);
	const unsigned int signalStartOffset = static_cast<unsigned int>(config.samplingFrequency * (config.deadZoneM * 2.0f) / config.propagationSpeed + 0.5f); // rounds
	if (signalStartOffset >= signalLength - 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: signalStartOffset (" << signalStartOffset <<
							") >= signalLength (" << signalLength << ").");
	}
	LOG_DEBUG << "origSignalLength: " << origSignalLength << " requestedSignalLength: " << requestedSignalLength <<
			" signalLength: " << signalLength << " signalStartOffset: " << signalStartOffset;

	std::vector<std::pair<TFloat, TFloat>> pointPositionList; // x, z
	pointPositionList.reserve(baseElemList.size() * (numActiveRxElements - 1));

	// Store the acquisition data in memory.
	Tensor3<TFloat> acqDataList(numBaseElemSteps, numActiveRxElements, signalLength - signalStartOffset);
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		acquisition->prepare(baseElement);
		acquisition->execute(txElem, acqData);
		for (unsigned int rxElem = firstRxElem; rxElem <= lastRxElem; ++rxElem) {
			auto srcRange = acqData.range2(rxElem);
			auto destRange = acqDataList.range3(baseElemIdx, rxElem - firstRxElem);
			std::copy(
				srcRange.begin() + signalStartOffset,
				srcRange.end() - (origSignalLength - signalLength),
				destRange.begin());
		}
	}

	//=====================================================================
	// Processing.

	auto p = std::make_unique<CCBFPitchCatchProcessor<TFloat>>(
			config,
			acqDataList,
			upsamplingFactor,
			txElem,
			firstRxElem,
			lastRxElem,
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_PITCH_CATCH_USE_CROSS_CORRELATION
			refPulse,
#else
			peakOffset,
#endif
			signalStartOffset);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tProcess;
	for (unsigned int n = 0; n < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++n) {
	if (n <= 1U) { // n = 0: initial reset, n = 1: ignores the first iteration
		tProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tPartialPrepareData.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * baseElemList.size());
		p->tPartialProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * baseElemList.size());
	}
	Timer procTimer;
	pointPositionList.clear();
#endif

	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[POINT] ---------- baseElement = " << baseElement;

		p->process(baseElemIdx, baseElement, pointPositionList);
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcess.put(procTimer.getTime());
	}

	EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:", p->tPartialPrepareData, baseElemList.size());
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tProcess:    ", p->tPartialProcess, baseElemList.size());
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "process:     ", tProcess);
#endif

	std::vector<TFloat> figPointX(pointPositionList.size()), figPointZ(pointPositionList.size());
	for (unsigned int i = 0; i < pointPositionList.size(); ++i) {
		figPointX[i] = pointPositionList[i].first;
		figPointZ[i] = pointPositionList[i].second;
	}
	project_.showFigure2D(0, "Points", figPointX, figPointZ, true, true);

	std::vector<unsigned int> rxElemLimits = {firstRxElem, lastRxElem};
	project_.saveHDF5(rxElemLimits, outputDir + "/rx_elem_limits", "limits");

	LOG_DEBUG << "Saving the positions of the points...";
	project_.saveHDF5(pointPositionList, outputDir + "/point_pos_xz", "pointPositionList");
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPointsInArcs()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir          = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir            = methodPM->value<std::string>( "output_dir");
	const auto upsamplingFactor     = methodPM->value<unsigned int>("upsampling_factor"     , 1, 128);
	const auto highUpsamplingFactor = methodPM->value<unsigned int>("high_upsampling_factor", upsamplingFactor, 128);
	const auto numBaseElemSteps     = methodPM->value<unsigned int>("num_base_elem_steps"   , 3, config.numElementsMux - config.numElements + 1);
	const auto firstRxElem          = methodPM->value<unsigned int>("first_rx_element"      , 0, config.numElements - 1);
	const auto lastRxElem           = methodPM->value<unsigned int>("last_rx_element"       , firstRxElem + 2, config.numElements - 1);
	const auto txElem               = methodPM->value<unsigned int>("tx_element"            , firstRxElem + 1, lastRxElem - 1);
	const auto arcStepDiv           = methodPM->value<TFloat>(      "arc_step_div"          , 0.25, 32.0);
	const auto minArcAngle          = methodPM->value<TFloat>(      "min_arc_angle"         , 1.0, 179.0);
	const auto maxArcAngle          = methodPM->value<TFloat>(      "max_arc_angle"         , minArcAngle + 1.0, 179.0);
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_ARCS_USE_CROSS_CORRELATION
	const auto referencePulseFile   = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset           = methodPM->value<TFloat>(      "peak_offset");
#endif

	project_.createDirectory(outputDir, false);
	const unsigned int numActiveRxElements = lastRxElem - firstRxElem + 1;

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_ARC_USE_ANALYTIC_SIGNAL
	AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
#else
	CoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
#endif

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_ARCS_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;

	// This acquisition is done to obtain the signal length.
	acquisition->prepare(baseElemList[0]);
	acquisition->execute(txElem, acqData);
	const unsigned int origSignalLength = acqData.n2();
	const unsigned int requestedSignalLength = static_cast<unsigned int>(config.acquisitionTime * config.samplingFrequency + 0.5f); // rounds
	const unsigned int signalLength = std::min(origSignalLength, requestedSignalLength);
	const unsigned int signalStartOffset = static_cast<unsigned int>(config.samplingFrequency * (config.deadZoneM * 2.0f) / config.propagationSpeed + 0.5f); // rounds
	if (signalStartOffset >= signalLength - 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: signalStartOffset (" << signalStartOffset <<
							") >= signalLength (" << signalLength << ").");
	}
	LOG_DEBUG << "origSignalLength: " << origSignalLength << " requestedSignalLength: " << requestedSignalLength <<
			" signalLength: " << signalLength << " signalStartOffset: " << signalStartOffset;

	// Store the acquisition data in memory.
	Tensor3<TFloat> acqDataList(numBaseElemSteps, numActiveRxElements, signalLength - signalStartOffset);
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		acquisition->prepare(baseElement);
		acquisition->execute(txElem, acqData);
		for (unsigned int rxElem = firstRxElem; rxElem <= lastRxElem; ++rxElem) {
			auto srcRange = acqData.range2(rxElem);
			auto destRange = acqDataList.range3(baseElemIdx, rxElem - firstRxElem);
			std::copy(
				srcRange.begin() + signalStartOffset,
				srcRange.end() - (origSignalLength - signalLength),
				destRange.begin());
		}
	}

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_ARCS_USE_CROSS_CORRELATION
	auto p = std::make_unique<CrossCorrArcCylindricalWaveProcessor<TFloat>>(
			config,
			acqDataList,
			upsamplingFactor,
			highUpsamplingFactor,
			txElem,
			firstRxElem,
			lastRxElem,
			arcStepDiv,
			minArcAngle,
			maxArcAngle,
			coherenceFactor,
			refPulse,
			signalStartOffset);
#else
	auto p = std::make_unique<ArcCylindricalWaveProcessor<TFloat>>(
			config,
			acqDataList,
			upsamplingFactor,
			highUpsamplingFactor,
			txElem,
			firstRxElem,
			lastRxElem,
			arcStepDiv,
			minArcAngle,
			maxArcAngle,
			coherenceFactor,
			peakOffset,
			signalStartOffset);
#endif

	std::vector<XZValueFactor<TFloat>> arcData;
	std::vector<std::pair<TFloat, TFloat>> pointPositionList(baseElemList.size()); // x, z

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_ARCS
	std::vector<TFloat> figArcX, figArcY, angle;
	enum {
		PEAK_VALUE,
		PEAK_ANGLE,
		MAIN_LOBE_WIDTH,
		MAX_SIDE_LOBE_LEVEL,
		STAT_SIZE
	};
	Matrix<TFloat> statistics(baseElemList.size(), STAT_SIZE);
#endif

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tProcess;
	for (unsigned int n = 0; n < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++n) {
	if (n <= 1U) { // n = 0: initial reset, n = 1: ignores the first iteration
		tProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tPartialPrepareData.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * baseElemList.size());
		p->tPartialProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * baseElemList.size());
	}
	Timer procTimer;
#endif

	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[POINT] ---------- baseElement = " << baseElement;

		// Obtain the arc.
		p->process(baseElemIdx, arcData);

		if (coherenceFactor.enabled()) {
			// Apply the CF.
			for (auto& point : arcData) {
				point.value *= point.factor;
				point.factor = 1;
			}
		}

		TFloat maxValue = 0;
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_ARC_USE_LOBE_CENTER
		// Find the main lobe center.
		for (unsigned int i = 0, size = arcData.size(); i < size; ++i) {
			if (arcData[i].value > maxValue) {
				maxValue = arcData[i].value;
			}
		}
		const TFloat threshold = CYLINDER_DETECTION_AND_FERMAT_METHOD_ARC_LOBE_THRESHOLD * maxValue;
		unsigned int idx1 = 0;
		for (unsigned int i = 0, size = arcData.size(); i < size; ++i) {
			if (arcData[i].value > threshold) {
				idx1 = i;
				break;
			}
		}

		unsigned int idx2 = arcData.size() - 1;
		for (int i = idx2; i >= 0; --i) {
			if (arcData[i].value > threshold) {
				idx2 = i;
				break;
			}
		}
		const unsigned int idxLobe = static_cast<unsigned int>(0.5 * (idx1 + idx2) + 0.5);
#else
		// Find the peak.
		unsigned int idxLobe = 0;
		for (unsigned int i = 0; i < arcData.size(); ++i) {
			if (arcData[i].value > maxValue) {
				maxValue = arcData[i].value;
				idxLobe = i;
			}
		}
#endif

		// Converts to the coordinate system centered on the first element.
		pointPositionList[baseElemIdx].first = arcData[idxLobe].x + (baseElement + (config.numElements - 1) / TFloat(2)) * config.pitch;
		pointPositionList[baseElemIdx].second = arcData[idxLobe].z;

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_ARCS
#define CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_ARC_MAIN_LOBE_HALF_RANGE 15.0
		const TFloat xArcCenter = (-(config.numElements - 1.0) / 2.0 + txElem) * config.pitch;
		figArcX.resize(arcData.size());
		figArcY.resize(arcData.size());
		angle.resize(arcData.size());
		for (unsigned int i = 0; i < arcData.size(); ++i) {
			angle[i] = Util::radianToDegree(std::atan2(arcData[i].z, arcData[i].x - xArcCenter));

			figArcX[i] = angle[i];
			figArcY[i] = arcData[i].value;
		}
		std::ostringstream windowName;
		windowName << "Arc " << baseElemIdx;
		project_.showFigure2D(baseElemIdx + 1, windowName.str().c_str(), figArcX, figArcY, true);

		LOG_DEBUG << "##### Peak value: " << maxValue << " angle: " << angle[idxLobe];
		statistics(baseElemIdx, PEAK_VALUE) = maxValue;
		statistics(baseElemIdx, PEAK_ANGLE) = angle[idxLobe];

		const TFloat halfPeakValue = maxValue / 2;
		unsigned int minMainLobeIdx = 0;
		unsigned int maxMainLobeIdx = 0;
		TFloat maxSideLobeLevel = 0;
		for (unsigned int i = 0; i < arcData.size(); ++i) {
			if (std::abs(angle[i] - angle[idxLobe]) < TFloat(CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_ARC_MAIN_LOBE_HALF_RANGE)) {
				if (arcData[i].value >= halfPeakValue) {
					if (minMainLobeIdx == 0) {
						minMainLobeIdx = i;
					}
					maxMainLobeIdx = i;
				}
			} else {
				if (arcData[i].value > maxSideLobeLevel) {
					maxSideLobeLevel = arcData[i].value;
				}
			}
			angle[i] = Util::radianToDegree(std::atan2(arcData[i].z, arcData[i].x - xArcCenter));

			figArcX[i] = angle[i];
			figArcY[i] = arcData[i].value;
		}

		LOG_DEBUG << "##### Main lobe width (-6 dB): " << angle[maxMainLobeIdx] - angle[minMainLobeIdx];
		LOG_DEBUG << "##### Max side lobe level: " << maxSideLobeLevel;
		statistics(baseElemIdx, MAIN_LOBE_WIDTH) = angle[maxMainLobeIdx] - angle[minMainLobeIdx];
		statistics(baseElemIdx, MAX_SIDE_LOBE_LEVEL) = maxSideLobeLevel;

		project_.saveHDF5(angle, String::Begin() << outputDir << "/arc_angle_base" << baseElement << String::End(), "angle");
		project_.saveHDF5(arcData, String::Begin() << outputDir << "/arc_value_base" << baseElement << String::End(), "value", Util::CopyValueOp());
#endif
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcess.put(procTimer.getTime());
	}

	EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:", p->tPartialPrepareData, baseElemList.size());
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tProcess:    ", p->tPartialProcess, baseElemList.size());
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "process:     ", tProcess);
#endif

#ifdef CYL_DETECT_AND_FERMAT_METHOD_POINT_DETECTION_SHOW_ARCS
	std::ostringstream out;
	out << "[";
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		if (baseElemIdx != 0) {
			out << ' ';
		}
		out << '[' << baseElement << ',';
		for (unsigned int i = 0; i < STAT_SIZE; ++i) {
			out << statistics(baseElemIdx, i);
			if (i != STAT_SIZE - 1) {
				out << ',';
			}
		}
		out << ']';
		if (baseElemIdx != baseElemList.size() - 1U) {
			out << ",\n";
		}
	}
	out << "]\n";
	LOG_DEBUG << "##### Statistics:\n" << out.str();
#endif

	std::vector<TFloat> figPointX(pointPositionList.size()), figPointZ(pointPositionList.size());
	for (unsigned int i = 0; i < pointPositionList.size(); ++i) {
		figPointX[i] = pointPositionList[i].first;
		figPointZ[i] = pointPositionList[i].second;
	}
	project_.showFigure2D(0, "Points", figPointX, figPointZ, true, true);

	std::vector<unsigned int> rxElemLimits = {firstRxElem, lastRxElem};
	project_.saveHDF5(rxElemLimits, outputDir + "/rx_elem_limits", "limits");

	LOG_DEBUG << "Saving the positions of the points...";
	project_.saveHDF5(pointPositionList, outputDir + "/point_pos_xz", "pointPositionList");
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::detectPointsUsingTangentCurveGeometry()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"  , 1, 128);
	const auto numBaseElemSteps   = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);
	const auto firstRxElem        = methodPM->value<unsigned int>("first_rx_element"   , 0, config.numElements - 1);
	const auto lastRxElem         = methodPM->value<unsigned int>("last_rx_element"    , firstRxElem + 2, config.numElements - 1);
	const auto txElem             = methodPM->value<unsigned int>("tx_element"         , firstRxElem + 1, lastRxElem - 1);
	const auto arcStepDiv         = methodPM->value<TFloat>(      "arc_step_div"       , 0.25, 32.0);
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_TANGENT_CURVE_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset");
#endif

	project_.createDirectory(outputDir, false);

	//=====================================================================
	// Initialization.

	const unsigned int numActiveRxElements = lastRxElem - firstRxElem + 1;
	const TFloat lambda = config.propagationSpeed / config.centerFrequency;
	const TFloat arcStep = lambda / arcStepDiv;

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_TANGENT_CURVE_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;

	// This acquisition is done to obtain the signal length.
	acquisition->prepare(baseElemList[0]);
	acquisition->execute(txElem, acqData);
	const unsigned int origSignalLength = acqData.n2();
	const unsigned int requestedSignalLength = static_cast<unsigned int>(config.acquisitionTime * config.samplingFrequency + 0.5f); // rounds
	const unsigned int signalLength = std::min(origSignalLength, requestedSignalLength);
	const unsigned int signalStartOffset = static_cast<unsigned int>(config.samplingFrequency * (config.deadZoneM * 2.0f) / config.propagationSpeed + 0.5f); // rounds
	if (signalStartOffset >= signalLength - 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: signalStartOffset (" << signalStartOffset <<
							") >= signalLength (" << signalLength << ").");
	}
	LOG_DEBUG << "origSignalLength: " << origSignalLength << " requestedSignalLength: " << requestedSignalLength <<
			" signalLength: " << signalLength << " signalStartOffset: " << signalStartOffset;

	std::vector<std::pair<TFloat, TFloat>> pointPositionList; // x, z

	// Store the acquisition data in memory.
	Tensor3<TFloat> acqDataList(numBaseElemSteps, numActiveRxElements, signalLength - signalStartOffset);
	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		acquisition->prepare(baseElement);
		acquisition->execute(txElem, acqData);
		for (unsigned int rxElem = firstRxElem; rxElem <= lastRxElem; ++rxElem) {
			auto srcRange = acqData.range2(rxElem);
			auto destRange = acqDataList.range3(baseElemIdx, rxElem - firstRxElem);
			std::copy(
				srcRange.begin() + signalStartOffset,
				srcRange.end() - (origSignalLength - signalLength),
				destRange.begin());
		}
	}

	//=====================================================================
	// Processing.

	auto p = std::make_unique<TangentCurveGeometryProcessor<TFloat>>(
			config,
			acqDataList,
			upsamplingFactor,
			txElem,
			firstRxElem,
			lastRxElem,
			arcStep,
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_TANGENT_CURVE_USE_CROSS_CORRELATION
			refPulse,
#else
			peakOffset,
#endif
			signalStartOffset);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tProcess;
	for (unsigned int n = 0; n < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++n) {
	if (n <= 1U) { // n = 0: initial reset, n = 1: ignores the first iteration
		tProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tPartialPrepareData.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * baseElemList.size());
		p->tPartialProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * baseElemList.size());
	}
	Timer procTimer;
	pointPositionList.clear();
#endif

	for (unsigned int baseElemIdx = 0; baseElemIdx < baseElemList.size(); ++baseElemIdx) {
		const unsigned int baseElement = baseElemList[baseElemIdx];
		LOG_DEBUG << "[POINT] ---------- baseElement = " << baseElement;

		p->process(baseElemIdx, baseElement, pointPositionList);
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcess.put(procTimer.getTime());
	}

	EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:", p->tPartialPrepareData, baseElemList.size());
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tProcess:    ", p->tPartialProcess, baseElemList.size());
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "process:     ", tProcess);
#endif

	std::vector<TFloat> figPointX(pointPositionList.size()), figPointZ(pointPositionList.size());
	for (unsigned int i = 0; i < pointPositionList.size(); ++i) {
		figPointX[i] = pointPositionList[i].first;
		figPointZ[i] = pointPositionList[i].second;
	}
	project_.showFigure2D(0, "Points", figPointX, figPointZ, true, true);

	std::vector<unsigned int> rxElemLimits = {firstRxElem, lastRxElem};
	project_.saveHDF5(rxElemLimits, outputDir + "/rx_elem_limits", "limits");

	LOG_DEBUG << "Saving the positions of the points...";
	project_.saveHDF5(pointPositionList, outputDir + "/point_pos_xz", "pointPositionList");
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::fitCircle()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const TwoMediumSTAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto outputDir        = methodPM->value<std::string>( "output_dir");
	const auto numBaseElemSteps = methodPM->value<unsigned int>("num_base_elem_steps", 3, config.numElementsMux - config.numElements + 1);

	project_.createDirectory(outputDir, false);

	std::vector<std::pair<TFloat, TFloat>> pointPositionList;
	project_.loadHDF5(outputDir + "/point_pos_xz", "pointPositionList", pointPositionList);

	TFloat centerInterfaceX, centerInterfaceZ, interfaceR;
	unsigned int numPointsWithNaN;

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tProcess;
	for (unsigned int n = 0; n < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++n) {
	if (n <= 1U) { // n = 0: initial reset, n = 1: ignores the first iteration
		tProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	Timer procTimer;
#endif

	// Execute a circle fitting.
	if (numBaseElemSteps > pointPositionList.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The parameter num_base_elem_steps must be <= " << pointPositionList.size() << '.');
	} else if (numBaseElemSteps == pointPositionList.size()) {
		Fitting::circleFittingByPratt(pointPositionList, centerInterfaceX, centerInterfaceZ, interfaceR, &numPointsWithNaN);
	} else {
		if (numBaseElemSteps < 9) { //TODO: Change this number if more groups are needed
			std::vector<std::pair<TFloat, TFloat>> pointPositionSubList(numBaseElemSteps);
			double step = (pointPositionList.size() - 1.0) / (numBaseElemSteps - 1.0);
			for (unsigned int i = 0; i < numBaseElemSteps; ++i) {
				unsigned int srcIndex = (i * step) + 0.5; // round
				pointPositionSubList[i] = pointPositionList[srcIndex];
				LOG_DEBUG << "pointPositionSubList: srcIndex=" << srcIndex << " i=" << i <<
					" x=" << pointPositionSubList[i].first << " z=" << pointPositionSubList[i].second;
			}
			//TODO: if numBaseElemSteps == 3, use another algorithm. With this function a "negative error" occurs.
			Fitting::circleFittingByPratt(pointPositionSubList, centerInterfaceX, centerInterfaceZ, interfaceR, &numPointsWithNaN);
		} else {
			// For CCBF pulse-echo.
			Fitting::circleFittingByPratt(pointPositionList, centerInterfaceX, centerInterfaceZ, interfaceR, &numPointsWithNaN);
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcess.put(procTimer.getTime());
	}

	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("process:", tProcess);
#endif

	LOG_DEBUG << "centerInterfaceX = " << centerInterfaceX << " centerInterfaceZ = " << centerInterfaceZ << " interfaceR = " << interfaceR;

	{
		std::vector<unsigned int> v = {numPointsWithNaN};
		project_.saveHDF5(v, outputDir + "/circle_fitting_num_points_with_nan", "n");

		std::vector<TFloat> xzr{centerInterfaceX, centerInterfaceZ, interfaceR};
		LOG_DEBUG << "Saving the circle center and radius...";
		project_.saveHDF5(xzr, outputDir + "/circle_xzr", "xzr");
	}
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::execTwoMediumImaging()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const TwoMediumSTAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto imageDir           = methodPM->value<std::string>( "image_dir");
	const auto finalImageBaseElem = methodPM->value<unsigned int>("final_image_base_element",         0, config.numElementsMux - config.numElements);
	const auto normalizationMinZ  = methodPM->value<TFloat>(      "normalization_min_z"     , -500.0e-3, 500.0e-3);
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"       ,         1,      128);
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset"             ,       0.0,     50.0);
	const auto maxFermatBlockSize = methodPM->value<TFloat>(      "max_fermat_block_size"   ,    1.0e-3,    100.0);
	const auto interfaceStepDiv   = methodPM->value<TFloat>(      "interface_step_div"      ,      0.25,     32.0);
	const auto interfaceAngle1    = methodPM->value<TFloat>(      "interface_angle_1"       ,      -170,    -10.0);
	const auto interfaceAngle2    = methodPM->value<TFloat>(      "interface_angle_2"       ,      -170,    -10.0);

	project_.createDirectory(imageDir, false);
	Matrix<XZValueFactor<TFloat>> gridData;

	const TFloat lambda2 = config.propagationSpeed2 / config.centerFrequency;
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), lambda2, gridData);

	std::vector<TFloat> xzr;
	project_.loadHDF5(outputDir + "/circle_xzr", "xzr", xzr);
	const TFloat centerInterfaceX = xzr[0];
	const TFloat centerInterfaceZ = xzr[1];
	const TFloat interfaceR = xzr[2];

	std::vector<TFloat> interfaceAngleList;
	Util::fillSequenceFromStartToEndWithMaximumStep(
		interfaceAngleList,
		Util::degreeToRadian(interfaceAngle1),
		Util::degreeToRadian(interfaceAngle2),
		lambda2 / (interfaceStepDiv * interfaceR));
	LOG_DEBUG << "interfaceAngleList.size() = " << interfaceAngleList.size();

	// Generate points in the interface.
	std::vector<std::pair<TFloat, TFloat>> interfacePointList(interfaceAngleList.size()); // x, z
	for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
		interfacePointList[i] = std::make_pair(
				centerInterfaceX + interfaceR * std::cos(interfaceAngleList[i]),
				centerInterfaceZ + interfaceR * std::sin(interfaceAngleList[i]));
	}
	LOG_DEBUG << "Saving the interface points...";
	project_.saveHDF5(interfacePointList, imageDir + "/interface_xz", "interfacePointList");

	// Convert the interface points to the element group local coordinate system.
	std::vector<XZ<TFloat>> localInterfacePointList(interfacePointList.size()); // x, z
	const TFloat offset = (finalImageBaseElem + (config.numElements - 1) / 2.0) * config.pitch;
	for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
		localInterfacePointList[i].x = interfacePointList[i].first - offset;
		localInterfacePointList[i].z = interfacePointList[i].second;
	}

	bool coherenceFactorEnabled = false;

	// Obtain the final image.
	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(
					project_,
					config.numElements,
					savedAcqDir + '/');

	if (project_.method() == MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_vectorial_dp ||
			project_.method() == MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_vectorial_sp) {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
		coherenceFactorEnabled = coherenceFactor.enabled();
		VectorialTwoMediumSTAProcessor<TFloat> p(config, *acquisition, upsamplingFactor, coherenceFactor, maxFermatBlockSize, peakOffset);
		p.prepare(finalImageBaseElem);
		p.process(localInterfacePointList, gridData);
	} else {
		CoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
		coherenceFactorEnabled = coherenceFactor.enabled();
		TwoMediumSTAProcessor<TFloat> p(config, *acquisition, coherenceFactor, peakOffset);
		p.prepare(finalImageBaseElem);
		p.process(localInterfacePointList, gridData);
	}

	TFloat maxValue = 0;
	for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
		if (iter->z >= normalizationMinZ && iter->value > maxValue) {
			maxValue = iter->value;
		}
	}

	// Normalization with saturation.
	const TFloat normCoef = 1 / maxValue;
	for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
		if (iter->value > maxValue) { // saturation
			iter->value = 1;
		} else {
			iter->value *= normCoef;
		}
	}

#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
	LOG_DEBUG << "Saving the final raw image... baseElement = " << finalImageBaseElem;
	project_.saveHDF5(gridData, String::Begin() << imageDir << "/final_raw_image_" << finalImageBaseElem << String::End(), "image", Util::CopyValueOp());
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, imageDir + "/final_image_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, imageDir + "/final_image_z", "z", Util::CopyZOp());
#endif

	project_.showFigure3D(1, "Raw", &gridData, &localInterfacePointList,
				true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS);

	if (coherenceFactorEnabled) {
#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
		LOG_DEBUG << "Saving the CF factors...";
		project_.saveHDF5(gridData, String::Begin() << imageDir << "/final_cf_factors_" << finalImageBaseElem << String::End(), "factors", Util::CopyFactorOp());
#endif
		// Apply the CF.
		for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->value *= iter->factor;
			iter->factor = 1.0;
		}
#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
		LOG_DEBUG << "Saving the final image (after CF)... baseElement = " << finalImageBaseElem;
		project_.saveHDF5(gridData, String::Begin() << imageDir << "/final_cf_image_" << finalImageBaseElem << String::End(), "image", Util::CopyValueOp());
#endif
		project_.showFigure3D(2, "CF", &gridData, &localInterfacePointList,
					true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS);
	}
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::execCombinedTwoMediumImaging()
{
#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
	const auto imageDir = project_.taskParamMap().value<std::string>("image_dir");
	project_.createDirectory(imageDir, false);
#endif
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const TwoMediumSTAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto rxApodFile         = methodPM->value<std::string>( "rx_apod_file");
	const auto normalizationMinZ  = methodPM->value<TFloat>(      "normalization_min_z"  , -500.0e-3, 500.0e-3);

	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"    ,         1,      128);
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset"          ,       0.0,     50.0);
	const auto maxFermatBlockSize = methodPM->value<TFloat>(      "max_fermat_block_size",    1.0e-3,    100.0);
	const auto numBaseElemSteps   = methodPM->value<unsigned int>("num_base_elem_steps"  ,         3, config.numElementsMux - config.numElements + 1);
	const auto interfaceStepDiv   = methodPM->value<TFloat>(      "interface_step_div"   ,      0.25,     32.0);
	const auto interfaceAngle1    = methodPM->value<TFloat>(      "interface_angle_1"    ,      -170,    -10.0);
	const auto interfaceAngle2    = methodPM->value<TFloat>(      "interface_angle_2"    ,      -170,    -10.0);
	const auto txElem             = methodPM->value<unsigned int>("tx_element"           ,         0, config.numElements - 1);

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	Matrix<XZComplexValue<TFloat>> gridData;
	Matrix<std::complex<TFloat>> gridValue;
#else
	Matrix<XZValue<TFloat>> gridData;
	Matrix<TFloat> gridValue;
#endif
	Matrix<XZ<TFloat>> gridXZ;

	const TFloat lambda2 = config.propagationSpeed2 / config.centerFrequency;
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), lambda2, gridXZ);
	gridValue.resize(gridXZ.n1(), gridXZ.n2());
	gridData.resize(gridXZ.n1(), gridXZ.n2());

	Value::copyXZSequence(gridXZ.begin(), gridXZ.end(), gridData.begin());

	std::vector<TFloat> rxApod;
	project_.loadHDF5(rxApodFile, "apod", rxApod);

	std::vector<TFloat> xzr;
	project_.loadHDF5(outputDir + "/circle_xzr", "xzr", xzr);
	const TFloat centerInterfaceX = xzr[0];
	const TFloat centerInterfaceZ = xzr[1];
	const TFloat interfaceR = xzr[2];

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
#else
	CoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
#endif

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::vector<typename VectorialCombinedTwoMediumImagingProcessor<TFloat>::StepConfiguration> stepConfigList;
	unsigned int numActiveTxElements = 0;
	switch (project_.method()) {
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_dp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_sp:
		for (unsigned int i = 0; i < baseElemList.size(); ++i) {
			typename VectorialCombinedTwoMediumImagingProcessor<TFloat>::StepConfiguration conf{
				i,
				baseElemList[i],
				txElem, txElem
			};
			stepConfigList.push_back(conf);
			LOG_DEBUG << "Step config: " << conf.baseElem << ' ' << conf.firstTxElem << ' ' << conf.lastTxElem;
		}
		numActiveTxElements = 1;
		break;
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_sta_dp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_sta_sp:
		for (unsigned int i = 0; i < baseElemList.size(); ++i) {
			typename VectorialCombinedTwoMediumImagingProcessor<TFloat>::StepConfiguration conf{
				i,
				baseElemList[i],
				0, config.numElements - 1
			};
			stepConfigList.push_back(conf);
		}
		numActiveTxElements = config.numElements;
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;

	// This acquisition is done to obtain the signal length.
	acquisition->prepare(baseElemList[0]);
	acquisition->execute(txElem, acqData);
	const unsigned int origSignalLength = acqData.n2();
	const unsigned int requestedSignalLength = static_cast<unsigned int>(config.acquisitionTime * config.samplingFrequency + 0.5f); // rounds
	const unsigned int signalLength = std::min(origSignalLength, requestedSignalLength);
	const unsigned int signalStartOffset = static_cast<unsigned int>(config.samplingFrequency * (config.deadZoneM * 2.0f) / config.propagationSpeed1 + 0.5f); // rounds
	if (signalStartOffset >= signalLength - 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: signalStartOffset (" << signalStartOffset <<
							") >= signalLength (" << signalLength << ").");
	}
	LOG_DEBUG << "origSignalLength: " << origSignalLength << " requestedSignalLength: " << requestedSignalLength <<
			" signalLength: " << signalLength << " signalStartOffset: " << signalStartOffset;

	// Store the acquisition data in memory.
	std::vector<Tensor3<TFloat>> acqDataList(numBaseElemSteps);
	for (auto& d : acqDataList) {
		d.resize(numActiveTxElements, config.numElements, signalLength - signalStartOffset);
	}
	for (const auto& stepConfig : stepConfigList) {
		acquisition->prepare(stepConfig.baseElem);
		for (unsigned int txElem = stepConfig.firstTxElem; txElem <= stepConfig.lastTxElem; ++txElem) {
			acquisition->execute(txElem, acqData);
			for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
				auto srcRange = acqData.range2(rxElem);
				auto destRange = acqDataList[stepConfig.baseElemIdx].range3(txElem - stepConfig.firstTxElem, rxElem);
				std::copy(
					srcRange.begin() + signalStartOffset,
					srcRange.end() - (origSignalLength - signalLength),
					destRange.begin());
			}
		}
	}

	auto p = std::make_unique<VectorialCombinedTwoMediumImagingProcessor<TFloat>>(
			config,
			acqDataList,
			upsamplingFactor,
			coherenceFactor,
			maxFermatBlockSize,
			peakOffset,
			signalStartOffset);

	std::vector<TFloat> interfaceAngleList;
	std::vector<XZ<TFloat>> interfacePointList;

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tProcess;
	for (unsigned int n = 0; n < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++n) {
	if (n <= 1U) { // n = 0: initial reset, n = 1: ignores the first iteration
		tProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tMinRowIdx.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tMedium1DelayMatrix.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tCalculateDelays.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tPrepareData.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tProcessColumn.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	Timer procTimer;
#endif

	Util::fillSequenceFromStartToEndWithMaximumStep(
		interfaceAngleList,
		Util::degreeToRadian(interfaceAngle1),
		Util::degreeToRadian(interfaceAngle2),
		lambda2 / (interfaceStepDiv * interfaceR));
	LOG_DEBUG << "interfaceAngleList.size() = " << interfaceAngleList.size();

	// Generate points in the interface.
	interfacePointList.resize(interfaceAngleList.size()); // x, z
	const TFloat centerOffset = ((config.numElementsMux - 1) / 2.0) * config.pitch;
	for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
		interfacePointList[i] = XZ<TFloat>(
				centerInterfaceX + interfaceR * std::cos(interfaceAngleList[i]) - centerOffset,
				centerInterfaceZ + interfaceR * std::sin(interfaceAngleList[i]));
	}

	p->process(stepConfigList, interfacePointList, rxApod, gridXZ, gridValue);

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	Value::transformSequence(gridValue.begin(), gridValue.end(), gridData.begin(), Value::ScalarToValueFieldOp<std::complex<TFloat>, XZComplexValue<TFloat>>());
#else
	Value::transformSequence(gridValue.begin(), gridValue.end(), gridData.begin(), Value::ScalarToValueFieldOp<TFloat, XZValue<TFloat>>());
#endif

	TFloat maxValue = 0;
	for (auto& point : gridData) {
		const TFloat absValue = std::abs(point.value);
		if (point.z >= normalizationMinZ && absValue > maxValue) {
			maxValue = absValue;
		}
	}

	if (maxValue > 0) {
		// Normalization with saturation.
		const TFloat normCoef = 1 / maxValue;
		for (auto& point : gridData) {
			if (std::abs(point.value) > maxValue) { // saturation
				point.value = 1;
			} else {
				point.value *= normCoef;
			}
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcess.put(procTimer.getTime());
	}

	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMinRowIdx:         ", p->tMinRowIdx);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMedium1DelayMatrix:", p->tMedium1DelayMatrix);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tCalculateDelays:   ", p->tCalculateDelays);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tPrepareData:       ", p->tPrepareData);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tProcessColumn:     ", p->tProcessColumn);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("process:            ", tProcess);

	std::ostringstream out;
	out << std::setprecision(15) << "t = [";
	for (auto v: tProcess.list()) {
		out << v << ", ";
	}
	out << "]";
	LOG_INFO << out.str();
#endif

#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
	LOG_DEBUG << "Saving the final raw image...";
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	project_.saveHDF5(gridData, imageDir + "/final_raw_image", "image", Util::CopyAbsValueOp());
# else
	project_.saveHDF5(gridData, imageDir + "/final_raw_image", "image", Util::CopyValueOp());
# endif
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, imageDir + "/final_image_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, imageDir + "/final_image_z", "z", Util::CopyZOp());
#endif

	project_.showFigure3D(1, "Raw", &gridData, &interfacePointList,
				true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_INVERTED_GRAY);
}

template<typename TFloat>
template<typename Proc>
void
CylinderDetectionAndFermatMethod<TFloat>::execCombinedTwoMediumImagingCyl()
{
#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
	const auto imageDir = project_.taskParamMap().value<std::string>("image_dir");
	project_.createDirectory(imageDir, false);
#endif
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const TwoMediumSTAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto outputDir          = methodPM->value<std::string>( "output_dir");
	const auto rxApodFile         = methodPM->value<std::string>( "rx_apod_file");
	const auto normalizationMinZ  = methodPM->value<TFloat>(      "normalization_min_z"  , -500.0e-3, 500.0e-3);

	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor"    ,         1,      128);
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset"          ,       0.0,     50.0);
	const auto maxFermatBlockSize = methodPM->value<TFloat>(      "max_fermat_block_size",    1.0e-3,    100.0);
	const auto numBaseElemSteps   = methodPM->value<unsigned int>("num_base_elem_steps"  ,         3, config.numElementsMux - config.numElements + 1);
	const auto interfaceStepDiv   = methodPM->value<TFloat>(      "interface_step_div"   ,      0.25,     32.0);
	const auto interfaceAngle1    = methodPM->value<TFloat>(      "interface_angle_1"    ,      -170,    -10.0);
	const auto interfaceAngle2    = methodPM->value<TFloat>(      "interface_angle_2"    ,      -170,    -10.0);
	const auto txElem             = methodPM->value<unsigned int>("tx_element"           ,         0, config.numElements - 1);

	Matrix<XZComplexValue<TFloat>> gridData;
	Matrix<std::complex<TFloat>> gridValue;
	Matrix<XZ<TFloat>> gridXZ;

	const TFloat lambda2 = config.propagationSpeed2 / config.centerFrequency;
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), lambda2, gridXZ);
	gridValue.resize(gridXZ.n1(), gridXZ.n2());
	gridData.resize(gridXZ.n1(), gridXZ.n2());

	Value::copyXZSequence(gridXZ.begin(), gridXZ.end(), gridData.begin());

	std::vector<TFloat> rxApod;
	project_.loadHDF5(rxApodFile, "apod", rxApod);

	std::vector<TFloat> xzr;
	project_.loadHDF5(outputDir + "/circle_xzr", "xzr", xzr);
	const TFloat centerInterfaceX = xzr[0];
	const TFloat centerInterfaceZ = xzr[1];
	const TFloat interfaceR = xzr[2];

	AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));

	std::vector<unsigned int> baseElemList;
	Util::fillSequenceFromStartToEndWithSize(baseElemList, 0U, config.numElementsMux - config.numElements, numBaseElemSteps);

	std::vector<typename Proc::StepConfiguration> stepConfigList;
	for (unsigned int i = 0; i < baseElemList.size(); ++i) {
		typename Proc::StepConfiguration conf{i, baseElemList[i], txElem};
		stepConfigList.push_back(conf);
		LOG_DEBUG << "Step config: " << conf.baseElem << ' ' << conf.txElem;
	}

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;

	// This acquisition is done to obtain the signal length.
	acquisition->prepare(baseElemList[0]);
	acquisition->execute(txElem, acqData);
	const unsigned int origSignalLength = acqData.n2();
	const unsigned int requestedSignalLength = static_cast<unsigned int>(config.acquisitionTime * config.samplingFrequency + 0.5f); // rounds
	const unsigned int signalLength = std::min(origSignalLength, requestedSignalLength);
	const unsigned int signalStartOffset = static_cast<unsigned int>(config.samplingFrequency * (config.deadZoneM * 2.0f) / config.propagationSpeed1 + 0.5f); // rounds
	if (signalStartOffset >= signalLength - 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: signalStartOffset (" << signalStartOffset <<
							") >= signalLength (" << signalLength << ").");
	}
	LOG_DEBUG << "origSignalLength: " << origSignalLength << " requestedSignalLength: " << requestedSignalLength <<
			" signalLength: " << signalLength << " signalStartOffset: " << signalStartOffset;

	// Store the acquisition data in memory.
	std::vector<Matrix<TFloat>> acqDataList(numBaseElemSteps);
	for (auto& d : acqDataList) {
		d.resize(config.numElements, signalLength - signalStartOffset);
	}
	for (const auto& stepConfig : stepConfigList) {
		const unsigned int txElem = stepConfig.txElem;
		acquisition->prepare(stepConfig.baseElem);
		acquisition->execute(txElem, acqData);
		for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
			auto srcRange = acqData.range2(rxElem);
			auto destRange = acqDataList[stepConfig.baseElemIdx].range2(rxElem);
			std::copy(
				srcRange.begin() + signalStartOffset,
				srcRange.end() - (origSignalLength - signalLength),
				destRange.begin());
		}
	}

	auto p = std::make_unique<Proc>(
			config,
			acqDataList,
			upsamplingFactor,
			coherenceFactor,
			maxFermatBlockSize,
			peakOffset,
			signalStartOffset);

	std::vector<TFloat> interfaceAngleList;
	std::vector<XZ<TFloat>> interfacePointList;

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tProcess;
	for (unsigned int n = 0; n < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++n) {
	if (n <= 1U) { // n = 0: initial reset, n = 1: ignores the first iteration
		tProcess.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tMinRowIdx.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tMedium1DelayMatrix.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tCalculateDelays.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tPrepareData.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		p->tProcessColumn.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	Timer procTimer;
#endif

	Util::fillSequenceFromStartToEndWithMaximumStep(
		interfaceAngleList,
		Util::degreeToRadian(interfaceAngle1),
		Util::degreeToRadian(interfaceAngle2),
		lambda2 / (interfaceStepDiv * interfaceR));
	LOG_DEBUG << "interfaceAngleList.size() = " << interfaceAngleList.size();

	// Generate points in the interface.
	interfacePointList.resize(interfaceAngleList.size()); // x, z
	const TFloat centerOffset = ((config.numElementsMux - 1) / 2.0) * config.pitch;
	for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
		interfacePointList[i] = XZ<TFloat>(
				centerInterfaceX + interfaceR * std::cos(interfaceAngleList[i]) - centerOffset,
				centerInterfaceZ + interfaceR * std::sin(interfaceAngleList[i]));
	}

	p->process(stepConfigList, interfacePointList, rxApod, gridXZ, gridValue);

	Value::transformSequence(gridValue.begin(), gridValue.end(), gridData.begin(), Value::ScalarToValueFieldOp<std::complex<TFloat>, XZComplexValue<TFloat>>());

	TFloat maxValue = 0;
	for (auto& point : gridData) {
		const TFloat absValue = std::abs(point.value);
		if (point.z >= normalizationMinZ && absValue > maxValue) {
			maxValue = absValue;
		}
	}

	if (maxValue > 0) {
		// Normalization with saturation.
		const TFloat normCoef = 1 / maxValue;
		for (auto& point : gridData) {
			if (std::abs(point.value) > maxValue) { // saturation
				point.value = 1;
			} else {
				point.value *= normCoef;
			}
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcess.put(procTimer.getTime());
	}

	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMinRowIdx:         ", p->tMinRowIdx);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMedium1DelayMatrix:", p->tMedium1DelayMatrix);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tCalculateDelays:   ", p->tCalculateDelays);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tPrepareData:       ", p->tPrepareData);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tProcessColumn:     ", p->tProcessColumn);
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("process:            ", tProcess);

	std::ostringstream out;
	out << std::setprecision(15) << "t = [";
	for (auto v: tProcess.list()) {
		out << v << ", ";
	}
	out << "]";
	LOG_INFO << out.str();
#endif

#ifdef CYL_DETECT_AND_FERMAT_METHOD_IMAGING_SAVE_DATA
	LOG_DEBUG << "Saving the final raw image...";
	project_.saveHDF5(gridData, imageDir + "/final_raw_image", "image", Util::CopyAbsValueOp());
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, imageDir + "/final_image_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, imageDir + "/final_image_z", "z", Util::CopyZOp());
#endif

	project_.showFigure3D(1, "Raw", &gridData, &interfacePointList,
				true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_INVERTED_GRAY);
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::measureSpeed1()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor",      1,      128);
	// This is the base element of the central group.
	const auto baseElement        = methodPM->value<unsigned int>("base_element"     ,      0, config.numElementsMux - config.numElements);
	// elem1 and elem2 are the elements nearest to the group center.
	const auto elem1              = methodPM->value<unsigned int>("element_1"        ,      0, config.numElements - 1);
	const auto elem2              = methodPM->value<unsigned int>("element_2"        ,      0, config.numElements - 1);
	// Distance from the array to the object.
	const auto distance           = methodPM->value<TFloat>(      "distance"         , 1.0e-3, 200.0e-3);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset");
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	Interpolator<TFloat> interp;
	interp.prepare(upsamplingFactor, CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);

	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
	std::vector<TFloat> signal;
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> tempSignal;
	FFTWFilter<TFloat> revRefPulseFilter;

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor);
	interp.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
	std::reverse(revRefPulse.begin(), revRefPulse.end());
	revRefPulseFilter.setCoefficients(revRefPulse);

	const TFloat signalOffset = revRefPulse.size() - 1; // cross-correlation using convolution (revRefPulseFilter)
#else
	HilbertEnvelope<TFloat> envelope;

	const TFloat signalOffset = (config.samplingFrequency * upsamplingFactor) * peakOffset / config.centerFrequency;
#endif
	LOG_DEBUG << "signalOffset: " << signalOffset;

	unsigned int deadZoneSamplesUp = ((config.samplingFrequency * upsamplingFactor) * 2.0 * config.deadZoneM / config.propagationSpeed);

	// Acquire an A-scan to obtain its length.
	acquisition->prepare(baseElement);
	acquisition->execute(elem1, acqData);
	const std::size_t signalLength = acqData.n2() * upsamplingFactor; // when using cross-correlation, the signal is longer than this
	if (deadZoneSamplesUp >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp (" << deadZoneSamplesUp <<
							") >= signalLength (" << signalLength << ").");
	}

	std::vector<unsigned int> txElemList{elem1, elem2};
	std::vector<unsigned int> rxElemList{elem2, elem1};
	std::vector<TFloat> speedList(txElemList.size());
	const TFloat dx = (elem2 - elem1) * config.pitch * 0.5f;
	const TFloat d = std::sqrt(dx * dx + distance * distance) * 2;

	for (unsigned int i = 0; i < txElemList.size(); ++i) {
		LOG_DEBUG << "tx: " << txElemList[i] << " rx: " << rxElemList[i];

		acquisition->execute(txElemList[i], acqData);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
		tempSignal.resize(acqData.n2() * upsamplingFactor);
		interp.interpolate(&acqData(rxElemList[i], 0), acqData.n2(), &tempSignal[0]);

		// Cross-correlation using convolution.
		revRefPulseFilter.filter(tempSignal, signal);
#else
		signal.resize(acqData.n2() * upsamplingFactor);
		interp.interpolate(&acqData(rxElemList[i], 0), acqData.n2(), &signal[0]);

		Util::removeDC(&signal[0], signal.size(), deadZoneSamplesUp);

		envelope.calculate(&signal[0], signal.size());
#endif

		TFloat maxValue = 0;
		unsigned int idxMax = 0;
		for (unsigned int j = deadZoneSamplesUp; j < signal.size(); ++j) {
			if (signal[j] > maxValue) {
				maxValue = signal[j];
				idxMax = j;
			}
		}

		const TFloat tPulseEcho = (idxMax - signalOffset) / (config.samplingFrequency * upsamplingFactor);
		if (tPulseEcho <= 0) {
			THROW_EXCEPTION(InvalidValueException, "Invalid value: tPulseEcho (" << tPulseEcho << ") <= 0.");
		}
		speedList[i] = d / tPulseEcho;
	}

	TFloat speedSum = 0;
	for (unsigned int i = 0; i < txElemList.size(); ++i) {
		LOG_DEBUG << "speed " << i << ": " << std::setprecision(10) << speedList[i];
		speedSum += speedList[i];
	}
	LOG_DEBUG << "mean speed: " << std::setprecision(10) << speedSum / txElemList.size();
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::measureSpeed1AndDistanceError()
{
	//=====================================================================
	// The results are too noisy.

	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor",      1,      128);
	// This is the base element of the central group.
	const auto baseElement        = methodPM->value<unsigned int>("base_element"     ,      0, config.numElementsMux - config.numElements);
	// elem1 and elem2 are the elements nearest to the group center.
	const auto elem1              = methodPM->value<unsigned int>("element_1"        ,      0, config.numElements - 1);
	const auto elem2              = methodPM->value<unsigned int>("element_2"        ,      0, config.numElements - 1);
	const auto elem1a             = methodPM->value<unsigned int>("element_1a"       ,      0, config.numElements - 1);
	const auto elem2a             = methodPM->value<unsigned int>("element_2a"       ,      0, config.numElements - 1);
	// Distance from the array to the object.
	const auto distance           = methodPM->value<TFloat>(      "distance"         , 1.0e-3, 200.0e-3);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset");
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	Interpolator<TFloat> interp;
	interp.prepare(upsamplingFactor, CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);

	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
	std::vector<TFloat> signal;
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> tempSignal;
	FFTWFilter<TFloat> revRefPulseFilter;

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor);
	interp.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
	std::reverse(revRefPulse.begin(), revRefPulse.end());
	revRefPulseFilter.setCoefficients(revRefPulse);

	const unsigned int signalOffset = revRefPulse.size() - 1; // cross-correlation using convolution (revRefPulseFilter)
#else
	HilbertEnvelope<TFloat> envelope;

	const unsigned int signalOffset = (config.samplingFrequency * upsamplingFactor) * peakOffset / config.centerFrequency;
#endif
	LOG_DEBUG << "signalOffset: " << signalOffset;

	unsigned int deadZoneSamplesUp = ((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed);

	// Acquire an A-scan to obtain its length.
	acquisition->prepare(baseElement);
	acquisition->execute(elem1, acqData);
	const std::size_t signalLength = acqData.n2() * upsamplingFactor; // when using cross-correlation, the signal is longer than this
	if (deadZoneSamplesUp >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp (" << deadZoneSamplesUp <<
							") >= signalLength (" << signalLength << ").");
	}

	std::vector<unsigned int> txElemList{elem1, elem2, elem1a, elem2a};
	std::vector<unsigned int> rxElemList{elem2, elem1, elem2a, elem1a};
	const TFloat dx1 = (elem2 - elem1) * config.pitch * 0.5f;
	const TFloat dx2 = (elem2a - elem1a) * config.pitch * 0.5f;

	std::vector<TFloat> tList(txElemList.size());

	for (unsigned int i = 0; i < txElemList.size(); ++i) {
		LOG_DEBUG << "tx: " << txElemList[i] << " rx: " << rxElemList[i];

		acquisition->execute(txElemList[i], acqData);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
		tempSignal.resize(acqData.n2() * upsamplingFactor);
		interp.interpolate(&acqData(rxElemList[i], 0), acqData.n2(), &tempSignal[0]);

		Util::removeDC(&tempSignal[0], tempSignal.size(), deadZoneSamplesUp);

		// Cross-correlation using convolution.
		revRefPulseFilter.filter(tempSignal, signal);
#else
		signal.resize(acqData.n2() * upsamplingFactor);
		interp.interpolate(&acqData(rxElemList[i], 0), acqData.n2(), &signal[0]);

		Util::removeDC(&signal[0], signal.size(), deadZoneSamplesUp);

		envelope.calculate(&signal[0], signal.size());
#endif

		TFloat maxValue = 0;
		unsigned int idxMax = 0;
		for (unsigned int j = deadZoneSamplesUp; j < signal.size(); ++j) {
			if (signal[j] > maxValue) {
				maxValue = signal[j];
				idxMax = j;
			}
		}

		tList[i] = (idxMax > signalOffset) ? (idxMax - signalOffset) / (config.samplingFrequency * upsamplingFactor) : 0;
	}

	const TFloat t1 = (tList[0] + tList[1]) / 2.0f;
	const TFloat t2 = (tList[2] + tList[3]) / 2.0f;
	const TFloat t1_2 = t1 * t1;
	const TFloat t2_2 = t2 * t2;

//	const TFloat a = (dx2 * dx2) / t2_2 - (dx1 * dx1) / t1_2;
//	const TFloat b = 1.0f / t1_2 - 1.0f / t2_2;
//	const TFloat dError = std::sqrt(a / b) - distance;
	const TFloat a = (dx2 * dx2) * t1_2 - (dx1 * dx1) * t2_2;
	const TFloat b = t2_2 - t1_2;
	const TFloat dError = std::sqrt(a / b) - distance;

	const TFloat dz = distance + dError;
	const TFloat c1 = 2.0f * std::sqrt(dz * dz + dx1 * dx1) / t1;
	const TFloat c2 = 2.0f * std::sqrt(dz * dz + dx2 * dx2) / t2;

	LOG_DEBUG << std::setprecision(10) << "distance error: " << dError << " c1: " << c1 << " c2-c1: " << c2 - c1;
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::measureSpeed1AndDistanceError2()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor",      1,      128);
	// This is the base element of the central group.
	const auto baseElement        = methodPM->value<unsigned int>("base_element"     ,      0, config.numElementsMux - config.numElements);
	// elem1 and elem2 are the elements nearest to the group center.
	const auto elem1              = methodPM->value<unsigned int>("element_1"        ,      0, config.numElements - 1);
	const auto elem2              = methodPM->value<unsigned int>("element_2"        ,      0, config.numElements - 1);
	const auto elem1a             = methodPM->value<unsigned int>("element_1a"       ,      0, config.numElements - 1);
	const auto elem2a             = methodPM->value<unsigned int>("element_2a"       ,      0, config.numElements - 1);
	// Distance from the array to the object.
	const auto distance           = methodPM->value<TFloat>(      "distance"         , 1.0e-3, 200.0e-3);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset");
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	Interpolator<TFloat> interp;
	interp.prepare(upsamplingFactor, CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);

	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
	std::vector<TFloat> signal;
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> tempSignal;
	FFTWFilter<TFloat> revRefPulseFilter;

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor);
	interp.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
	std::reverse(revRefPulse.begin(), revRefPulse.end());
	revRefPulseFilter.setCoefficients(revRefPulse);

	const unsigned int signalOffset = revRefPulse.size() - 1; // cross-correlation using convolution (revRefPulseFilter)
#else
	HilbertEnvelope<TFloat> envelope;

	const unsigned int signalOffset = (config.samplingFrequency * upsamplingFactor) * peakOffset / config.centerFrequency;
#endif
	LOG_DEBUG << "signalOffset: " << signalOffset;

	unsigned int deadZoneSamplesUp = ((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed);

	// Acquire an A-scan to obtain its length.
	acquisition->prepare(baseElement);
	acquisition->execute(elem1, acqData);
	const std::size_t signalLength = acqData.n2() * upsamplingFactor; // when using cross-correlation, the signal is longer than this
	if (deadZoneSamplesUp >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp (" << deadZoneSamplesUp <<
							") >= signalLength (" << signalLength << ").");
	}

	std::vector<unsigned int> elem1List, elem2List;
	Util::fillSequence(elem1List, elem1, elem1a, -1);
	Util::fillSequence(elem2List, elem2, elem2a);

	std::vector<TFloat> tList(elem1List.size());
	std::vector<TFloat> dxList(elem1List.size());

	for (unsigned int i = 0; i < elem1List.size(); ++i) {
		LOG_DEBUG << "elem1: " << elem1List[i] << " elem2: " << elem2List[i];
		unsigned int txElemList[2] = {elem1List[i], elem2List[i]};
		unsigned int rxElemList[2] = {elem2List[i], elem1List[i]};
		TFloat t[2];

		for (unsigned int j = 0; j < 2; ++j) {
			acquisition->execute(txElemList[j], acqData);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_SPEED_1_MEASUREMENT_USE_CROSS_CORRELATION
			tempSignal.resize(acqData.n2() * upsamplingFactor);
			interp.interpolate(&acqData(rxElemList[j], 0), acqData.n2(), &tempSignal[0]);

			Util::removeDC(&tempSignal[0], tempSignal.size(), deadZoneSamplesUp);

			// Cross-correlation using convolution.
			revRefPulseFilter.filter(tempSignal, signal);
#else
			signal.resize(acqData.n2() * upsamplingFactor);
			interp.interpolate(&acqData(rxElemList[j], 0), acqData.n2(), &signal[0]);

			Util::removeDC(&signal[0], signal.size(), deadZoneSamplesUp);

			envelope.calculate(&signal[0], signal.size());
#endif

			TFloat maxValue = 0;
			unsigned int idxMax = 0;
			for (unsigned int k = deadZoneSamplesUp; k < signal.size(); ++k) {
				if (signal[k] > maxValue) {
					maxValue = signal[k];
					idxMax = k;
				}
			}

			t[j] = (idxMax > signalOffset) ? (idxMax - signalOffset) / (config.samplingFrequency * upsamplingFactor) : 0;
		}
		tList[i] = (t[0] + t[1]) / 2.0f;
		dxList[i] = (elem2List[i] - elem1List[i]) * config.pitch * 0.5f;
	}

	LOG_DEBUG << std::setprecision(10) << "d = " << distance;
	LOG_DEBUG << "dx_list = [";
	for (unsigned int i = 0; i < elem1List.size(); ++i) {
		LOG_DEBUG << std::setprecision(10) << dxList[i] << ',';
	}
	LOG_DEBUG << ']';
	LOG_DEBUG << "t_list = [";
	for (unsigned int i = 0; i < elem1List.size(); ++i) {
		LOG_DEBUG << std::setprecision(10) << tList[i] << ',';
	}
	LOG_DEBUG << ']';
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::measureDistance()
{
	const ParamMapPtr methodPM = project_.getSubParamMap("method_config_file");
	const STAConfiguration<TFloat> config(*project_.getSubParamMap("main_config_file"));
	const auto savedAcqDir        = methodPM->value<std::string>( "saved_acquisition_dir");
	const auto upsamplingFactor   = methodPM->value<unsigned int>("upsampling_factor", 1, 128);
	// This is the base element of the central group.
	const auto baseElement        = methodPM->value<unsigned int>("base_element"     , 0, config.numElementsMux - config.numElements);
	// elem1 and elem2 are the elements nearest to the group center.
	const auto elem1              = methodPM->value<unsigned int>("element_1"        , 0, config.numElements - 1);
	const auto elem2              = methodPM->value<unsigned int>("element_2"        , 0, config.numElements - 1);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_USE_CROSS_CORRELATION
	const auto referencePulseFile = methodPM->value<std::string>( "reference_pulse_file");
#else
	const auto peakOffset         = methodPM->value<TFloat>(      "peak_offset");
#endif

	std::unique_ptr<STAAcquisition<TFloat>> acquisition =
		std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, savedAcqDir + '/');

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> refPulse;
	project_.loadHDF5(referencePulseFile, "ascan", refPulse);
#endif

	Interpolator<TFloat> interp;
	interp.prepare(upsamplingFactor, CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);

	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
	std::vector<TFloat> signal;
#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_USE_CROSS_CORRELATION
	std::vector<TFloat> tempSignal;
	FFTWFilter<TFloat> revRefPulseFilter;

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor);
	interp.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
	std::reverse(revRefPulse.begin(), revRefPulse.end());
	revRefPulseFilter.setCoefficients(revRefPulse);

	const TFloat signalOffset = revRefPulse.size() - 1; // cross-correlation using convolution (revRefPulseFilter)
#else
	HilbertEnvelope<TFloat> envelope;

	const TFloat signalOffset = (config.samplingFrequency * upsamplingFactor) * peakOffset / config.centerFrequency;
#endif
	LOG_DEBUG << "signalOffset: " << signalOffset;

	unsigned int deadZoneSamplesUp = ((config.samplingFrequency * upsamplingFactor) * 2.0 * config.deadZoneM / config.propagationSpeed);

	// Acquire an A-scan to obtain its length.
	acquisition->prepare(baseElement);
	acquisition->execute(elem1, acqData);
	const std::size_t signalLength = acqData.n2() * upsamplingFactor; // when using cross-correlation, the signal is longer than this
	if (deadZoneSamplesUp >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp (" << deadZoneSamplesUp <<
							") >= signalLength (" << signalLength << ").");
	}

	std::vector<unsigned int> txElemList{elem1, elem2};
	std::vector<unsigned int> rxElemList{elem2, elem1};
	std::vector<TFloat> tList(txElemList.size());
	const TFloat dx = (elem2 - elem1) * config.pitch * 0.5f;

	for (unsigned int i = 0; i < txElemList.size(); ++i) {
		LOG_DEBUG << "tx: " << txElemList[i] << " rx: " << rxElemList[i];

		acquisition->execute(txElemList[i], acqData);

#ifdef CYLINDER_DETECTION_AND_FERMAT_METHOD_DISTANCE_MEASUREMENT_USE_CROSS_CORRELATION
		tempSignal.resize(acqData.n2() * upsamplingFactor);
		interp.interpolate(&acqData(rxElemList[i], 0), acqData.n2(), &tempSignal[0]);

		// Cross-correlation using convolution.
		revRefPulseFilter.filter(tempSignal, signal);
#else
		signal.resize(acqData.n2() * upsamplingFactor);
		interp.interpolate(&acqData(rxElemList[i], 0), acqData.n2(), &signal[0]);

		Util::removeDC(&signal[0], signal.size(), deadZoneSamplesUp);

		envelope.calculate(&signal[0], signal.size());
#endif

		TFloat maxValue = 0;
		unsigned int idxMax = 0;
		for (unsigned int j = deadZoneSamplesUp; j < signal.size(); ++j) {
			if (signal[j] > maxValue) {
				maxValue = signal[j];
				idxMax = j;
			}
		}

		tList[i] = (idxMax - signalOffset) / (config.samplingFrequency * upsamplingFactor);
		if (tList[i] <= 0) {
			THROW_EXCEPTION(InvalidValueException, "Invalid value: tList[i] (" << tList[i] << ") <= 0.");
		}
	}

	for (unsigned int i = 0; i < txElemList.size(); ++i) {
		LOG_DEBUG << "[measureDistance] t" << i << ": " << std::setprecision(10) << tList[i];
	}
	const TFloat t = (tList[0] + tList[1]) * 0.5f;
	const TFloat ct_2 = config.propagationSpeed * t * 0.5f;
	const TFloat d = std::sqrt(ct_2 * ct_2 - dx * dx);

	LOG_DEBUG << "[measureDistance] distance: " << std::setprecision(10) << d;
}

template<typename TFloat>
void
CylinderDetectionAndFermatMethod<TFloat>::execute()
{
	switch (project_.method()) {
	case MethodEnum::cylinder_detection_and_fermat_acquisition_simulated_dp:
	case MethodEnum::cylinder_detection_and_fermat_acquisition_simulated_sp:
	case MethodEnum::cylinder_detection_and_fermat_acquisition_network_dp:
	case MethodEnum::cylinder_detection_and_fermat_acquisition_network_sp:
		acquireSignals();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_sp:
		detectPoints();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_vectorial_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_vectorial_sp:
		detectPointsUsingVectorSum();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_cc_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_cc_sp:
		detectPointsUsingCrossCorrelation();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_ccbf_pulse_echo_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_ccbf_pulse_echo_sp:
		detectPointsUsingCCBFPulseEcho();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_ccbf_pitch_catch_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_ccbf_pitch_catch_sp:
		detectPointsUsingCCBFPitchCatch();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_arc_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_arc_sp:
		detectPointsInArcs();
		break;
	case MethodEnum::cylinder_detection_and_fermat_point_detection_tangent_curve_dp:
	case MethodEnum::cylinder_detection_and_fermat_point_detection_tangent_curve_sp:
		detectPointsUsingTangentCurveGeometry();
		break;
	case MethodEnum::cylinder_detection_and_fermat_circle_fitting_dp:
	case MethodEnum::cylinder_detection_and_fermat_circle_fitting_sp:
		fitCircle();
		break;
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_dp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_sp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_vectorial_dp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_vectorial_sp:
		execTwoMediumImaging();
		break;
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_dp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_sp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_sta_dp:
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_sta_sp:
		execCombinedTwoMediumImaging();
		break;
#ifdef USE_OPENCL
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_ocl_sp:
		if constexpr (std::is_same<TFloat, float>::value) {
			execCombinedTwoMediumImagingCyl<VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>>();
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
		break;
#endif
#ifdef USE_CUDA
	case MethodEnum::cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_cuda_sp:
		if constexpr (std::is_same<TFloat, float>::value) {
			execCombinedTwoMediumImagingCyl<VectorialCombinedTwoMediumImagingCUDAProcessor>();
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
		break;
#endif
	case MethodEnum::cylinder_detection_and_fermat_speed_1_measurement:
		measureSpeed1();
		measureSpeed1AndDistanceError();
		measureSpeed1AndDistanceError2();
		break;
	case MethodEnum::cylinder_detection_and_fermat_distance_measurement:
		measureDistance();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

} // namespace Lab

#endif /* CYLINDERDETECTIONANDFERMATMETHOD_H_ */
