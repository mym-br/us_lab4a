#ifndef STAMETHOD_H_
#define STAMETHOD_H_

#include <algorithm>

#include <boost/scoped_ptr.hpp>

#include "CoherenceFactor.h"
#include "DefaultSTAProcessor.h"
#include "HilbertEnvelope.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Method.h"
#include "NetworkSTAAcquisition.h"
#include "ParallelHilbertEnvelope.h"
#include "Project.h"
#include "SavedSTAAcquisition.h"
#include "SimpleSTAProcessor.h"
#include "SimulatedSTAAcquisition.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "STAProcessor.h"
#include "Timer.h"
#include "VectorialSTAProcessor.h"
#include "Util.h"
#include "XZValue.h"
#include "XZValueFactor.h"



namespace Lab {

template<typename FloatType>
class STAMethod : public Method {
public:
	STAMethod(Project& project);
	virtual ~STAMethod();

	virtual void execute();

private:
	STAMethod(const STAMethod&);
	STAMethod& operator=(const STAMethod&);

	void saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq, unsigned int baseElement, const std::string& dataDir);

	Project& project_;
};



template<typename FloatType>
STAMethod<FloatType>::STAMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
STAMethod<FloatType>::~STAMethod()
{
}

template<typename FloatType>
void
STAMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const STAConfiguration<FloatType> config(project_.loadChildParameterMap(taskPM, "sta_config_file"));
	const unsigned int baseElement = taskPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);

	boost::scoped_ptr<STAAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case MethodType::sta_sectorial_simple_simulated: // falls through
	case MethodType::sta_sectorial_simulated:
		acquisition.reset(new SimulatedSTAAcquisition<FloatType>(project_, config));
		break;
	case MethodType::sta_sectorial_dp_network: // falls through
	case MethodType::sta_save_signals:
		acquisition.reset(new NetworkSTAAcquisition<FloatType>(project_, config));
		break;
	case MethodType::sta_sectorial_simple_saved:       // falls through
	case MethodType::sta_sectorial_dp_saved:           // falls through
	case MethodType::sta_sectorial_vectorial_dp_saved: // falls through
	case MethodType::sta_sectorial_sp_saved:           // falls through
	case MethodType::sta_sectorial_vectorial_sp_saved:
		acquisition.reset(new SavedSTAAcquisition<FloatType>(project_, config.numElements));
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodType::sta_save_signals) {
		std::string dataDir = taskPM->value<std::string>("data_dir");
		saveSignals(config, *acquisition, baseElement, dataDir);
		return;
	}

	Matrix2<XZValueFactor<FloatType>> gridData;

	const FloatType lambda = config.propagationSpeed / config.centerFrequency;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), lambda, gridData);

	bool coherenceFactorEnabled = false;
	bool vectorialProcessingWithEnvelope = false;
	if (project_.method() == MethodType::sta_sectorial_vectorial_dp_saved ||
			project_.method() == MethodType::sta_sectorial_vectorial_sp_saved) {
		vectorialProcessingWithEnvelope = taskPM->value<bool>("calculate_envelope_in_processing");
	}

	const FloatType peakOffset  = taskPM->value<FloatType>(   "peak_offset", 0.0, 50.0);
	const std::string outputDir = taskPM->value<std::string>( "output_dir");

	switch (project_.method()) {
	case MethodType::sta_sectorial_simple_simulated: // falls through
	case MethodType::sta_sectorial_simple_saved:
		{
			boost::scoped_ptr<STAProcessor<FloatType> > processor;
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			processor.reset(new SimpleSTAProcessor<FloatType>(config, *acquisition, peakOffset));

			Timer tProc;
			processor->process(baseElement, gridData);
			LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
		}
		break;
	case MethodType::sta_sectorial_vectorial_dp_saved: // falls through
	case MethodType::sta_sectorial_vectorial_sp_saved:
		{
			const unsigned int upsamplingFactor = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
			boost::scoped_ptr<STAProcessor<FloatType> > processor;
			AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			processor.reset(new VectorialSTAProcessor<FloatType>(config, *acquisition, upsamplingFactor, coherenceFactor, peakOffset, vectorialProcessingWithEnvelope));

			Timer tProc;
			processor->process(baseElement, gridData);
			LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
		}
		break;
	default:
		{
			boost::scoped_ptr<STAProcessor<FloatType> > processor;
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			processor.reset(new DefaultSTAProcessor<FloatType>(config, *acquisition, coherenceFactor, peakOffset));

			Timer tProc;
			processor->process(baseElement, gridData);
			LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
		}
	}

	{
		Matrix2<FloatType> aux;
		Util::copyValueToSimpleMatrix(gridData, aux);
		LOG_DEBUG << "Saving the image...";
		project_.saveHDF5(aux, outputDir + "/raw_image", "image");

		Matrix2<double> gridX, gridZ;
		Util::copyXZToSimpleMatrices(gridData, gridX, gridZ);
		LOG_DEBUG << "Saving the X coordinates...";
		project_.saveHDF5(gridX, outputDir + "/image_x", "x");
		LOG_DEBUG << "Saving the Z coordinates...";
		project_.saveHDF5(gridZ, outputDir + "/image_z", "z");
	}

	Figure::Visualization visual;
	if (vectorialProcessingWithEnvelope) {
		visual = Figure::VISUALIZATION_RECTIFIED_LOG;
	} else {
		visual = Figure::VISUALIZATION_ENVELOPE_LOG;
	}

	std::vector<XZ<float>> pointList = {{0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXZValue(gridData, projGridData);
	project_.showFigure3D(1, "Raw", &projGridData, &pointList,
				true, visual, Figure::COLORMAP_VIRIDIS);

	if (coherenceFactorEnabled) {
		{
			Matrix2<FloatType> aux;
			Util::copyFactorToSimpleMatrix(gridData, aux);
			LOG_DEBUG << "Saving the CF factors...";
			project_.saveHDF5(aux, outputDir + "/cf_factors", "image");
		}

		ParallelHilbertEnvelope<FloatType>::calculateDim2(gridData);

		// Applies the coherence factor method.
		for (typename Matrix2<XZValueFactor<FloatType> >::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->value *= iter->factor;
			iter->factor = 1.0;
		}
		{
			Matrix2<FloatType> aux;
			Util::copyValueToSimpleMatrix(gridData, aux);
			LOG_DEBUG << "Saving the CF image...";
			project_.saveHDF5(aux, outputDir + "/cf_image", "image");
		}

		Util::copyXZValue(gridData, projGridData);
		project_.showFigure3D(2, "CF", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
	}
}

template<typename FloatType>
void
STAMethod<FloatType>::saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq, unsigned int baseElement, const std::string& dataDir)
{
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData;

	for (unsigned int txElem = 0; txElem < config.numElements; ++txElem) {
		acq.execute(baseElement, txElem, acqData);

		std::ostringstream filePath;
		filePath << dataDir << std::setfill('0') << "/signal-base" << std::setw(4) << baseElement << "-tx" << std::setw(4) << txElem;
		std::string filePathStr = filePath.str();
		LOG_DEBUG << "Saving " << filePathStr << "...";
		project_.saveHDF5(acqData, filePathStr, "signal");
	}
}

} // namespace Lab

#endif /* STAMETHOD_H_ */
