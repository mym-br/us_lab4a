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

	const STAConfiguration<FloatType> config(taskPM);
	const std::string outputDir    = taskPM->value<std::string>( "output_dir");
	const unsigned int baseElement = taskPM->value<unsigned int>("base_element",   0, config.numElementsMux - config.numElements);
	const FloatType peakOffset     = taskPM->value<FloatType>(   "peak_offset" , 0.0, 50.0);

	boost::scoped_ptr<STAAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case Method::STA_SECTORIAL_SIMPLE_SIMULATED:   // falls through
	case Method::STA_SECTORIAL_SIMULATED:
		acquisition.reset(new SimulatedSTAAcquisition<FloatType>(project_, config));
		break;
	case Method::STA_SECTORIAL_DP_NETWORK:
		acquisition.reset(new NetworkSTAAcquisition<FloatType>(project_, config));
		break;
	case Method::STA_SECTORIAL_SIMPLE_SAVED:       // falls through
	case Method::STA_SECTORIAL_DP_SAVED:           // falls through
	case Method::STA_SECTORIAL_VECTORIAL_DP_SAVED: // falls through
	case Method::STA_SECTORIAL_SP_SAVED:           // falls through
	case Method::STA_SECTORIAL_VECTORIAL_SP_SAVED:
		acquisition.reset(new SavedSTAAcquisition<FloatType>(project_, config.numElements));
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << project_.method() << '.');
	}

	Matrix2<XZValueFactor<FloatType>> gridData;

	const FloatType lambda = config.propagationSpeed / config.centerFrequency;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), lambda, gridData);

	bool coherenceFactorEnabled;
	bool vectorialProcessingWithEnvelope = false;
	if (project_.method() == Method::STA_SECTORIAL_VECTORIAL_DP_SAVED ||
			project_.method() == Method::STA_SECTORIAL_VECTORIAL_SP_SAVED) {
		vectorialProcessingWithEnvelope = taskPM->value<bool>("calculate_envelope_in_processing");
	}

	switch (project_.method()) {
	case Method::STA_SECTORIAL_SIMPLE_SIMULATED: // falls through
	case Method::STA_SECTORIAL_SIMPLE_SAVED:
		{
			boost::scoped_ptr<STAProcessor<FloatType> > processor;
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			processor.reset(new SimpleSTAProcessor<FloatType>(config, *acquisition, peakOffset));
			processor->process(baseElement, gridData);
		}
		break;
	case Method::STA_SECTORIAL_VECTORIAL_DP_SAVED:
	case Method::STA_SECTORIAL_VECTORIAL_SP_SAVED:
		{
			const unsigned int upsamplingFactor = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
			boost::scoped_ptr<STAProcessor<FloatType> > processor;
			AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			processor.reset(new VectorialSTAProcessor<FloatType>(config, *acquisition, upsamplingFactor, coherenceFactor, peakOffset, vectorialProcessingWithEnvelope));
			processor->process(baseElement, gridData);
		}
		break;
	default:
		{
			boost::scoped_ptr<STAProcessor<FloatType> > processor;
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			processor.reset(new DefaultSTAProcessor<FloatType>(config, *acquisition, coherenceFactor, peakOffset));
			processor->process(baseElement, gridData);
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

	Project::GridDataType projGridData;
	Util::copyXZValue(gridData, projGridData);
	project_.showFigure3D(1, "Raw", &projGridData, Project::emptyPointList,
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
		project_.showFigure3D(2, "CF", &projGridData, Project::emptyPointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
	}
}

} // namespace Lab

#endif /* STAMETHOD_H_ */
