/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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

#ifndef SIMRECTANGULARFLATSOURCEMETHOD_H
#define SIMRECTANGULARFLATSOURCEMETHOD_H

#include <cmath>
#include <memory>

#include "Exception.h"
#include "Log.h"
#include "Matrix2.h"
#include "Method.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SimTransientAcousticBeam.h"
#include "Timer.h"
#include "Util.h"
#include "Waveform.h"
#include "XYZ.h"
#include "XZ.h"
#include "XZValue.h"



namespace Lab {

template<typename FloatType>
class SimRectangularFlatSourceMethod : public Method {
public:
	SimRectangularFlatSourceMethod(Project& project);
	virtual ~SimRectangularFlatSourceMethod();

	virtual void execute();
private:
	SimRectangularFlatSourceMethod(const SimRectangularFlatSourceMethod&) = delete;
	SimRectangularFlatSourceMethod& operator=(const SimRectangularFlatSourceMethod&) = delete;

	void execTransientAcousticBeam();

	Project& project_;
};



template<typename FloatType>
SimRectangularFlatSourceMethod<FloatType>::SimRectangularFlatSourceMethod(Project& project)
		: project_{project}
{
}

template<typename FloatType>
SimRectangularFlatSourceMethod<FloatType>::~SimRectangularFlatSourceMethod()
{
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientAcousticBeam()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");

	const FloatType beamDistance     = taskPM->value<FloatType>("beam_distance", 0.0, 100.0);
	const FloatType beamThetaYStep   = taskPM->value<FloatType>("beam_theta_y_step", 0.0, 10.0);
	const FloatType beamThetaYMax    = taskPM->value<FloatType>("beam_theta_y_max" , 0.0, 90.0);
	const FloatType beamThetaXStep   = taskPM->value<FloatType>("beam_theta_x_step", 0.0, 10.0);
	const FloatType beamThetaXMax    = taskPM->value<FloatType>("beam_theta_x_max" , 0.0, 90.0);
	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency", 0.0, 1.0e12);
	const FloatType subElemSize      = taskPM->value<FloatType>("sub_elem_size", 0.0, 1.0);
	const std::string excitationType = taskPM->value<std::string>("excitation_type");

	std::vector<FloatType> exc;
	if (excitationType == "1") {
		Waveform::getType1(centerFreq, samplingFreq, exc,
					taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0));
	} else if (excitationType == "2a") {
		Waveform::getType2a(centerFreq, samplingFreq, exc);
	} else if (excitationType == "2b") {
		Waveform::getType2b(centerFreq, samplingFreq, exc);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType << '.');
	}

	const FloatType dt = 1.0 / samplingFreq;
	std::vector<FloatType> tExc(exc.size());
	for (unsigned int i = 0; i < tExc.size(); ++i) {
		tExc[i] = dt * i;
	}
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, 1.0 / samplingFreq, dvdt);

	std::vector<FloatType> thetaXList;
	Util::fillSequenceWithSize(thetaXList, 0.0, beamThetaXMax, std::ceil(beamThetaXMax / beamThetaXStep) + 1);
	std::vector<FloatType> thetaYList;
	Util::fillSequenceWithSize(thetaYList, 0.0, beamThetaYMax, std::ceil(beamThetaYMax / beamThetaYStep) + 1);

	Matrix2<XZValue<FloatType>> gridData{thetaXList.size(), thetaYList.size()};
	Matrix2<XYZ<FloatType>> inputData{thetaXList.size(), thetaYList.size()};

	for (unsigned int ix = 0, xSize = thetaXList.size(); ix < xSize; ++ix) {
		const FloatType tX = Util::degreeToRadian(thetaXList[ix]);
		const FloatType sinTX = std::sin(tX);
		const FloatType cosTX = std::cos(tX);
		for (unsigned int iy = 0, ySize = thetaYList.size(); iy < ySize; ++iy) {
			const FloatType tY = Util::degreeToRadian(thetaYList[iy]);
			const FloatType x = beamDistance * std::sin(tY);
			const FloatType rx = beamDistance * std::cos(tY);
			const FloatType y = rx * sinTX;
			const FloatType z = rx * cosTX;
			XZValue<FloatType>& gd = gridData(ix, iy);
			gd.x = thetaYList[iy];
			gd.z = thetaXList[ix];
			gd.value = 0.0;
			XYZ<FloatType>& id = inputData(ix, iy);
			id.x = x;
			id.y = y;
			id.z = z;
		}
	}

	auto acBeam = std::make_unique<SimTransientAcousticBeam<FloatType>>();
	acBeam->getRectangularFlatSourceAcousticBeam(
				sourceWidth,
				sourceHeight,
				samplingFreq,
				propagationSpeed,
				subElemSize,
				dvdt,
				inputData,
				gridData);

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XZ<float>> pointList = {{0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXZValue(gridData, projGridData);
	project_.showFigure3D(1, "Beam", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execute()
{
	Timer tProc;

	switch (project_.method()) {
	case MethodType::sim_acoustic_beam_rectangular_flat_source_transient:
		execTransientAcousticBeam();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();
}

} // namespace Lab

#endif // SIMRECTANGULARFLATSOURCEMETHOD_H
