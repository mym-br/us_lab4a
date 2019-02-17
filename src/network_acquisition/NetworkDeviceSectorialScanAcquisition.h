/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
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

#ifndef NETWORKDEVICESECTORIALSCANACQUISITION_H
#define NETWORKDEVICESECTORIALSCANACQUISITION_H

#include <cmath>
#include <memory>
#include <string>

#include <boost/cstdint.hpp>

#include "DeviceSectorialScanAcquisition.h"
#include "DeviceSectorialScanConfiguration.h"
#include "global.h"
#include "ParameterMap.h"
#include "PhasedArrayAcqClient.h"
#include "Project.h"
#include "Util.h"



namespace Lab {

template<typename FloatType>
class NetworkDeviceSectorialScanAcquisition : public DeviceSectorialScanAcquisition<FloatType> {
public:
	NetworkDeviceSectorialScanAcquisition(
		const Project& project,
		const DeviceSectorialScanConfiguration<FloatType>& config);
	virtual ~NetworkDeviceSectorialScanAcquisition();

	virtual void execute(typename DeviceSectorialScanAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	NetworkDeviceSectorialScanAcquisition(const NetworkDeviceSectorialScanAcquisition&);
	NetworkDeviceSectorialScanAcquisition& operator=(const NetworkDeviceSectorialScanAcquisition&);

	const Project& project_;
	const DeviceSectorialScanConfiguration<FloatType>& config_;
	std::unique_ptr<PhasedArrayAcqClient> acq_;
	std::vector<float> imageBuffer_;
	std::vector<float> lineStartX_;
	std::vector<float> lineStartZ_;
	std::vector<float> lineAngle_;
	float valueFactor_;
};



template<typename FloatType>
NetworkDeviceSectorialScanAcquisition<FloatType>::NetworkDeviceSectorialScanAcquisition(const Project& project, const DeviceSectorialScanConfiguration<FloatType>& config)
		: project_(project)
		, config_(config)
{
	ConstParameterMapPtr pm = project_.loadParameterMap(NETWORK_AQUISITION_CONFIG_FILE);
	const std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	const unsigned short portNumber   = pm->value<unsigned short>("server_port_number",   49152,  65535);
	valueFactor_               = 1.0f / pm->value<float>(         "value_scale"       , 1.0e-30, 1.0e30);

	acq_ = std::make_unique<PhasedArrayAcqClient>(serverIpAddress.c_str(), portNumber);

	acq_->execPreConfiguration();

	acq_->setAcquisitionDelay(config_.acquisitionDelay);
	acq_->setSignalMode(config_.signalMode);
	acq_->setMaterialVelocity(config_.propagationSpeed);
	if (config_.enableFocusing) {
		acq_->setFocalPoint(config_.focalEmissionDistance, config_.focalReceptionDistance);
	} else {
		acq_->setFocalPoint(0.0, 0.0);
	}
	acq_->setPhasedArrayConfiguration(config_.numElementsMux, config_.pitch, config_.centerFrequency);
	acq_->setRange(config_.rangeStart, config_.rangeEnd);
	acq_->setSectorialScan(config_.baseElement, config_.numElements, config_.startAngle, config_.endAngle, config_.angleStep);
	acq_->setGain(config_.gain);
	LOG_DEBUG << "[NetworkSectorialScanAcquisition] sampling frequency: " << acq_->getSamplingFrequency();

	acq_->execPostConfiguration();
}

template<typename FloatType>
NetworkDeviceSectorialScanAcquisition<FloatType>::~NetworkDeviceSectorialScanAcquisition()
{
}

template<typename FloatType>
void
NetworkDeviceSectorialScanAcquisition<FloatType>::execute(typename DeviceSectorialScanAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ";

	acq_->getImageBuffer(imageBuffer_);
	Util::multiply(imageBuffer_, valueFactor_);

	boost::uint32_t numRows = acq_->getImageNumRows();
	boost::uint32_t numCols = acq_->getImageNumCols();
	acqData.resize(numCols, numRows); // image lines are stored in each row of acqData

	acq_->getImageLineGeometry(lineStartX_, lineStartZ_, lineAngle_);

	for (unsigned int j = 0; j < numCols; ++j) {
		const FloatType angle = lineAngle_[j];
		const FloatType sa = std::sin(angle);
		const FloatType ca = std::cos(angle);
		for (unsigned int i = 0; i < numRows; ++i) {
			const FloatType coef = i / (numRows - FloatType(1.0));
			const FloatType r = (config_.rangeEnd - config_.rangeStart) * coef + config_.rangeStart;
			XZValue<FloatType>& point = acqData(j, i);
			point.x = lineStartX_[j] * 1.0e-3f + sa * r;
			point.z = lineStartZ_[j] * 1.0e-3f + ca * r;
			point.value = imageBuffer_[i + j * numRows];
		}
	}
}

} // namespace Lab

#endif // NETWORKDEVICESECTORIALSCANACQUISITION_H
