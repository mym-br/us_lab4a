#ifndef NETWORKSECTORIALSCANACQUISITION_H
#define NETWORKSECTORIALSCANACQUISITION_H

#include <cmath>
#include <memory>
#include <string>

#include <boost/cstdint.hpp>

#include "ParameterMap.h"
#include "PhasedArrayAcqClient.h"
#include "Project.h"
#include "SectorialScanAcquisition.h"
#include "SectorialScanConfiguration.h"



namespace Lab {

template<typename FloatType>
class NetworkSectorialScanAcquisition : public SectorialScanAcquisition<FloatType> {
public:
	NetworkSectorialScanAcquisition(
		const Project& project,
		const SectorialScanConfiguration<FloatType>& config);
	virtual ~NetworkSectorialScanAcquisition();

	virtual void execute(typename SectorialScanAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	NetworkSectorialScanAcquisition(const NetworkSectorialScanAcquisition&);
	NetworkSectorialScanAcquisition& operator=(const NetworkSectorialScanAcquisition&);

	const Project& project_;
	const SectorialScanConfiguration<FloatType>& config_;
	std::unique_ptr<PhasedArrayAcqClient> acq_;
	std::vector<float> imageBuffer_;
	std::vector<float> lineStartX_;
	std::vector<float> lineStartZ_;
	std::vector<float> lineAngle_;
};



template<typename FloatType>
NetworkSectorialScanAcquisition<FloatType>::NetworkSectorialScanAcquisition(const Project& project, const SectorialScanConfiguration<FloatType>& config)
		: project_(project)
		, config_(config)
{
	ConstParameterMapPtr pm = project_.loadParameterMap("config-network_acquisition.txt");
	std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	unsigned short portNumber   = pm->value<unsigned short>("server_port_number", 49152, 65535);

	acq_.reset(new PhasedArrayAcqClient(serverIpAddress.c_str(), portNumber));

	acq_->execPreConfiguration();

	acq_->setAcquisitionDelay(config_.acquisitionDelay);
	acq_->setSignalMode(config_.signalMode);
	acq_->setMaterialVelocity(config_.propagationSpeed);
	acq_->setFocalPoint(config_.focalEmissionDistance, config_.focalReceptionDistance);
	acq_->setPhasedArrayConfiguration(config_.numElementsMux, config_.pitch, config_.centerFrequency);
	acq_->setRange(config_.rangeStart, config_.rangeEnd);
	acq_->setSectorialScan(config_.baseElement, config_.numElements, config_.startAngle, config_.endAngle, config_.angleStep);
	acq_->setGain(config_.gain);
	LOG_DEBUG << "[NetworkSectorialScanAcquisition] sampling frequency: " << acq_->getSamplingFrequency();

	acq_->execPostConfiguration();
}

template<typename FloatType>
NetworkSectorialScanAcquisition<FloatType>::~NetworkSectorialScanAcquisition()
{
}

template<typename FloatType>
void
NetworkSectorialScanAcquisition<FloatType>::execute(typename SectorialScanAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ";

	acq_->getImageBuffer(imageBuffer_);

	boost::uint32_t numRows = acq_->getImageNumRows();
	boost::uint32_t numCols = acq_->getImageNumCols();
	acqData.resize(numCols, numRows); // image lines are stored in each row of acqData

	acq_->getImageLineGeometry(lineStartX_, lineStartZ_, lineAngle_);

	for (unsigned int j = 0; j < numCols; ++j) {
		const FloatType angle = lineAngle_[j];
		const FloatType sa = std::sin(angle);
		const FloatType ca = std::cos(angle);
		for (unsigned int i = 0; i < numRows; ++i) {
			const FloatType coef = i / (numRows - FloatType{1.0});
			const FloatType r = (config_.rangeEnd - config_.rangeStart) * coef + config_.rangeStart;
			XZValue<FloatType>& point = acqData(j, i);
			point.x = lineStartX_[j] * 1.0e-3f + sa * r;
			point.z = lineStartZ_[j] * 1.0e-3f + ca * r;
			point.value = imageBuffer_[i + j * numRows];
		}
	}
}

} // namespace Lab

#endif // NETWORKSECTORIALSCANACQUISITION_H
