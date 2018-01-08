#include "SingleAcquisitionMethod.h"

#include <cstddef> /* std::size_t */
#include <sstream>
#include <vector>

#include "ArrayAcqClient.h"
#include "global.h"
#include "Log.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"

//#define TEST_MODE 1



namespace Lab {

void
SingleAcquisitionMethod::fillConfiguration()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	config_.samplingFrequency = taskPM->value<double>(      "sampling_frequency",    100.0,  200.0e6);
	config_.centerFrequency   = taskPM->value<double>(      "center_frequency"  ,    100.0,  100.0e6);
	config_.acquisitionTime   = taskPM->value<double>(      "acquisition_time"  ,   1.0e-6,      1.0);
	config_.minGain           = taskPM->value<double>(      "min_gain"          ,      0.0,     48.0);
	//config_.maxGain           = taskPM->value<double>(      "max_gain"          ,      0.0,     48.0);
	config_.numElementsMux    = taskPM->value<unsigned int>("num_elements_mux"  ,        8,     1024);
	config_.numElements       = taskPM->value<unsigned int>("num_elements"      ,        8, config_.numElementsMux);
	config_.baseElement       = taskPM->value<unsigned int>("base_element"      ,        0, config_.numElementsMux - config_.numElements);
	config_.txGroupElement    = taskPM->value<unsigned int>("tx_group_element"  ,        0, config_.numElements);
	config_.rxGroupElement    = taskPM->value<unsigned int>("rx_group_element"  ,        0, config_.numElements);
	config_.averageN          = taskPM->value<unsigned int>("average_n"         ,        1,      256);
	config_.savedAcqDir       = taskPM->value<std::string>( "saved_acquisition_dir");
}

SingleAcquisitionMethod::SingleAcquisitionMethod(Project& project)
		: project_(project)
{
	fillConfiguration();

	ConstParameterMapPtr pm = project_.loadParameterMap("network_acquisition.lab4config");
	std::string serverIpAddress = pm->value<std::string>("server_ip_address");
	unsigned short portNumber = pm->value<unsigned short>("server_port_number", 49152, 65535);

#ifndef TEST_MODE
	acq_.reset(new ArrayAcqClient(serverIpAddress.c_str(), portNumber));

	LOG_DEBUG << "max. sample value = " << acq_->getMaxSampleValue();
	LOG_DEBUG << "min. sample value = " << acq_->getMinSampleValue();

	acq_->execPreConfiguration();

	acq_->setSamplingFrequency(config_.samplingFrequency);
	acq_->setCenterFrequency(config_.centerFrequency);

	std::string rxMask(config_.numElements, '1');
	acq_->setActiveReceiveElements(rxMask);

	acq_->setGain(config_.minGain);
	acq_->setAcquisitionTime(config_.acquisitionTime);

	acq_->execPostConfiguration();
#endif
}

SingleAcquisitionMethod::~SingleAcquisitionMethod()
{
}

void
SingleAcquisitionMethod::execute()
{
#ifndef TEST_MODE
	const std::size_t ascanLength = acq_->getAscanLength();
	std::vector<double> ascan(ascanLength);
	std::vector<float> ascanBuffer;

	acq_->setBaseElement(config_.baseElement);

	std::string txMask(config_.numElements, '0');
	txMask[config_.txGroupElement] = '1';
	acq_->setActiveTransmitElements(txMask);

	for (unsigned int i = 0; i < config_.averageN; ++i) {
		LOG_DEBUG << "ACQ " << i;

		acq_->getAscan(ascanBuffer);

		const float* ascanBase = &ascanBuffer[config_.rxGroupElement * ascanLength];
		Util::addElements(ascanBase, ascan.begin(), ascan.end());
	}

	if (config_.averageN > 1) {
		double factor = 1.0 / config_.averageN;
		Util::multiply(ascan, factor);
	}

	std::vector<double> t(ascanLength);
	double dt = 1.0 / config_.samplingFrequency;
	for (std::size_t i = 0; i < ascanLength; ++i) {
		t[i] = i * dt;
	}

#else
# define ASCAN_LENGTH 1000000
	std::vector<double> ascan(ASCAN_LENGTH), t(ASCAN_LENGTH), b(ASCAN_LENGTH);
	for (int i = 0; i < ASCAN_LENGTH / 4; ++i) {
		ascan[i] = i;
		t[i] = i * 1e-12;
		b[i] = -ascan[i] * 0.5;
	}
	for (int i = ASCAN_LENGTH / 4; i < ASCAN_LENGTH * 3 / 4; ++i) {
		ascan[i] = ASCAN_LENGTH / 2 - i;
		t[i] = i * 1e-12;
		b[i] = -ascan[i] * 0.5;
	}
	for (int i = ASCAN_LENGTH * 3 / 4; i < ASCAN_LENGTH; ++i) {
		ascan[i] = i - ASCAN_LENGTH;
		t[i] = i * 1e-12;
		b[i] = -ascan[i] * 0.5;
	}
#endif

	{
		LOG_DEBUG << "Saving the A-scan...";
		std::ostringstream fileName;
		fileName << config_.savedAcqDir << "/base" << config_.baseElement << "_tx" << config_.txGroupElement << "_rx" << config_.rxGroupElement;
		project_.saveHDF5(ascan, fileName.str(), "ascan");
	}

	project_.showFigure2D(0, "A-scan", t, ascan);

#if 0
	LOG_DEBUG << "Executing the plot script...";
	std::vector<std::string> args;
	std::ostringstream arg1;
	arg1 << config_.savedAcqDir << '/' << config_.ascanFile << ".h5";
	args.push_back(arg1.str());
	std::ostringstream arg2;
	arg2.precision(FLOAT_SCIENTIFIC_NOTATION_NUM_DIGITS_AFTER_DECIMAL_POINT);
	arg2 << std::scientific << config_.samplingFrequency;
	args.push_back(arg2.str());
	project_.executeProgram(config_.plotScript, args);
#endif
}

} // namespace Lab
