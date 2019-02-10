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

#include "SingleAcquisitionMethod.h"

#include <cstddef> /* std::size_t */
#include <iomanip>
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
	config_.minGain           = taskPM->value<double>(      "min_gain"          ,      0.0,     100.0);
	//config_.maxGain           = taskPM->value<double>(      "max_gain"          ,      0.0,     100.0);
	config_.numElementsMux    = taskPM->value<unsigned int>("num_elements_mux"  ,        8,     1024);
	config_.numElements       = taskPM->value<unsigned int>("num_elements"      ,        8, config_.numElementsMux);
	config_.baseElement       = taskPM->value<unsigned int>("base_element"      ,        0, config_.numElementsMux - config_.numElements);
	config_.txGroupElement    = taskPM->value<unsigned int>("tx_group_element"  ,        0, config_.numElements);
	config_.rxGroupElement    = taskPM->value<unsigned int>("rx_group_element"  ,        0, config_.numElements);
	config_.numPulses         = taskPM->value<unsigned int>("num_pulses"        ,        1, 100);
	config_.savedAcqDir       = taskPM->value<std::string>( "saved_acquisition_dir");
}

SingleAcquisitionMethod::SingleAcquisitionMethod(Project& project)
		: project_(project)
{
	fillConfiguration();

	project_.createDirectory(config_.savedAcqDir, false);

	ConstParameterMapPtr pm = project_.loadParameterMap(NETWORK_AQUISITION_CONFIG_FILE);
	const std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	const unsigned short portNumber   = pm->value<unsigned short>("server_port_number",   49152,  65535);
	valueFactor_                = 1.0 / pm->value<double>(        "value_scale"       , 1.0e-30, 1.0e30);
	averageN_                         = pm->value<unsigned int>(  "average_n"         ,       1,    256);

#ifndef TEST_MODE
	acq_ = std::make_unique<ArrayAcqClient>(serverIpAddress.c_str(), portNumber);

	LOG_DEBUG << "max. sample value = " << acq_->getMaxSampleValue();
	LOG_DEBUG << "min. sample value = " << acq_->getMinSampleValue();

	acq_->execPreConfiguration();

	acq_->setSamplingFrequency(config_.samplingFrequency);
	acq_->setCenterFrequency(config_.centerFrequency, config_.numPulses);

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
	const std::size_t signalLength = acq_->getSignalLength();
	std::vector<double> signal(signalLength);
	std::vector<float> signalBuffer;

	acq_->setBaseElement(config_.baseElement);

	std::string txMask(config_.numElements, '0');
	txMask[config_.txGroupElement] = '1';
	acq_->setActiveTransmitElements(txMask);

	for (unsigned int i = 0; i < averageN_; ++i) {
		LOG_DEBUG << "ACQ " << i;

		acq_->getSignal(signalBuffer);

		const float* signalBase = &signalBuffer[config_.rxGroupElement * signalLength];
		Util::addElements(signalBase, signal.begin(), signal.end());
	}

	if (averageN_ > 1U) {
		const double factor = valueFactor_ / averageN_;
		Util::multiply(signal, factor);
	} else {
		Util::multiply(signal, valueFactor_);
	}

	std::vector<double> t(signalLength);
	double dt = 1.0 / config_.samplingFrequency;
	for (std::size_t i = 0; i < signalLength; ++i) {
		t[i] = i * dt;
	}

#else
# define SIGNAL_LENGTH 1000000
	std::vector<double> signal(SIGNAL_LENGTH), t(SIGNAL_LENGTH), b(SIGNAL_LENGTH);
	for (int i = 0; i < SIGNAL_LENGTH / 4; ++i) {
		signal[i] = i;
		t[i] = i * 1e-12;
		b[i] = -signal[i] * 0.5;
	}
	for (int i = SIGNAL_LENGTH / 4; i < SIGNAL_LENGTH * 3 / 4; ++i) {
		signal[i] = SIGNAL_LENGTH / 2 - i;
		t[i] = i * 1e-12;
		b[i] = -signal[i] * 0.5;
	}
	for (int i = SIGNAL_LENGTH * 3 / 4; i < SIGNAL_LENGTH; ++i) {
		signal[i] = i - SIGNAL_LENGTH;
		t[i] = i * 1e-12;
		b[i] = -signal[i] * 0.5;
	}
#endif

	LOG_DEBUG << "Saving the signal...";
	project_.saveSignalToHDF5(signal, config_.savedAcqDir,
					0, config_.baseElement,
					config_.txGroupElement, config_.rxGroupElement);

	project_.showFigure2D(0, "A-scan", t, signal);

#if 0
	LOG_DEBUG << "Executing the plot script...";
	std::vector<std::string> args;
	std::ostringstream arg1;
	arg1 << config_.savedAcqDir << '/' << config_.signalFile << ".h5";
	args.push_back(arg1.str());
	std::ostringstream arg2;
	arg2.precision(FLOAT_SCIENTIFIC_NOTATION_NUM_DIGITS_AFTER_DECIMAL_POINT);
	arg2 << std::scientific << config_.samplingFrequency;
	args.push_back(arg2.str());
	project_.executeProgram(config_.plotScript, args);
#endif
}

} // namespace Lab
