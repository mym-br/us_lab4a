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

#ifndef SIMULATEDSTAACQUISITION_H_
#define SIMULATEDSTAACQUISITION_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "Log.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SimulatedAcquisitionDevice.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Util.h"
#include "XZ.h"



namespace Lab {

template<typename FloatType>
class SimulatedSTAAcquisition : public STAAcquisition<FloatType> {
public:
	SimulatedSTAAcquisition(Project& project, const STAConfiguration<FloatType>& config);
	virtual ~SimulatedSTAAcquisition();

	virtual void execute(unsigned int baseElement, unsigned int txElement,
				typename STAAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	SimulatedSTAAcquisition(const SimulatedSTAAcquisition&);
	SimulatedSTAAcquisition& operator=(const SimulatedSTAAcquisition&);

	Project& project_;
	const STAConfiguration<FloatType>& config_;
	FloatType maxAbsValue_;
	boost::scoped_ptr<SimulatedAcquisitionDevice<FloatType> > acqDevice_;
	std::vector<XZ<FloatType> > reflectorList_;
};



template<typename FloatType>
SimulatedSTAAcquisition<FloatType>::SimulatedSTAAcquisition(Project& project, const STAConfiguration<FloatType>& config)
		: project_(project)
		, config_(config)
		, maxAbsValue_(0.0)
		, acqDevice_(0)
{
	//TODO: check numChannels/numChannelsMux and other params

	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr pm = project_.loadChildParameterMap(taskPM, "simulated_sta_acquisition_config_file");
	const FloatType widthElem            = pm->value<FloatType>(   "element_width"      , 1.0e-6,  10.0e-3);
	const FloatType heightElem           = pm->value<FloatType>(   "element_height"     , 1.0e-3, 100.0e-3);
	const unsigned int numDivWidth       = pm->value<unsigned int>("num_width_div"      ,      1,     1000);
	const std::string reflectorsFileName = pm->value<std::string>( "reflectors_file_name");
	const FloatType reflectorsXOffset    = pm->value<FloatType>(   "reflectors_x_offset",   -0.2,      0.2);
	const FloatType noiseAmplitude       = pm->value<FloatType>(   "noise_amplitude"    ,    0.0,   1.0e10);

	std::vector<std::pair<FloatType, FloatType> > inputReflectorList;
	project_.loadHDF5(reflectorsFileName, "reflectors", inputReflectorList);

	reflectorList_.resize(inputReflectorList.size());
	for (std::size_t i = 0, end = inputReflectorList.size(); i < end; ++i) {
		std::pair<FloatType, FloatType>& orig = inputReflectorList[i];
		XZ<FloatType>& dest = reflectorList_[i];
		dest.x = orig.first + reflectorsXOffset;
		dest.z = orig.second;
	}

	acqDevice_.reset(new SimulatedAcquisitionDevice<FloatType>(
							config_.numElementsMux,
							config_.pitch,
							config_.samplingFrequency,
							widthElem,
							heightElem,
							numDivWidth,
							noiseAmplitude));
	acqDevice_->setAcquisitionTime(config_.acquisitionTime);
	acqDevice_->setExcitationWaveform(config_.centerFrequency);
	acqDevice_->setPropagationSpeed(config_.propagationSpeed);
	acqDevice_->setReflectorList(reflectorList_);
}

template<typename FloatType>
SimulatedSTAAcquisition<FloatType>::~SimulatedSTAAcquisition()
{
}

template<typename FloatType>
void
SimulatedSTAAcquisition<FloatType>::execute(unsigned int baseElement, unsigned int txElement,
						typename STAAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement << " txElement=" << txElement;

	const std::size_t ascanLength = acqDevice_->ascanLength();
	if (ascanLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "ascanLength = 0.");
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != ascanLength) {
		acqData.resize(config_.numElements, ascanLength);
		LOG_DEBUG << "RESIZE acqData(channels,ascanLength): channels=" << config_.numElements << " ascanLength=" << ascanLength;
	}

	std::vector<bool> txMask(config_.numElementsMux);
	txMask[baseElement + txElement] = true;
	acqDevice_->setActiveTxElements(txMask);

	std::vector<bool> rxMask(config_.numElementsMux);
	for (unsigned int i = baseElement, end = baseElement + config_.numElements; i < end; ++i) {
		rxMask[i] = true;
	}
	acqDevice_->setActiveRxElements(rxMask);

	const std::vector<FloatType>& ascanList = acqDevice_->getAscan();

	//TODO: copy to SimulatedComplSeqSTAAcquisition
	const FloatType mx = Util::maxAbsolute(ascanList);
	if (mx > maxAbsValue_) maxAbsValue_ = mx;
	LOG_INFO << "########## max(abs(acqData)) = " << mx << " global: " << maxAbsValue_;

	std::copy(ascanList.begin(), ascanList.end(), acqData.begin()); //TODO: memcpy? how to avoid this copy?
}

} // namespace Lab

#endif /* SIMULATEDSTAACQUISITION_H_ */
