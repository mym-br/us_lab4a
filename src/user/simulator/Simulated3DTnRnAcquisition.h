/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

#ifndef SIMULATED3DTNRNACQUISITION_H
#define SIMULATED3DTNRNACQUISITION_H

#include <memory>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Simulated3DAcquisitionDevice.h"
#include "TnRnAcquisition.h"
#include "TnRnConfiguration.h"
#include "Util.h"
#include "XYZValue.h"



namespace Lab {

template<typename TFloat>
class Simulated3DTnRnAcquisition : public TnRnAcquisition<TFloat> {
public:
	Simulated3DTnRnAcquisition(Project& project, const TnRnConfiguration<TFloat>& config);
	virtual ~Simulated3DTnRnAcquisition() = default;

	virtual void prepare(unsigned int baseElement, const std::vector<TFloat>& txDelays);
	virtual void execute(typename TnRnAcquisition<TFloat>::AcquisitionDataType& acqData);

	void modifyReflectorsOffset(TFloat offsetX, TFloat offsetY);
private:
	Simulated3DTnRnAcquisition(const Simulated3DTnRnAcquisition&) = delete;
	Simulated3DTnRnAcquisition& operator=(const Simulated3DTnRnAcquisition&) = delete;
	Simulated3DTnRnAcquisition(Simulated3DTnRnAcquisition&&) = delete;
	Simulated3DTnRnAcquisition& operator=(Simulated3DTnRnAcquisition&&) = delete;

	Project& project_;
	const TnRnConfiguration<TFloat>& config_;
	TFloat maxAbsValue_; // auxiliar
	std::unique_ptr<Simulated3DAcquisitionDevice<TFloat>> acqDevice_;
	std::vector<XYZValue<TFloat>> reflectorList_;
	TFloat reflectorsOffsetX_;
	TFloat reflectorsOffsetY_;
};



template<typename TFloat>
Simulated3DTnRnAcquisition<TFloat>::Simulated3DTnRnAcquisition(Project& project, const TnRnConfiguration<TFloat>& config)
		: project_(project)
		, config_(config)
		, maxAbsValue_()
{
	const ParamMapPtr pm      = project_.getSubParamMap("simulated_3d_acquisition_config_file");
	const ParamMapPtr arrayPM = project_.getSubParamMap("array_config_file");

	const auto reflectorsFileName = pm->value<std::string>("reflectors_file");
	pm->getValue(reflectorsOffsetX_, "reflectors_offset_x", -10000.0, 10000.0);
	pm->getValue(reflectorsOffsetY_, "reflectors_offset_y", -10000.0, 10000.0);

	Matrix<TFloat> inputReflectorList;
	project_.loadHDF5(reflectorsFileName, "reflectors", inputReflectorList);
	if (inputReflectorList.n2() != 4) {
		THROW_EXCEPTION(InvalidValueException, "Wrong number of columns (" << inputReflectorList.n2() <<
				") in the file \"" << reflectorsFileName << ".h5\" (should be 4). ");
	}

	reflectorList_.resize(inputReflectorList.n1());
	for (std::size_t i = 0, end = inputReflectorList.n1(); i < end; ++i) {
		XYZValue<TFloat>& data = reflectorList_[i];
		data.x     = inputReflectorList(i, 0);
		data.y     = inputReflectorList(i, 1);
		data.z     = inputReflectorList(i, 2);
		data.value = inputReflectorList(i, 3);
	}

	acqDevice_ = std::make_unique<Simulated3DAcquisitionDevice<TFloat>>(
								*pm,
								*arrayPM,
								config_.samplingFrequency,
								config_.propagationSpeed,
								config_.maxFrequency,
								project_.expDirectory());
	acqDevice_->setAcquisitionTime(config_.acquisitionTime);
	acqDevice_->setExcitationWaveform(config_.centerFrequency);
	acqDevice_->setReflectorList(reflectorList_);
	acqDevice_->setReflectorOffset(reflectorsOffsetX_, reflectorsOffsetY_);
	acqDevice_->setGain(config_.minGain);
}

template<typename TFloat>
void
Simulated3DTnRnAcquisition<TFloat>::prepare(unsigned int baseElement, const std::vector<TFloat>& txDelays)
{
	if (baseElement + config_.numElements > config_.numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid base element: " << baseElement << '.');
	}

	std::vector<bool> txMask(config_.txElemPos.size());
	for (unsigned int localElem = 0; localElem < config_.numElements; ++localElem) {
		txMask[baseElement + localElem] = true;
	}
	acqDevice_->setActiveTxElements(txMask);

	acqDevice_->setTxDelays(txDelays);

	std::vector<bool> rxMask(config_.rxElemPos.size());
	for (unsigned int localElem = 0; localElem < config_.numElements; ++localElem) {
		rxMask[baseElement + localElem] = true;
	}
	acqDevice_->setActiveRxElements(rxMask);
}

template<typename TFloat>
void
Simulated3DTnRnAcquisition<TFloat>::execute(typename TnRnAcquisition<TFloat>::AcquisitionDataType& acqData)
{
	const std::size_t signalLength = acqDevice_->signalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "signalLength = 0.");
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != signalLength) {
		acqData.resize(config_.numElements, signalLength);
		LOG_DEBUG << "RESIZE acqData(channels, signalLength): channels=" << config_.numElements
				<< " signalLength=" << signalLength;
	}

	const std::vector<TFloat>& signalList = acqDevice_->getSignalList();

	const TFloat mx = Util::maxAbsolute(signalList);
	if (mx > maxAbsValue_) maxAbsValue_ = mx;
	LOG_DEBUG << "########## max(abs(signalList)): " << mx << " global: " << maxAbsValue_;

	std::copy(signalList.begin(), signalList.end(), acqData.begin());
}

template<typename TFloat>
void
Simulated3DTnRnAcquisition<TFloat>::modifyReflectorsOffset(TFloat offsetX, TFloat offsetY)
{
	acqDevice_->setReflectorOffset(reflectorsOffsetX_ + offsetX, reflectorsOffsetY_ + offsetY);
}

} // namespace Lab

#endif // SIMULATED3DTNRNACQUISITION_H
