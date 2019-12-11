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
#ifndef SIMULATED3DSTAACQUISITION_H
#define SIMULATED3DSTAACQUISITION_H

#include <memory>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Simulated3DAcquisitionDevice.h"
#include "SA3DConfiguration.h"
#include "STAAcquisition.h"
#include "Util.h"
#include "XYZValue.h"



namespace Lab {

template<typename FloatType>
class Simulated3DSTAAcquisition : public STAAcquisition<FloatType> {
public:
	Simulated3DSTAAcquisition(Project& project, const SA3DConfiguration<FloatType>& config);
	virtual ~Simulated3DSTAAcquisition() = default;

	virtual void prepare(unsigned int baseElement);
	virtual void execute(unsigned int txElement,
				typename STAAcquisition<FloatType>::AcquisitionDataType& acqData);

	void modifyReflectorsOffset(FloatType offsetX, FloatType offsetY);
private:
	Simulated3DSTAAcquisition(const Simulated3DSTAAcquisition&) = delete;
	Simulated3DSTAAcquisition& operator=(const Simulated3DSTAAcquisition&) = delete;
	Simulated3DSTAAcquisition(Simulated3DSTAAcquisition&&) = delete;
	Simulated3DSTAAcquisition& operator=(Simulated3DSTAAcquisition&&) = delete;

	Project& project_;
	const SA3DConfiguration<FloatType>& config_;
	FloatType maxAbsValue_; // auxiliar
	std::unique_ptr<Simulated3DAcquisitionDevice<FloatType>> acqDevice_;
	std::vector<XYZValue<FloatType>> reflectorList_;
	FloatType reflectorsOffsetX_;
	FloatType reflectorsOffsetY_;
	unsigned int baseElement_;
};



template<typename FloatType>
Simulated3DSTAAcquisition<FloatType>::Simulated3DSTAAcquisition(Project& project, const SA3DConfiguration<FloatType>& config)
		: project_(project)
		, config_(config)
		, maxAbsValue_()
		, baseElement_()
{
	const ParamMapPtr pm      = project_.loadChildParameterMap("simulated_3d_acquisition_config_file");
	const ParamMapPtr arrayPM = project_.loadChildParameterMap("array_config_file");

	const auto reflectorsFileName = pm->value<std::string>("reflectors_file");
	pm->getValue(reflectorsOffsetX_, "reflectors_offset_x", -10000.0, 10000.0);
	pm->getValue(reflectorsOffsetY_, "reflectors_offset_y", -10000.0, 10000.0);

	Matrix<FloatType> inputReflectorList;
	project_.loadHDF5(reflectorsFileName, "reflectors", inputReflectorList);
	if (inputReflectorList.n2() != 4) {
		THROW_EXCEPTION(InvalidValueException, "Wrong number of columns (" << inputReflectorList.n2() <<
				") in the file \"" << reflectorsFileName << ".h5\" (should be 4). ");
	}

	reflectorList_.resize(inputReflectorList.n1());
	for (std::size_t i = 0, end = inputReflectorList.n1(); i < end; ++i) {
		XYZValue<FloatType>& data = reflectorList_[i];
		data.x     = inputReflectorList(i, 0);
		data.y     = inputReflectorList(i, 1);
		data.z     = inputReflectorList(i, 2);
		data.value = inputReflectorList(i, 3);
	}

	acqDevice_ = std::make_unique<Simulated3DAcquisitionDevice<FloatType>>(
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

template<typename FloatType>
void
Simulated3DSTAAcquisition<FloatType>::prepare(unsigned int baseElement)
{
	std::vector<bool> rxMask(config_.rxElemPos.size());
	for (unsigned int localElem : config_.activeRxElem) {
		const unsigned int elem = baseElement + localElem;
		if (elem >= rxMask.size()) {
			THROW_EXCEPTION(InvalidValueException, "Invalid active rx element: " << elem << '.');
		}
		rxMask[elem] = true;
	}
	acqDevice_->setActiveRxElements(rxMask);

	baseElement_ = baseElement;
}

template<typename FloatType>
void
Simulated3DSTAAcquisition<FloatType>::execute(unsigned int txElement,
						typename STAAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement_ << " txElement=" << txElement;

	const std::size_t signalLength = acqDevice_->signalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "signalLength = 0.");
	}
	if (acqData.n1() != config_.activeRxElem.size() || acqData.n2() != signalLength) {
		acqData.resize(config_.activeRxElem.size(), signalLength);
		LOG_DEBUG << "RESIZE acqData(channels, signalLength): channels=" << config_.activeRxElem.size()
				<< " signalLength=" << signalLength;
	}

	std::vector<bool> txMask(config_.txElemPos.size());
	bool found = false;
	for (unsigned int localElem : config_.activeTxElem) {
		if (txElement == localElem) {
			found = true;
			const unsigned int elem = baseElement_ + localElem;
			if (elem >= txMask.size()) {
				THROW_EXCEPTION(InvalidValueException, "Invalid active tx element: " << elem << '.');
			}
			break;
		}
	}
	if (!found) THROW_EXCEPTION(InvalidValueException, "Invalid tx element: " << baseElement_ + txElement << '.');
	txMask[baseElement_ + txElement] = true;
	acqDevice_->setActiveTxElements(txMask);

	const std::vector<FloatType>& signalList = acqDevice_->getSignalList();

	const FloatType mx = Util::maxAbsolute(signalList);
	if (mx > maxAbsValue_) maxAbsValue_ = mx;
	LOG_DEBUG << "########## max(abs(signalList)): " << mx << " global: " << maxAbsValue_;

	std::copy(signalList.begin(), signalList.end(), acqData.begin());
}

template<typename FloatType>
void
Simulated3DSTAAcquisition<FloatType>::modifyReflectorsOffset(FloatType offsetX, FloatType offsetY)
{
	acqDevice_->setReflectorOffset(reflectorsOffsetX_ + offsetX, reflectorsOffsetY_ + offsetY);
}

} // namespace Lab

#endif // SIMULATED3DSTAACQUISITION_H
