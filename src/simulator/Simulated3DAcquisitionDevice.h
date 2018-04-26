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
#ifndef SIMULATED3DACQUISITIONDEVICE_H
#define SIMULATED3DACQUISITIONDEVICE_H

#include <algorithm> /* copy, fill */
#include <complex>
#include <cstddef> /* std::size_t */
#include <memory>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayUtil.h"
#include "Decimator.h"
#include "Exception.h"
#include "FFTWFilter2.h"
#include "Log.h"
#include "MinstdPseudorandomNumberGenerator.h"
#include "NumericArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "ParameterMap.h"
#include "Util.h"
#include "Waveform.h"
#include "XY.h"
#include "XYZValue.h"

#define SIMULATED_3D_ACQUISITION_DEVICE_PSEUDORANDOM_NUMBER_GENERATOR_SEED 1

// Fraction of the destination bandwidth.
#define SIMULATED_3D_ACQUISITION_DEVICE_DECIMATOR_LP_FILTER_TRANSITION_WIDTH (0.3)

#define SIMULATED_3D_ACQUISITION_DEVICE_MAX_SIGNAL_LENGTH (std::numeric_limits<unsigned int>::max() / sizeof(FloatType))



//TODO: Clean.
namespace Lab {

class Project;

// x = 0, y = 0 is at the center of the array.
// z = 0 is at the surface of the array.
template<typename FloatType>
class Simulated3DAcquisitionDevice {
public:
	struct ThreadData {
		ThreadData(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType rxElemWidth,
			FloatType rxElemHeight,
			FloatType rxSubElemSize,
			const Decimator<FloatType>& dec)
				: rxImpResp{samplingFreq, propagationSpeed, rxElemWidth, rxElemHeight, rxSubElemSize}
				, decimator{dec}
		{
		}
		NumericRectangularFlatSourceImpulseResponse<FloatType> rxImpResp;
		std::vector<FloatType> hRx;
		std::vector<FloatType> p;
		std::vector<FloatType> pDown;
		FFTWFilter2<FloatType> convDadtHTxFilter;
		std::vector<std::complex<FloatType>> dadtFilterFreqCoeff;
		Decimator<FloatType> decimator;
	};

	Simulated3DAcquisitionDevice(const ParameterMap& pm, FloatType outputSamplingFreq,
					FloatType propagationSpeed, FloatType maxFrequency);
	~Simulated3DAcquisitionDevice();

	void setAcquisitionTime(FloatType acqTime);

	std::size_t signalLength() const { return signalLength_; }
	const std::vector<FloatType>& getSignalList();

	void setActiveTxElements(const std::vector<bool>& mask);
	void setActiveRxElements(const std::vector<bool>& mask);

	// Must be preceded by a call to setActiveTxElements(const std::vector<bool>& mask).
	void setTxDelays(const std::vector<FloatType>& delays);

	// Must be preceded by a call to setActiveTxElements(const std::vector<bool>& mask).
	void setTxFocalPoint(FloatType xf, FloatType yf, FloatType zf);

	void setReflectorList(const std::vector<XYZValue<FloatType>>& reflectorList) {
		reflectorList_ = reflectorList;
	}
	void setReflectorOffset(FloatType reflectorsOffsetX, FloatType reflectorsOffsetY) {
		reflectorsOffsetX_ = reflectorsOffsetX;
		reflectorsOffsetY_ = reflectorsOffsetY;
	}

	void setExcitationWaveform(FloatType centerFrequency /* Hz */);
	void setExcitationWaveform(const std::vector<FloatType>& dadtExc);

private:
	Simulated3DAcquisitionDevice(const Simulated3DAcquisitionDevice&) = delete;
	Simulated3DAcquisitionDevice& operator=(const Simulated3DAcquisitionDevice&) = delete;

	void prepareExcitationDadt(const std::vector<FloatType>& vExc);

	const FloatType c_; // propagation speed (m/s)
	const FloatType invC_;
	FloatType simFs_; // simulation sampling frequency (Hz)
	FloatType outFs_; // output sampling frequency (Hz)
	FloatType elemWidth_; // width of each element (m)
	FloatType elemHeight_; // height of each element (m)
	FloatType subElemSize_; // target size of sub-elements (m)
	FloatType noiseAmplitude_;
	std::string excitationType_;
	FloatType excNumPeriods_;
	std::vector<FloatType> dadt_;
	std::vector<FloatType> hTx_;
	std::vector<unsigned int> activeTxElem_;
	std::vector<unsigned int> activeRxElem_;
	std::size_t signalLength_;
	std::vector<FloatType> signalList_;
	std::vector<XYZValue<FloatType>> reflectorList_; // coordinates of the reflectors
	FloatType reflectorsOffsetX_;
	FloatType reflectorsOffsetY_;
	MinstdPseudorandomNumberGenerator prng_;
	std::vector<std::complex<FloatType>> dadtFilterFreqCoeff_;
	std::vector<FloatType> convDadtHTx_;
	std::vector<std::complex<FloatType>> convDadtHTxFilterFreqCoeff_;
	std::vector<XY<FloatType>> elemPos_;
	std::vector<FloatType> txDelays_;
	std::unique_ptr<NumericArrayOfRectangularFlatSourcesImpulseResponse<FloatType>> txImpResp_;
	std::unique_ptr<Decimator<FloatType>> decimator_;
};



template<typename FloatType>
Simulated3DAcquisitionDevice<FloatType>::Simulated3DAcquisitionDevice(
		const ParameterMap& pm,
		FloatType outputSamplingFreq,
		FloatType propagationSpeed,
		FloatType maxFrequency)
			: c_{propagationSpeed}
			, invC_{1 / c_}
			, simFs_{}
			, outFs_{}
			, elemWidth_{}
			, elemHeight_{}
			, subElemSize_{}
			, noiseAmplitude_{}
			, excNumPeriods_{}
			, signalLength_{}
			, reflectorsOffsetX_{}
			, reflectorsOffsetY_{}
			, prng_{SIMULATED_3D_ACQUISITION_DEVICE_PSEUDORANDOM_NUMBER_GENERATOR_SEED}
{
	elemWidth_      = pm.value<FloatType>(  "element_width"         , 1.0e-6, 1000.0e-3);
	elemHeight_     = pm.value<FloatType>(  "element_height"        , 1.0e-6, 1000.0e-3);
	noiseAmplitude_ = pm.value<FloatType>(  "noise_amplitude"       ,    0.0,   1.0e100);
	excitationType_ = pm.value<std::string>("excitation_type");
	excNumPeriods_  = pm.value<FloatType>(  "excitation_num_periods",    0.0,     100.0);
	const FloatType nyquistRate = 2.0 * maxFrequency;
	subElemSize_ = propagationSpeed / (nyquistRate * pm.value<FloatType>("sub_elem_size_factor", 0.0, 1000.0));
	simFs_ = outputSamplingFreq * pm.value<unsigned int>("sim_sampling_frequency_factor", 1, 10000);
	outFs_ = outputSamplingFreq;

	if (subElemSize_ == 0.0) {
		THROW_EXCEPTION(InvalidParameterException, "The size of sub-elements is equal to zero.");
	}

	// Calculate the coordinates of the centers of the elements.
	ArrayUtil::calculateXYArrayParameters(pm, propagationSpeed, simFs_, elemPos_, txDelays_, false);

	txImpResp_ = std::make_unique<NumericArrayOfRectangularFlatSourcesImpulseResponse<FloatType>>(
				simFs_, propagationSpeed, elemWidth_, elemHeight_, subElemSize_,
				elemPos_, txDelays_);

	decimator_ = std::make_unique<Decimator<FloatType>>();
	decimator_->prepare(simFs_ / outFs_, SIMULATED_3D_ACQUISITION_DEVICE_DECIMATOR_LP_FILTER_TRANSITION_WIDTH);
}

template<typename FloatType>
Simulated3DAcquisitionDevice<FloatType>::~Simulated3DAcquisitionDevice()
{
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setAcquisitionTime(FloatType acqTime)
{
	FloatType len = std::ceil(acqTime * outFs_);
	if (len > SIMULATED_3D_ACQUISITION_DEVICE_MAX_SIGNAL_LENGTH) {
		THROW_EXCEPTION(InvalidValueException, "The acquisition time is too long.");
	}

	signalLength_ = static_cast<std::size_t>(len);
}

// Convert the excitation velocity to da/dt.
template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::prepareExcitationDadt(const std::vector<FloatType>& vExc)
{
	const FloatType period = 1.0 / simFs_;

	// Calculates the excitation acceleration.
	std::vector<FloatType> aExc;
	Util::centralDiff(vExc, period, aExc); // adds delay of 1 sample

	// Calculates the excitation d(acceleration)/dt.
	std::vector<FloatType> unfilteredDadt;
	Util::centralDiff(aExc, period, unfilteredDadt); // adds delay of 1 sample
	std::vector<std::complex<FloatType>> filterFreqCoeff;
	FFTWFilter2<FloatType> f;
	f.setCoefficients(decimator_->lowPassFIRFilter(), filterFreqCoeff);
	f.filter(filterFreqCoeff, unfilteredDadt, dadt_); // adds delay of (filter size - 1) / 2 samples

	// Normalize.
	const FloatType coeff = 1 / Util::maxAbsolute(dadt_);
	Util::multiply(dadt_, coeff);
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setExcitationWaveform(FloatType centerFrequency)
{
	std::vector<FloatType> vExc;
	if (excitationType_ == "1") {
		Waveform::getType1(centerFrequency, simFs_, vExc, excNumPeriods_);
	} else if (excitationType_ == "2a") {
		Waveform::getType2a(centerFrequency, simFs_, vExc, excNumPeriods_);
	} else if (excitationType_ == "2b") {
		Waveform::getType2b(centerFrequency, simFs_, vExc, excNumPeriods_);
	} else if (excitationType_ == "2c") {

	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType_ << '.');
	}

	prepareExcitationDadt(vExc);
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setExcitationWaveform(const std::vector<FloatType>& dadtExc)
{
	dadt_ = dadtExc;
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setActiveTxElements(const std::vector<bool>& mask)
{
	if (mask.size() != elemPos_.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong size of the tx mask: " << mask.size()
				<< " (should be " << elemPos_.size() << ").");
	}

	activeTxElem_.clear();
	for (unsigned int i = 0, end = mask.size(); i < end; ++i) {
		if (mask[i]) {
			activeTxElem_.push_back(i);
		}
	}
	if (activeTxElem_.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "There is no active transmit element.");
	}
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setActiveRxElements(const std::vector<bool>& mask)
{
	if (signalLength_ == 0) {
		THROW_EXCEPTION(InvalidStateException, "The acquisition time has not been set.");
	}
	if (mask.size() != elemPos_.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong size of the rx mask: " << mask.size()
				<< " (should be " << elemPos_.size() << ").");
	}

	activeRxElem_.clear();
	for (unsigned int i = 0, end = mask.size(); i < end; ++i) {
		if (mask[i]) {
			activeRxElem_.push_back(i);
		}
	}
	if (activeRxElem_.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "There is no active receive element.");
	}

	signalList_.resize(activeRxElem_.size() * signalLength_);
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setTxDelays(const std::vector<FloatType>& delays)
{
	if (delays.size() != txDelays_.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong size of the delay vector: " << delays.size()
				<< " (should be " << txDelays_.size() << ").");
	}
	std::copy(delays.begin(), delays.end(), txDelays_.begin());
}

template<typename FloatType>
void
Simulated3DAcquisitionDevice<FloatType>::setTxFocalPoint(FloatType xf, FloatType yf, FloatType zf)
{
	if (c_ == 0.0) {
		THROW_EXCEPTION(InvalidStateException, "The propagation speed has not been set.");
	}

	FloatType maxDt = 0.0;

	// For each active transmit element:
	for (auto iter = activeTxElem_.begin(); iter != activeTxElem_.end(); ++iter) {
		const FloatType dx = xf - elemPos_[*iter].x;
		const FloatType dy = yf - elemPos_[*iter].y;
		const FloatType r = std::sqrt(dx * dx + dy * dy + zf * zf);
		const FloatType dt = r * invC_; // transit time from the element to the focal point
		if (dt > maxDt) maxDt = dt;
		txDelays_[*iter] = dt;
	}

	// For each active transmit element:
	for (auto iter = activeTxElem_.begin(); iter != activeTxElem_.end(); ++iter) {
		txDelays_[*iter] = maxDt - txDelays_[*iter];
	}
}

template<typename FloatType>
const std::vector<FloatType>&
Simulated3DAcquisitionDevice<FloatType>::getSignalList()
{
	if (signalList_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The active receive elements have not been defined.");
	}
	if (dadt_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The excitation waveform has not been set.");
	}
	if (c_ == 0.0) {
		THROW_EXCEPTION(InvalidStateException, "The propagation speed has not been set.");
	}
	if (reflectorList_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The list of reflectors is empty.");
	}

	std::fill(signalList_.begin(), signalList_.end(), FloatType{0});
	FFTWFilter2<FloatType> dadtFilter;
	dadtFilter.setCoefficients(dadt_, dadtFilterFreqCoeff_);

	// For each reflector:
	for (std::size_t iRef = 0, iRefEnd = reflectorList_.size(); iRef < iRefEnd; ++iRef) {
		LOG_INFO << "ACQ Reflector: " << iRef << " < " << iRefEnd;

		const XYZValue<FloatType>& reflector = reflectorList_[iRef];
		const FloatType refX = reflector.x;
		const FloatType refY = reflector.y;
		const FloatType refZ = reflector.z;
		const FloatType refCoeff = reflector.value;

		// Calculate the impulse response in transmission (all active elements).
		std::size_t hTxOffset;
		txImpResp_->getImpulseResponse(refX, refY, refZ, hTxOffset, hTx_, &activeTxElem_);//TODO: multi-threading

		// dadt * hTx
		dadtFilter.filter(dadtFilterFreqCoeff_, hTx_, convDadtHTx_);

		ThreadData threadData{
			simFs_,
			c_,
			elemWidth_,
			elemHeight_,
			subElemSize_,
			*decimator_
		};
		threadData.convDadtHTxFilter.setCoefficients(convDadtHTx_, convDadtHTxFilterFreqCoeff_);
		tbb::enumerable_thread_specific<ThreadData> tls{threadData};

		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, activeRxElem_.size()),
		[&, refX, refY, refZ, hTxOffset](const tbb::blocked_range<std::size_t>& r) {
			auto& local = tls.local();

			// For each active receive element:
			for (std::size_t iActiveRx = r.begin(); iActiveRx != r.end(); ++iActiveRx) {
				const unsigned int activeRxElem = activeRxElem_[iActiveRx];
				const XY<FloatType>& pos = elemPos_[activeRxElem];

				// Calculate the impulse response in reception (only for the active element).
				std::size_t hRxOffset;
				local.rxImpResp.getImpulseResponse(refX - pos.x, refY - pos.y, refZ, hRxOffset, local.hRx);

				local.convDadtHTxFilter.filter(convDadtHTxFilterFreqCoeff_, local.hRx, local.p);

				const std::size_t pOffset = hTxOffset + hRxOffset;
				// Offset due to the low-pass filter applied to the excitation waveform.
				const std::size_t lpOffset = (decimator_->lowPassFIRFilter().size() - 1) / 2;

				// Downsample to the output sampling frequency.
				std::size_t pDownOffset;
				decimator_->downsample(lpOffset, pOffset, local.p, pDownOffset, local.pDown);

				const std::size_t signalListOffset = iActiveRx * signalLength_;
				for (std::size_t i = 0, iEnd = local.pDown.size(), j = pDownOffset;
						(i < iEnd) && (j < signalLength_);
						++i, ++j) {
					signalList_[signalListOffset + j] += local.pDown[i] * refCoeff;
				}
			}
		});
	}

	if (noiseAmplitude_ != 0.0) {
		// Add noise.
		const FloatType a = 2.0 * noiseAmplitude_;
		for (auto& sample : signalList_) {
			sample += (prng_.get() - 0.5) * a;
		}
	}

	return signalList_;
}

} // namespace Lab

#endif // SIMULATED3DACQUISITIONDEVICE_H
