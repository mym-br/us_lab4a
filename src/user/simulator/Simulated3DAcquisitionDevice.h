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
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "AnalyticRectangularSourceImpulseResponse.h"
#include "ArrayOfRectangularSourcesImpulseResponse.h"
#include "ArrayUtil.h"
#include "Decimator.h"
#include "Exception.h"
#include "FFTWFilter2.h"
#include "Log.h"
#include "PseudorandomNumberGenerator.h"
#include "NumericRectangularSourceImpulseResponse.h"
#include "ParameterMap.h"
#include "Util.h"
#include "Waveform.h"
#include "WavefrontObjFileWriter.h"
#include "XY.h"
#include "XYZValue.h"



namespace Lab {

class Project;

// x = 0, y = 0 is at the center of the array.
// z = 0 is at the surface of the array.
template<typename TFloat>
class Simulated3DAcquisitionDevice {
public:
	template<typename ImpulseResponse>
	struct ThreadData {
		ThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat rxElemWidth,
			TFloat rxElemHeight,
			TFloat rxDiscretization)
				: rxImpResp(samplingFreq, propagationSpeed, rxElemWidth, rxElemHeight, rxDiscretization) {}
		ImpulseResponse rxImpResp;
		std::vector<TFloat> hRx;
		std::vector<TFloat> p;
		std::vector<TFloat> pDown;
		FFTWFilter2<TFloat> convDadtHTxFilter;
		std::vector<std::complex<TFloat>> dadtFilterFreqCoeff;
	};

	Simulated3DAcquisitionDevice(const ParameterMap& pm, const ParameterMap& arrayPM,
					TFloat outputSamplingFreq, TFloat propagationSpeed, TFloat maxFrequency,
					const std::string& expDirectory);
	~Simulated3DAcquisitionDevice() = default;

	void setAcquisitionTime(TFloat acqTime);

	std::size_t signalLength() const { return signalLength_; }
	const std::vector<TFloat>& getSignalList();

	void setActiveTxElements(const std::vector<bool>& mask);
	void setActiveRxElements(const std::vector<bool>& mask);

	// Must be preceded by a call to setActiveTxElements(const std::vector<bool>& mask).
	void setTxDelays(const std::vector<TFloat>& delays);

	// Must be preceded by a call to setActiveTxElements(const std::vector<bool>& mask).
	void setTxFocalPoint(TFloat xf, TFloat yf, TFloat zf);

	void setReflectorList(const std::vector<XYZValue<TFloat>>& reflectorList) {
		reflectorList_ = reflectorList;
	}
	void setReflectorOffset(TFloat reflectorsOffsetX, TFloat reflectorsOffsetY) {
		reflectorsOffsetX_ = reflectorsOffsetX;
		reflectorsOffsetY_ = reflectorsOffsetY;
	}

	void setExcitationWaveform(TFloat centerFrequency /* Hz */);
	void setExcitationWaveform(const std::vector<TFloat>& dadtExc);

	void setGain(TFloat gain /* dB */);

private:
	// Fraction of the destination bandwidth.
	static constexpr TFloat decimatorLPFilterTransitionWidth = 0.3;
	static constexpr TFloat maxSignalLength = std::numeric_limits<unsigned int>::max() / sizeof(TFloat);

	Simulated3DAcquisitionDevice(const Simulated3DAcquisitionDevice&) = delete;
	Simulated3DAcquisitionDevice& operator=(const Simulated3DAcquisitionDevice&) = delete;
	Simulated3DAcquisitionDevice(Simulated3DAcquisitionDevice&&) = delete;
	Simulated3DAcquisitionDevice& operator=(Simulated3DAcquisitionDevice&&) = delete;

	void prepareExcitationDadt(const std::vector<TFloat>& vExc);
	template<typename ImpulseResponse> void processReflector(const XYZValue<TFloat>& reflector,
									FFTWFilter2<TFloat>& dadtFilter);

	const TFloat c_; // propagation speed (m/s)
	const TFloat invC_;
	TFloat simFs_; // simulation sampling frequency (Hz)
	TFloat outFs_; // output sampling frequency (Hz)
	TFloat txElemWidth_; // width of each element (m)
	TFloat txElemHeight_; // height of each element (m)
	TFloat rxElemWidth_; // width of each element (m)
	TFloat rxElemHeight_; // height of each element (m)
	bool useNumericMethod_;
	TFloat txDiscretization_;
	TFloat rxDiscretization_;
	TFloat noiseAmplitude_;
	std::string excitationType_;
	TFloat excNumPeriods_;
	std::vector<TFloat> dadt_;
	std::vector<TFloat> hTx_;
	std::vector<unsigned int> activeTxElem_;
	std::vector<unsigned int> activeRxElem_;
	std::size_t signalLength_;
	std::vector<TFloat> signalList_;
	std::vector<XYZValue<TFloat>> reflectorList_; // coordinates of the reflectors
	TFloat reflectorsOffsetX_;
	TFloat reflectorsOffsetY_;
	TFloat signalCoeff_;
	PseudorandomNumberGenerator prng_;
	std::vector<std::complex<TFloat>> dadtFilterFreqCoeff_;
	std::vector<TFloat> convDadtHTx_;
	std::vector<std::complex<TFloat>> convDadtHTxFilterFreqCoeff_;
	std::vector<XY<TFloat>> txElemPos_;
	std::vector<XY<TFloat>> rxElemPos_;
	std::vector<TFloat> txDelays_;
	std::unique_ptr<Decimator<TFloat>> decimator_;
};



template<typename TFloat>
Simulated3DAcquisitionDevice<TFloat>::Simulated3DAcquisitionDevice(
		const ParameterMap& pm,
		const ParameterMap& arrayPM,
		TFloat outputSamplingFreq,
		TFloat propagationSpeed,
		TFloat maxFrequency,
		const std::string& expDirectory)
			: c_(propagationSpeed)
			, invC_(1 / c_)
			, simFs_()
			, outFs_()
			, txElemWidth_()
			, txElemHeight_()
			, rxElemWidth_()
			, rxElemHeight_()
			, useNumericMethod_()
			, txDiscretization_()
			, rxDiscretization_()
			, noiseAmplitude_()
			, excNumPeriods_()
			, signalLength_()
			, reflectorsOffsetX_()
			, reflectorsOffsetY_()
			, signalCoeff_(1.0)
{
	arrayPM.getValue(txElemWidth_ , "tx_element_width" , 1.0e-6, 1000.0e-3);
	arrayPM.getValue(txElemHeight_, "tx_element_height", 1.0e-6, 1000.0e-3);
	arrayPM.getValue(rxElemWidth_ , "rx_element_width" , 1.0e-6, 1000.0e-3);
	arrayPM.getValue(rxElemHeight_, "rx_element_height", 1.0e-6, 1000.0e-3);

	pm.getValue(noiseAmplitude_, "noise_amplitude"       , 0.0, 1.0e100);
	pm.getValue(excitationType_, "excitation_type");
	pm.getValue(excNumPeriods_ , "excitation_num_periods", 0.0, 100.0);

	const auto irMethod = pm.value<std::string>("impulse_response_method");
	const TFloat nyquistRate = Util::nyquistRate(maxFrequency);
	if (irMethod == "numeric") {
		useNumericMethod_ = true;
		const TFloat txSubElemSize = propagationSpeed /
						(nyquistRate * pm.value<TFloat>("tx_sub_elem_size_factor", 0.0, 1.0e3));
		const TFloat rxSubElemSize = propagationSpeed /
						(nyquistRate * pm.value<TFloat>("rx_sub_elem_size_factor", 0.0, 1.0e3));
		if (txSubElemSize == 0.0) {
			THROW_EXCEPTION(InvalidParameterException, "The size of the transmit sub-elements is equal to zero.");
		}
		if (rxSubElemSize == 0.0) {
			THROW_EXCEPTION(InvalidParameterException, "The size of the receive sub-elements is equal to zero.");
		}
		txDiscretization_ = txSubElemSize;
		rxDiscretization_ = rxSubElemSize;
	} else if (irMethod == "analytic") {
		const auto txMinEdgeDivisor = pm.value<TFloat>("tx_min_edge_divisor", 0.0, 1.0e6);
		const auto rxMinEdgeDivisor = pm.value<TFloat>("rx_min_edge_divisor", 0.0, 1.0e6);
		txDiscretization_ = txMinEdgeDivisor;
		rxDiscretization_ = rxMinEdgeDivisor;
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << irMethod << '.');
	}

	simFs_ = outputSamplingFreq * pm.value<unsigned int>("sim_sampling_frequency_factor", 1, 10000);
	outFs_ = outputSamplingFreq;

	// Calculate the coordinates of the centers of the elements.
	ArrayUtil::calculateTxElementPositions(arrayPM, txElemPos_);
	ArrayUtil::calculateRxElementPositions(arrayPM, rxElemPos_);

	decimator_ = std::make_unique<Decimator<TFloat>>();
	decimator_->prepare(simFs_ / outFs_, decimatorLPFilterTransitionWidth);

	{
		WavefrontObjFileWriter<TFloat> fw((expDirectory + "/tx_geometry.obj").c_str());
		const TFloat hw = 0.5 * txElemWidth_;
		const TFloat hh = 0.5 * txElemHeight_;
		for (const auto& pos : txElemPos_) {
			fw.addPoint(pos.x - hw, pos.y + hh, 0.0);
			fw.addPoint(pos.x + hw, pos.y + hh, 0.0);
			fw.addPoint(pos.x + hw, pos.y - hh, 0.0);
			fw.addPoint(pos.x - hw, pos.y - hh, 0.0);
			fw.addQuad(-4, -3, -2, -1);
		}
		fw.write();
	}
	{
		WavefrontObjFileWriter<TFloat> fw((expDirectory + "/rx_geometry.obj").c_str());
		const TFloat hw = 0.5 * rxElemWidth_;
		const TFloat hh = 0.5 * rxElemHeight_;
		for (const auto& pos : rxElemPos_) {
			fw.addPoint(pos.x - hw, pos.y + hh, 0.0);
			fw.addPoint(pos.x + hw, pos.y + hh, 0.0);
			fw.addPoint(pos.x + hw, pos.y - hh, 0.0);
			fw.addPoint(pos.x - hw, pos.y - hh, 0.0);
			fw.addQuad(-4, -3, -2, -1);
		}
		fw.write();
	}
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setAcquisitionTime(TFloat acqTime)
{
	TFloat len = std::ceil(acqTime * outFs_);
	if (len > maxSignalLength) {
		THROW_EXCEPTION(InvalidValueException, "The acquisition time is too long.");
	}

	signalLength_ = static_cast<std::size_t>(len);
}

// Convert the excitation velocity to da/dt.
template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::prepareExcitationDadt(const std::vector<TFloat>& vExc)
{
	const TFloat period = 1.0 / simFs_;

	// Calculate the excitation acceleration.
	std::vector<TFloat> aExc;
	Util::centralDiff(vExc, period, aExc); // adds delay of 1 sample

	// Calculate the excitation d(acceleration)/dt.
	std::vector<TFloat> unfilteredDadt;
	Util::centralDiff(aExc, period, unfilteredDadt); // adds delay of 1 sample
	std::vector<std::complex<TFloat>> filterFreqCoeff;
	FFTWFilter2<TFloat> f;
	f.setCoefficients(decimator_->lowPassFIRFilter(), filterFreqCoeff);
	f.filter(filterFreqCoeff, unfilteredDadt, dadt_); // adds delay of (filter size - 1) / 2 samples

	Util::normalizeBySumOfAbs(dadt_);
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setExcitationWaveform(TFloat centerFrequency)
{
	std::vector<TFloat> vExc;
	Waveform::get(excitationType_, centerFrequency, simFs_, excNumPeriods_, vExc);

	prepareExcitationDadt(vExc);
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setExcitationWaveform(const std::vector<TFloat>& dadtExc)
{
	dadt_ = dadtExc;
	Util::normalizeBySumOfAbs(dadt_);
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setGain(TFloat gain)
{
	signalCoeff_ = Util::decibelsToLinear(gain);
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setActiveTxElements(const std::vector<bool>& mask)
{
	if (mask.size() != txElemPos_.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong size of the tx mask: " << mask.size()
				<< " (should be " << txElemPos_.size() << ").");
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

	txDelays_.assign(activeTxElem_.size(), 0.0);
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setActiveRxElements(const std::vector<bool>& mask)
{
	if (signalLength_ == 0) {
		THROW_EXCEPTION(InvalidStateException, "The acquisition time has not been set.");
	}
	if (mask.size() != rxElemPos_.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong size of the rx mask: " << mask.size()
				<< " (should be " << rxElemPos_.size() << ").");
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

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setTxDelays(const std::vector<TFloat>& delays)
{
	if (delays.size() != txDelays_.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong size of the delay vector: " << delays.size()
				<< " (should be " << txDelays_.size() << ").");
	}
	std::copy(delays.begin(), delays.end(), txDelays_.begin());
}

template<typename TFloat>
void
Simulated3DAcquisitionDevice<TFloat>::setTxFocalPoint(TFloat xf, TFloat yf, TFloat zf)
{
	if (c_ == 0.0) {
		THROW_EXCEPTION(InvalidStateException, "The propagation speed has not been set.");
	}
	if (zf <= 0.0) {
		THROW_EXCEPTION(InvalidStateException, "Focus z must be > 0.");
	}

	TFloat maxDt = 0.0;
	txDelays_.resize(activeTxElem_.size());

	// For each active transmit element:
	for (unsigned int i = 0, iEnd = activeTxElem_.size(); i < iEnd; ++i) {
		const unsigned int txElem = activeTxElem_[i];
		const TFloat dt = Geometry::distance3DZ0(txElemPos_[txElem].x, txElemPos_[txElem].y,
								xf, yf, zf) * invC_;
		if (dt > maxDt) maxDt = dt;
		txDelays_[i] = dt;
	}

	// For each active transmit element:
	for (unsigned int i = 0, iEnd = activeTxElem_.size(); i < iEnd; ++i) {
		txDelays_[i] = maxDt - txDelays_[i];
	}
}

template<typename TFloat>
template<typename ImpulseResponse>
void
Simulated3DAcquisitionDevice<TFloat>::processReflector(const XYZValue<TFloat>& reflector, FFTWFilter2<TFloat>& dadtFilter)
{
	const TFloat refX = reflector.x + reflectorsOffsetX_;
	const TFloat refY = reflector.y + reflectorsOffsetY_;
	const TFloat refZ = reflector.z;
	const TFloat refCoeff = reflector.value;

	auto txImpResp = std::make_unique<ArrayOfRectangularSourcesImpulseResponse<TFloat, ImpulseResponse>>(
				simFs_, c_, txElemWidth_, txElemHeight_, txDiscretization_, txElemPos_, txDelays_);

	// Calculate the impulse response in transmission (all active elements).
	std::size_t hTxOffset;
	txImpResp->getImpulseResponse(refX, refY, refZ, hTxOffset, hTx_, &activeTxElem_);

	// dadt * hTx
	dadtFilter.filter(dadtFilterFreqCoeff_, hTx_, convDadtHTx_);
	Util::multiply(convDadtHTx_, 1 / simFs_);

	ThreadData<ImpulseResponse> threadData{
		simFs_,
		c_,
		rxElemWidth_,
		rxElemHeight_,
		rxDiscretization_
	};
	threadData.convDadtHTxFilter.setCoefficients(convDadtHTx_, convDadtHTxFilterFreqCoeff_);
	tbb::enumerable_thread_specific<ThreadData<ImpulseResponse>> tls(threadData);

	// If activeRxElem_.size() == 1, only one thread will be used.
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, activeRxElem_.size()),
	[&, refX, refY, refZ, hTxOffset](const tbb::blocked_range<std::size_t>& r) {
		auto& local = tls.local();

		// For each active receive element:
		for (std::size_t iActiveRx = r.begin(); iActiveRx != r.end(); ++iActiveRx) {
			const unsigned int activeRxElem = activeRxElem_[iActiveRx];
			const XY<TFloat>& pos = rxElemPos_[activeRxElem];

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

template<typename TFloat>
const std::vector<TFloat>&
Simulated3DAcquisitionDevice<TFloat>::getSignalList()
{
	if (activeTxElem_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The active transmit elements have not been defined.");
	}
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

	std::fill(signalList_.begin(), signalList_.end(), TFloat(0));
	FFTWFilter2<TFloat> dadtFilter;
	dadtFilter.setCoefficients(dadt_, dadtFilterFreqCoeff_);

	TFloat refCoeffSum = 0.0;

	// For each reflector:
	for (std::size_t iRef = 0, iRefEnd = reflectorList_.size(); iRef < iRefEnd; ++iRef) {
		LOG_DEBUG << "ACQ Reflector: " << iRef << " < " << iRefEnd;

		refCoeffSum += std::abs(reflectorList_[iRef].value);

		if (useNumericMethod_) {
			processReflector<NumericRectangularSourceImpulseResponse<TFloat>>(reflectorList_[iRef], dadtFilter);
		} else {
			processReflector<AnalyticRectangularSourceImpulseResponse<TFloat>>(reflectorList_[iRef], dadtFilter);
		}
	}

	Util::multiply(signalList_, signalCoeff_ / refCoeffSum);

	// Add noise.
	if (noiseAmplitude_ == 0.0) {
		noiseAmplitude_ = 1.0e-4; // -80 dB
	}
	const TFloat a = 2.0 * noiseAmplitude_;
	for (auto& sample : signalList_) {
		sample += (prng_.get() - 0.5) * a;
	}

	return signalList_;
}

} // namespace Lab

#endif // SIMULATED3DACQUISITIONDEVICE_H
