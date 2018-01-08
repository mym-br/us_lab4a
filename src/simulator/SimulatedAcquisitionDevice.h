#ifndef SIMULATEDACQUISITIONDEVICE_H_
#define SIMULATEDACQUISITIONDEVICE_H_

#include <algorithm>
#include <complex>
#include <cstddef> /* std::size_t */
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "Exception.h"
#include "FFTWFilter2.h"
#include "IndexValue.h"
#include "Log.h"
#include "Matrix2.h"
#include "MinstdPseudorandomNumberGenerator.h"
#include "Util.h"
#include "Waveform.h"
#include "XY.h"
#include "XZ.h"

#define SIMULATED_ACQUISITION_DEVICE_PSEUDORANDOM_NUMBER_GENERATOR_SEED 1



namespace Lab {

class Project;

// x = 0 is at the center of the first element.
// y = 0
// z = 0 is at the surface of the array.
template<typename FloatType>
class SimulatedAcquisitionDevice {
public:
	struct ThreadData {
		std::vector<FloatType> convHtHr;
		std::vector<FloatType> hr;
		std::vector<FloatType> p;
		std::vector<IndexValue<unsigned int, FloatType> > idxCoeffList;
		FFTWFilter2<FloatType> filter;
	};

	SimulatedAcquisitionDevice(
		unsigned int numElements,
		FloatType pitch,
		FloatType samplingFrequency,
		FloatType widthElem,
		FloatType heightElem,
		unsigned int numDivWidth,
		FloatType noiseAmplitude);

	~SimulatedAcquisitionDevice();

	FloatType upsampledSamplingFrequency() const { return fsUp_; }

	void setAcquisitionTime(FloatType acqTime);
	void setPropagationSpeed(FloatType c) {
		c_ = c;
		invC_ = 1.0 / c;
	}

	const std::vector<FloatType>& getAscan();
	std::size_t ascanLength() const { return ascanLength_; }

	void setActiveRxElements(const std::vector<bool>& mask);
	void setActiveTxElements(const std::vector<bool>& mask);

	// Must be preceded by a call to setActiveTxElements(const std::vector<bool>& mask).
	void setTxDelays(const std::vector<FloatType>& delays);

	// Must be preceded by a call to setActiveTxElements(const std::vector<bool>& mask).
	void setTxFocalPoint(FloatType xf, FloatType zf);

	void setReflectorList(const std::vector<XZ<FloatType> >& reflectorList) { reflectorList_ = reflectorList; }

	void prepareExcitationDadt(const std::vector<FloatType>& vExc, std::vector<FloatType>& dadtExc);
	void setExcitationWaveform(FloatType centerFrequency /* Hz */);
	void setExcitationWaveform(const std::vector<FloatType>& dadtExc);

private:
	enum {
		USE_LOWPASS_FIR_FILTER = 1, // 0 (disabled) or 1 (enabled)
		UPSAMPLING_FACTOR = 4       // 4 or 8
	};

	class CalculateAscan;

	SimulatedAcquisitionDevice(const SimulatedAcquisitionDevice&);
	SimulatedAcquisitionDevice& operator=(const SimulatedAcquisitionDevice&);

	void loadFilter();

	const unsigned int numElements_; // number of elements
	const FloatType pitch_; // distance between the elements' centers (m)
	const FloatType fs_; // sampling frequency (sample/s)
	const FloatType widthElem_; // width of each element (m)
	const FloatType heightElem_; // height of each element (m)
	const unsigned int numDivWidth_; // number of divisions along the width of each element
	const FloatType noiseAmplitude_;
	const FloatType fsUp_; // upsampled frequency (sample/s)
	FloatType acqTime_; // acquisition time (s)
	FloatType c_; // propagation speed (m/s)
	FloatType invC_;
	const FloatType maxHIdx_; // maximum impulse response index
	std::size_t ascanLength_;
	std::vector<FloatType> xeList_; // x-coordinates of the elements' centers
	std::vector<FloatType> dadtExc_;
	std::vector<XY<FloatType> > xySubElemList_;
	std::vector<FloatType> convHtHr_;
	std::vector<FloatType> ht_;
	std::vector<FloatType> hr_;
	std::vector<FloatType> p_;
	std::vector<FloatType> deltaY2_; // auxiliary
	std::vector<FloatType> txDelays_;
	std::vector<unsigned int> activeRxElements_;
	std::vector<unsigned int> activeTxElements_;
	std::vector<std::complex<FloatType> > dadtExcFilterFreqCoeff_;
	std::vector<FloatType> ascanBuffer_;
	std::vector<XZ<FloatType> > reflectorList_; // coordinates of the reflectors
	Matrix2<IndexValue<unsigned int, FloatType> > idxCoeffTable_;
	MinstdPseudorandomNumberGenerator prng_;
	int lowPassFIRLength_;
	std::vector<FloatType> lowPassFIRFilter_;
};

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::loadFilter()
{
	if (UPSAMPLING_FACTOR == 4) {
		// Calculated in Octave using:
		//-----------------------------------------------------------------------------
		// % Coefficients of the low-pass FIR filter:
		// wp = 0.5 * (1 / 4);
		// ws = 1 / 4;
		// %dp = 0.001;
		// %ds = 0.001;
		// n = 54; % order
		// b = remez(n, [0 wp ws 1], [1 1 0 0]);
		//-----------------------------------------------------------------------------
		// UPSAMPLING_FACTOR must be = 4.
		lowPassFIRLength_ = 55; /* must be odd */
		FloatType coeffList[] = {
			-4.13312132792898e-04,
			 3.84910656352768e-04,
			 8.95384486596103e-04,
			 1.42658409818017e-03,
			 1.57267578839296e-03,
			 9.56099017098954e-04,
			-5.59378457342753e-04,
			-2.67821756822127e-03,
			-4.62997598283731e-03,
			-5.35858923838621e-03,
			-3.93311746409238e-03,
			-5.97100599219861e-05,
			 5.52131936388271e-03,
			 1.09834954784043e-02,
			 1.38409960829662e-02,
			 1.18173151063219e-02,
			 3.90528342502080e-03,
			-8.76884400969965e-03,
			-2.26822124005643e-02,
			-3.24980236871478e-02,
			-3.24567720471757e-02,
			-1.82256580858917e-02,
			 1.13866341566510e-02,
			 5.34565424400340e-02,
			 1.01168250947271e-01,
			 1.45263694388271e-01,
			 1.76384224234024e-01,
			 1.87607302744230e-01,
			 1.76384224234024e-01,
			 1.45263694388271e-01,
			 1.01168250947271e-01,
			 5.34565424400340e-02,
			 1.13866341566510e-02,
			-1.82256580858917e-02,
			-3.24567720471757e-02,
			-3.24980236871478e-02,
			-2.26822124005643e-02,
			-8.76884400969965e-03,
			 3.90528342502080e-03,
			 1.18173151063219e-02,
			 1.38409960829662e-02,
			 1.09834954784043e-02,
			 5.52131936388271e-03,
			-5.97100599219861e-05,
			-3.93311746409238e-03,
			-5.35858923838621e-03,
			-4.62997598283731e-03,
			-2.67821756822127e-03,
			-5.59378457342753e-04,
			 9.56099017098954e-04,
			 1.57267578839296e-03,
			 1.42658409818017e-03,
			 8.95384486596103e-04,
			 3.84910656352768e-04,
			-4.13312132792898e-04
		};
		lowPassFIRFilter_.assign(coeffList, coeffList + lowPassFIRLength_);
	} else if (UPSAMPLING_FACTOR == 8) {
		// Calculated in Octave using:
		//-----------------------------------------------------------------------------
		// % Coefficients of the low-pass FIR filter:
		// wp = 0.5 * (1 / 8);
		// ws = 1 / 8;
		// %dp = 0.001;
		// %ds = 0.001;
		// n = 108; % order
		// b = remez(n, [0 wp ws 1], [1 1 0 0]);
		//-----------------------------------------------------------------------------
		// UPSAMPLING_FACTOR must be = 8.
		lowPassFIRLength_ = 109; /* must be odd */
		FloatType coeffList[] = {
			-4.58133871546099e-04,
			 1.14241187981019e-04,
			 1.89870412476132e-04,
			 3.06423396932843e-04,
			 4.47240614047113e-04,
			 5.91376008764852e-04,
			 7.14049948441715e-04,
			 7.88338564201616e-04,
			 7.87893765663302e-04,
			 6.89955275199809e-04,
			 4.79905363111702e-04,
			 1.53828506658959e-04,
			-2.77956142513765e-04,
			-7.89314095258098e-04,
			-1.33824911672097e-03,
			-1.86889947897627e-03,
			-2.31553451135030e-03,
			-2.60840302903155e-03,
			-2.68125283836319e-03,
			-2.47964215194002e-03,
			-1.96948072278556e-03,
			-1.14467542878866e-03,
			-3.28601866821800e-05,
			 1.30171271181286e-03,
			 2.75878649015765e-03,
			 4.20680834377148e-03,
			 5.49159326180611e-03,
			 6.44942341386254e-03,
			 6.92215396139025e-03,
			 6.77467900217173e-03,
			 5.91222550500521e-03,
			 4.29623358262505e-03,
			 1.95721014089734e-03,
			-9.97805060396267e-04,
			-4.38033525230739e-03,
			-7.93440313173101e-03,
			-1.13394793353801e-02,
			-1.42362149085319e-02,
			-1.62500608716019e-02,
			-1.70205260908561e-02,
			-1.62321003691986e-02,
			-1.36441783574696e-02,
			-9.11804487852647e-03,
			-2.63521947820468e-03,
			 5.68811201210161e-03,
			 1.56009006058155e-02,
			 2.67251439972692e-02,
			 3.85747701508721e-02,
			 5.05837312500924e-02,
			 6.21416201844054e-02,
			 7.26341528274994e-02,
			 8.14850674492227e-02,
			 8.81967516142448e-02,
			 9.23855130352458e-02,
			 9.38093338655455e-02,
			 9.23855130352458e-02,
			 8.81967516142448e-02,
			 8.14850674492227e-02,
			 7.26341528274994e-02,
			 6.21416201844054e-02,
			 5.05837312500924e-02,
			 3.85747701508721e-02,
			 2.67251439972692e-02,
			 1.56009006058155e-02,
			 5.68811201210161e-03,
			-2.63521947820468e-03,
			-9.11804487852647e-03,
			-1.36441783574696e-02,
			-1.62321003691986e-02,
			-1.70205260908561e-02,
			-1.62500608716019e-02,
			-1.42362149085319e-02,
			-1.13394793353801e-02,
			-7.93440313173101e-03,
			-4.38033525230739e-03,
			-9.97805060396267e-04,
			 1.95721014089734e-03,
			 4.29623358262505e-03,
			 5.91222550500521e-03,
			 6.77467900217173e-03,
			 6.92215396139025e-03,
			 6.44942341386254e-03,
			 5.49159326180611e-03,
			 4.20680834377148e-03,
			 2.75878649015765e-03,
			 1.30171271181286e-03,
			-3.28601866821800e-05,
			-1.14467542878866e-03,
			-1.96948072278556e-03,
			-2.47964215194002e-03,
			-2.68125283836319e-03,
			-2.60840302903155e-03,
			-2.31553451135030e-03,
			-1.86889947897627e-03,
			-1.33824911672097e-03,
			-7.89314095258098e-04,
			-2.77956142513765e-04,
			 1.53828506658959e-04,
			 4.79905363111702e-04,
			 6.89955275199809e-04,
			 7.87893765663302e-04,
			 7.88338564201616e-04,
			 7.14049948441715e-04,
			 5.91376008764852e-04,
			 4.47240614047113e-04,
			 3.06423396932843e-04,
			 1.89870412476132e-04,
			 1.14241187981019e-04,
			-4.58133871546099e-04
		};
		lowPassFIRFilter_.assign(coeffList, coeffList + lowPassFIRLength_);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid upsampling factor: " << UPSAMPLING_FACTOR << '.');
	}
}

template<typename FloatType>
SimulatedAcquisitionDevice<FloatType>::SimulatedAcquisitionDevice(
	unsigned int numElements,
	FloatType pitch,
	FloatType samplingFrequency,
	FloatType widthElem,
	FloatType heightElem,
	unsigned int numDivWidth,
	FloatType noiseAmplitude)
		: numElements_(numElements)
		, pitch_(pitch)
		, fs_(samplingFrequency)
		, widthElem_(widthElem)
		, heightElem_(heightElem)
		, numDivWidth_(numDivWidth)
		, noiseAmplitude_(noiseAmplitude)
		, fsUp_(fs_ * UPSAMPLING_FACTOR)
		, acqTime_(0.0)
		, c_(0.0)
		, invC_(0.0)
		, maxHIdx_(static_cast<FloatType>(std::numeric_limits<unsigned int>::max() / sizeof(FloatType)))
		, ascanLength_(0)
		, xeList_(numElements_)
		, txDelays_(numElements_)
		, prng_(SIMULATED_ACQUISITION_DEVICE_PSEUDORANDOM_NUMBER_GENERATOR_SEED)
{
	if (numDivWidth_ == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of divisions in the element width is equal to zero.");
	}

	LOG_DEBUG << "[SimulatedAcquisitionDevice] maxHIdx_=" << maxHIdx_ << " defNumThreads=" << tbb::task_scheduler_init::default_num_threads();

	// Calculates the x-coordinates of the elements' centers.
	for (unsigned int i = 0; i < numElements_; ++i) {
		xeList_[i] = i * pitch_;
	}

	// Calculates the sub-elements' coordinates. The reference is the center of the element.
	const FloatType sizeDivWidth = widthElem_ / numDivWidth_;
	const unsigned int numDivHeight = static_cast<unsigned int>(std::ceil(heightElem_ / sizeDivWidth / 2.0)) * 2U; // must be even
	const FloatType sizeDivHeight = heightElem_ / numDivHeight;
	const FloatType halfW = sizeDivWidth * (numDivWidth_ - 1) / 2.0;
	std::vector<FloatType> xSubElemList;
	Util::fillSequenceWithSize(xSubElemList, -halfW, halfW, numDivWidth_);
	const FloatType halfH = sizeDivHeight * (numDivHeight - 1) / 2.0;
	std::vector<FloatType> ySubElemList;
	Util::fillSequenceWithSize<FloatType>(ySubElemList, sizeDivHeight / 2.0, halfH, numDivHeight / 2); // uses only half of the element due to the simetry (ye must be = 0)

	xySubElemList_.resize(xSubElemList.size() * ySubElemList.size());
	for (unsigned int iw = 0, sizeW = xSubElemList.size(); iw < sizeW; ++iw) {
		for (unsigned int ih = 0, sizeH = ySubElemList.size(); ih < sizeH; ++ih) {
			XY<FloatType> center;
			center.x = xSubElemList[iw];
			center.y = ySubElemList[ih];
			xySubElemList_[iw + ih * numDivWidth_] = center;
		}
	}

	deltaY2_.resize(xySubElemList_.size());
	// For each sub-element:
	for (unsigned int ise = 0, size = xySubElemList_.size(); ise < size; ++ise) {
		const FloatType deltaY = /*yr*/ - xySubElemList_[ise].y;
		deltaY2_[ise] = deltaY * deltaY;
	}

	if (USE_LOWPASS_FIR_FILTER == 1) {
		loadFilter();
	}
}

template<typename FloatType>
SimulatedAcquisitionDevice<FloatType>::~SimulatedAcquisitionDevice()
{
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setAcquisitionTime(FloatType acqTime)
{
	FloatType len = std::ceil(acqTime * fs_);
	if (len > maxHIdx_) {
		THROW_EXCEPTION(InvalidValueException, "The acquisition time is too long.");
	}

	ascanLength_ = static_cast<std::size_t>(len);
}

// Converts the excitation velocity to da/dt.
// If the low-pass filter is used, the offset caused by it will be returned.
template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::prepareExcitationDadt(const std::vector<FloatType>& vExc, std::vector<FloatType>& dadtExc)
{
	const FloatType period = 1.0 / fsUp_;

	// Calculates the excitation acceleration.
	std::vector<FloatType> aExc;
	Util::centralDiff(vExc, period, aExc);

	if (USE_LOWPASS_FIR_FILTER == 1) {
		// Calculates the excitation d(acceleration)/dt.
		std::vector<FloatType> unfilteredDadtExc;
		Util::centralDiff(aExc, period, unfilteredDadtExc);
		std::vector<std::complex<FloatType> > filterFreqCoeff;

		FFTWFilter2<FloatType> f;
		f.setCoefficients(lowPassFIRFilter_, filterFreqCoeff);
		f.filter(filterFreqCoeff, unfilteredDadtExc, dadtExc);
	} else {
		Util::centralDiff(aExc, period, dadtExc);
	}

#if 0
	// Creates the excitation velocity waveform.
	const FloatType end = 1.5 / centerFrequency;
	const FloatType period = 1.0 / fsUp_;
	const FloatType twoPi = 2.0 * PI;
	const unsigned int numVExcPoints = static_cast<unsigned int>(end / period);
	dadtExc.resize(numVExcPoints);
	for (unsigned int i = 0; i < numVExcPoints; ++i) {
		const FloatType t = period * i;
		dadtExc[i] = std::sin(twoPi * centerFrequency * t);
	}
#endif
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setExcitationWaveform(FloatType centerFrequency)
{
	std::vector<FloatType> vExc;
	//createExcitationVelocity(centerFrequency, vExc);
	Waveform::createPulseA(centerFrequency, fsUp_, vExc);

	prepareExcitationDadt(vExc, dadtExc_);
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setExcitationWaveform(const std::vector<FloatType>& dadtExc)
{
	dadtExc_ = dadtExc;
}

template<typename FloatType>
const std::vector<FloatType>&
SimulatedAcquisitionDevice<FloatType>::getAscan()
{
	if (ascanBuffer_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The active receive elements have not been defined.");
	}
	if (dadtExc_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The excitation waveform has not been set.");
	}
	if (c_ == 0.0) {
		THROW_EXCEPTION(InvalidStateException, "The propagation speed has not been set.");
	}
	if (reflectorList_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "The list of reflectors is empty.");
	}

	std::fill(ascanBuffer_.begin(), ascanBuffer_.end(), FloatType(0));
	std::vector<std::complex<FloatType> > convDadtExcHtFilterFreqCoeff;
	FFTWFilter2<FloatType> filter1;
	filter1.setCoefficients(dadtExc_, dadtExcFilterFreqCoeff_);

	// For each reflector:
	for (std::size_t idxRef = 0, size = reflectorList_.size(); idxRef < size; ++idxRef) {
		const XZ<FloatType>& reflector = reflectorList_[idxRef];
		const FloatType xr = reflector.x;
		const FloatType zr = reflector.z;
		const FloatType zr2 = zr * zr;

		unsigned int idxMinTx = std::numeric_limits<unsigned int>::max();
		unsigned int idxMaxTx = 0;

		// For each active transmit element:
		for (std::size_t idxActiveTx = 0, sizeActiveTx = activeTxElements_.size(); idxActiveTx < sizeActiveTx; ++idxActiveTx) {
			const unsigned int activeTxElem = activeTxElements_[idxActiveTx];
			const FloatType deltaXCenter = xr - xeList_[activeTxElem];

			// For each sub-element:
			for (std::size_t idxSubElem = 0, size = xySubElemList_.size(); idxSubElem < size; ++idxSubElem) {
				const FloatType deltaX = deltaXCenter - xySubElemList_[idxSubElem].x;
				const FloatType r = std::sqrt(deltaX * deltaX + deltaY2_[idxSubElem] + zr2);
				const FloatType delaySamples = (r * invC_ + txDelays_[activeTxElem]) * fsUp_;
				if (delaySamples > maxHIdx_) {
					THROW_EXCEPTION(InvalidValueException, "The delay is too long.");
				}
				const unsigned int idx = static_cast<unsigned int>(0.5 + delaySamples);
				const FloatType coeff = 1.0 / r; //TODO: window? amplitudeList = wSubElem(:) ./ r;
				const IndexValue<unsigned int, FloatType> iv = { idx, coeff };
				idxCoeffTable_(idxActiveTx, idxSubElem) = iv;
				if (idx < idxMinTx) idxMinTx = idx;
				if (idx > idxMaxTx) idxMaxTx = idx;
			}
		}

		// Calculates the impulse response in transmission.
		ht_.assign(idxMaxTx - idxMinTx + 1U, 0.0);
		// For each active transmit element:
		for (std::size_t idxActiveTx = 0, sizeActiveTx = activeTxElements_.size(); idxActiveTx < sizeActiveTx; ++idxActiveTx) {
			// Combines the impulse response for each sub-element.
			// Stores only the useful part of h (from idxMin to idxMax).
			typename Matrix2<IndexValue<unsigned int, FloatType> >::ConstDim2Interval interval = idxCoeffTable_.dim2Interval(idxActiveTx);
			for (typename Matrix2<IndexValue<unsigned int, FloatType> >::ConstDim2Iterator iter = interval.first; iter != interval.second; ++iter) {
				ht_[iter->index - idxMinTx] += iter->value; // adds the coefficient
			}
		}

		std::vector<FloatType> convDadtExcHt;
		filter1.filter(dadtExcFilterFreqCoeff_, ht_, convDadtExcHt);

		FFTWFilter2<FloatType> filter2;
		filter2.setCoefficients(convDadtExcHt, convDadtExcHtFilterFreqCoeff);

		ThreadData data;
		data.filter = filter2;
		tbb::enumerable_thread_specific<ThreadData> dataTLS(data);

		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, activeRxElements_.size()),
				CalculateAscan(
						dataTLS,
						idxMinTx,
						invC_,
						fsUp_,
						maxHIdx_,
						ascanLength_,
						xr,
						zr,
						convDadtExcHt.size(),
						UPSAMPLING_FACTOR,
						lowPassFIRLength_,
						activeRxElements_,
						xySubElemList_,
						xeList_,
						deltaY2_,
						convDadtExcHtFilterFreqCoeff,
						ascanBuffer_));
	}

	if (noiseAmplitude_ != 0.0) {
		// Adds noise.
		const FloatType a = 2.0 * noiseAmplitude_;
		for (typename std::vector<FloatType>::iterator iter = ascanBuffer_.begin(), end = ascanBuffer_.end(); iter != end; ++iter) {
			*iter += (prng_.get() - 0.5) * a;
		}
	}

	return ascanBuffer_;
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setActiveRxElements(const std::vector<bool>& mask)
{
	if (ascanLength_ == 0) {
		THROW_EXCEPTION(InvalidStateException, "The acquisition time has not been set.");
	}
	//TODO: check size(mask) == numElements
	//TODO: check size/limits / acqTime_ * fs_

	activeRxElements_.clear();
	for (unsigned int i = 0, size = mask.size(); i < size; ++i) {
		if (mask[i]) {
			activeRxElements_.push_back(i);
		}
	}
	if (activeRxElements_.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "There is no active receive element.");
	}

	ascanBuffer_.resize(activeRxElements_.size() * ascanLength_);//TODO: remove to avoid copies?
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setActiveTxElements(const std::vector<bool>& mask)
{
	//TODO: check size(mask) == numElements
	//TODO: check size/limits

	activeTxElements_.clear();
	for (unsigned int i = 0, size = mask.size(); i < size; ++i) {
		if (mask[i]) {
			activeTxElements_.push_back(i);
		}
	}
	if (activeTxElements_.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "There is no active transmit element.");
	}

	idxCoeffTable_.resize(activeTxElements_.size(), xySubElemList_.size());
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setTxDelays(const std::vector<FloatType>& delays)
{
	//TODO: check sizes

	std::copy(delays.begin(), delays.end(), txDelays_.begin());
}

template<typename FloatType>
void
SimulatedAcquisitionDevice<FloatType>::setTxFocalPoint(FloatType xf, FloatType zf)
{
	if (c_ == 0.0) {
		THROW_EXCEPTION(InvalidStateException, "The propagation speed has not been set.");
	}

	FloatType maxDt = -1.0;

	// For each active transmit element:
	for (std::vector<unsigned int>::const_iterator iter = activeTxElements_.begin(); iter != activeTxElements_.end(); ++iter) {
		const FloatType dx = xf - xeList_[*iter];
		const FloatType r = std::sqrt(dx * dx + zf * zf);
		const FloatType dt = r * invC_; // transit time from the element to the focal point
		if (dt > maxDt) maxDt = dt;
		txDelays_[*iter] = dt;
	}

	// For each active transmit element:
	for (std::vector<unsigned int>::const_iterator iter = activeTxElements_.begin(); iter != activeTxElements_.end(); ++iter) {
		txDelays_[*iter] = maxDt - txDelays_[*iter];
	}
}



template<typename FloatType>
class SimulatedAcquisitionDevice<FloatType>::CalculateAscan {
public:
	CalculateAscan(
		tbb::enumerable_thread_specific<ThreadData>& dataTLS,
		unsigned int idxMinTx,
		FloatType invC,
		FloatType fsUp,
		FloatType maxHIdx,
		unsigned int ascanLength,
		FloatType xr,
		FloatType zr,
		std::size_t sizeConvDadtExcHt,
		int upsamplingFactor,
		int lowPassFIRLength,
		const std::vector<unsigned int>& activeRxElements,
		const std::vector<XY<FloatType> >& xySubElemList,
		const std::vector<FloatType>& xeList,
		const std::vector<FloatType>& deltaY2,
		const std::vector<std::complex<FloatType> >& convDadtExcHtFilterFreqCoeff,
		std::vector<FloatType>& ascanBuffer)
			: dataTLS_(dataTLS)
			, idxMinTx_(idxMinTx)
			, invC_(invC)
			, fsUp_(fsUp)
			, maxHIdx_(maxHIdx)
			, ascanLength_(ascanLength)
			, xr_(xr)
			, zr_(zr)
			, sizeConvDadtExcHt_(sizeConvDadtExcHt)
			, upsamplingFactor_(upsamplingFactor)
			, lowPassFIRLength_(lowPassFIRLength)
			, activeRxElements_(activeRxElements)
			, xySubElemList_(xySubElemList)
			, xeList_(xeList)
			, deltaY2_(deltaY2)
			, convDadtExcHtFilterFreqCoeff_(convDadtExcHtFilterFreqCoeff)
			, ascanBuffer_(ascanBuffer)
	{
	}

	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		typename tbb::enumerable_thread_specific<ThreadData>::reference data = dataTLS_.local();

		if (data.idxCoeffList.size() != xySubElemList_.size()) {
			data.idxCoeffList.resize(xySubElemList_.size());
		}

		const FloatType zr2 = zr_ * zr_;
		const FloatType invCFsUp = invC_ * fsUp_;

		// For each active receive element:
		for (std::size_t idxActiveRx = r.begin(); idxActiveRx != r.end(); ++idxActiveRx) {
			const unsigned int activeRxElem = activeRxElements_[idxActiveRx];
			const FloatType deltaXCenter = xr_ - xeList_[activeRxElem];
			unsigned int idxMinRx = std::numeric_limits<unsigned int>::max();
			unsigned int idxMaxRx = 0;

			// For each sub-element:
			for (std::size_t idxSubElem = 0, size = xySubElemList_.size(); idxSubElem < size; ++idxSubElem) {
				const FloatType deltaX = deltaXCenter - xySubElemList_[idxSubElem].x;
				const FloatType r = std::sqrt(deltaX * deltaX + deltaY2_[idxSubElem] + zr2);
				const FloatType delaySamples = r * invCFsUp;
				if (delaySamples > maxHIdx_) {
					THROW_EXCEPTION(InvalidValueException, "The delay is too long.");
				}
				const unsigned int idx = static_cast<unsigned int>(0.5 + delaySamples);
				const FloatType coeff = 1.0 / r; //TODO: window? amplitudeList = wSubElem(:) ./ r;
				const IndexValue<unsigned int, FloatType> iv = { idx, coeff };
				data.idxCoeffList[idxSubElem] = iv;
				if (idx < idxMinRx) idxMinRx = idx;
				if (idx > idxMaxRx) idxMaxRx = idx;
			}

			// Calculates the impulse response in reception.
			data.hr.assign(idxMaxRx - idxMinRx + 1U, 0.0);
			// Combines the impulse response for each sub-element.
			// Stores only the useful part of h (from idxMin to idxMax).
			for (typename std::vector<IndexValue<unsigned int, FloatType> >::const_iterator iter = data.idxCoeffList.begin(), iterEnd = data.idxCoeffList.end(); iter != iterEnd; ++iter) {
				data.hr[iter->index - idxMinRx] += iter->value; // adds the coefficient
			}

			data.filter.filter(convDadtExcHtFilterFreqCoeff_, data.hr, data.p);

			const unsigned int offsetAscanElem = idxActiveRx * ascanLength_;

			// pOffset: offset of p (sampling frequency: fsUp_).
			unsigned int offsetP;
			if (USE_LOWPASS_FIR_FILTER == 1) {
				const unsigned int offsetP1 = idxMinTx_ + idxMinRx;
				const unsigned int offsetP2 = (lowPassFIRLength_ - 1) / 2;
				if (offsetP2 > offsetP1) {
					THROW_EXCEPTION(InvalidValueException, "Negative offset of p.");
				}
				offsetP = offsetP1 - offsetP2; // the offset of the result of the convolutions
			} else {
				offsetP = idxMinTx_ + idxMinRx; // the offset of the result of the convolutions
			}
			const unsigned int modOffsetP = offsetP % upsamplingFactor_;
			// offsetPAdj: adjusted offset of p, must be a multiple of SIMULATED_ACQUISITION_DEVICE_UPSAMPLING_FACTOR (sampling frequency: fsUp_).
			const unsigned int offsetPAdj = (modOffsetP != 0) ? offsetP - modOffsetP + upsamplingFactor_ : offsetP;
			const unsigned int offsetAscan = offsetAscanElem + offsetPAdj / upsamplingFactor_;

			// Accumulates the pressure, applying the offset and downsampling.
			for (std::size_t i = offsetAscan, iEnd = offsetAscanElem + ascanLength_, j = offsetPAdj - offsetP, jEnd = data.p.size();
					(i < iEnd) && (j < jEnd);
					++i, j += upsamplingFactor_) {
				ascanBuffer_[i] += data.p[j];
			}
		}
	}
private:
	tbb::enumerable_thread_specific<ThreadData>& dataTLS_;
	const unsigned int idxMinTx_;
	const FloatType invC_;
	const FloatType fsUp_;
	const FloatType maxHIdx_;
	const unsigned int ascanLength_;
	const FloatType xr_;
	const FloatType zr_;
	const std::size_t sizeConvDadtExcHt_;
	const int upsamplingFactor_;
	const int lowPassFIRLength_;
	const std::vector<unsigned int>& activeRxElements_;
	const std::vector<XY<FloatType> >& xySubElemList_;
	const std::vector<FloatType>& xeList_;
	const std::vector<FloatType>& deltaY2_;
	const std::vector<std::complex<FloatType> >& convDadtExcHtFilterFreqCoeff_;
	std::vector<FloatType>& ascanBuffer_;
};

} // namespace Lab

#endif /* SIMULATEDACQUISITIONDEVICE_H_ */
