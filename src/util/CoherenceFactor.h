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

#ifndef COHERENCEFACTOR_H
#define COHERENCEFACTOR_H

#include <cmath> /* abs, atan2, pow, sqrt */
#include <string>

#include <boost/scoped_ptr.hpp>

#include "MinstdPseudorandomNumberGenerator.h"
#include "ParameterMap.h"
#include "Statistics.h"
#include "Util.h"



namespace Lab {

template<typename FloatType> class SignCoherenceFactor;
template<typename FloatType> class PhaseCoherenceFactor;
template<typename FloatType> class PRNGPhaseCoherenceFactor;

template<typename FloatType>
class CoherenceFactor {
public:
	CoherenceFactor() { }
	virtual ~CoherenceFactor() { }

	virtual CoherenceFactor<FloatType>* clone() const = 0;
	virtual FloatType calculate(const FloatType* data, unsigned int size) = 0;

	static CoherenceFactor<FloatType>* get(ConstParameterMapPtr pm);
private:
	CoherenceFactor(const CoherenceFactor&);
	CoherenceFactor& operator=(const CoherenceFactor&);
};

template<typename FloatType>
CoherenceFactor<FloatType>*
CoherenceFactor<FloatType>::get(ConstParameterMapPtr pm)
{
	if (!pm) {
		THROW_EXCEPTION(InvalidParameterException, "The parameter map has not been initialized.");
	}
	std::string coherenceFactorMethod = pm->value<std::string>("coherence_factor_method");
	if (coherenceFactorMethod == "none") {
		return 0;
	} else if (coherenceFactorMethod == "sign_coherence_factor") {
		const FloatType p = pm->value<FloatType>("sign_coherence_factor_p", 0.0, 100.0);
		return new SignCoherenceFactor<FloatType>(p);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid coherence factor method: " << coherenceFactorMethod << '.');
	}
}

//=============================================================================

// Copy constructible and assignable.
// May be allocated on the stack.
template<typename FloatType>
class CoherenceFactorProcessor {
public:
	CoherenceFactorProcessor() : cf_(0) { }
	CoherenceFactorProcessor(ConstParameterMapPtr pm)
			: cf_(CoherenceFactor<FloatType>::get(pm)) {
	}
	~CoherenceFactorProcessor() { }
	CoherenceFactorProcessor(const CoherenceFactorProcessor& o)
			: cf_(o.cf_ ? o.cf_->clone() : 0) {
	}
	CoherenceFactorProcessor& operator=(const CoherenceFactorProcessor& o) {
		if (&o != this) {
			cf_.reset(o.cf_ ? o.cf_->clone() : 0);
		}
		return *this;
	}
	FloatType calculate(const FloatType* data, unsigned int size) {
		if (cf_) {
			return cf_->calculate(data, size);
		} else {
			THROW_EXCEPTION(InvalidStateException, "The object has not been initialized.");
		}
	}
	bool enabled() { return cf_.get() != nullptr; }
private:
	boost::scoped_ptr<CoherenceFactor<FloatType> > cf_;
};

//=============================================================================

template<typename FloatType>
class AnalyticSignalCoherenceFactor {
public:
	AnalyticSignalCoherenceFactor() { }
	virtual ~AnalyticSignalCoherenceFactor() { }

	virtual AnalyticSignalCoherenceFactor<FloatType>* clone() const = 0;
	virtual FloatType calculate(const std::complex<FloatType>* data, unsigned int size) = 0;
	virtual void getConstants(std::vector<FloatType>& list) const = 0;

	static AnalyticSignalCoherenceFactor<FloatType>* get(ConstParameterMapPtr pm);
private:
	AnalyticSignalCoherenceFactor(const AnalyticSignalCoherenceFactor&);
	AnalyticSignalCoherenceFactor& operator=(const AnalyticSignalCoherenceFactor&);
};

template<typename FloatType>
AnalyticSignalCoherenceFactor<FloatType>*
AnalyticSignalCoherenceFactor<FloatType>::get(ConstParameterMapPtr pm)
{
	if (!pm) {
		THROW_EXCEPTION(InvalidParameterException, "The parameter map has not been initialized.");
	}
	std::string coherenceFactorMethod = pm->value<std::string>("analytic_signal_coherence_factor_method");
	if (coherenceFactorMethod == "none") {
		return 0;
	} else if (coherenceFactorMethod == "phase_coherence_factor") {
		const FloatType gamma = pm->value<FloatType>("phase_coherence_factor_gamma", 0.0, 100.0);
		return new PhaseCoherenceFactor<FloatType>(gamma);
	} else if (coherenceFactorMethod == "prng_phase_coherence_factor") {
		const FloatType gamma = pm->value<FloatType>("phase_coherence_factor_gamma", 0.0, 100.0);
		return new PRNGPhaseCoherenceFactor<FloatType>(gamma);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid coherence factor method: " << coherenceFactorMethod << '.');
	}
}

//=============================================================================

// Copy constructible and assignable.
// May be allocated on the stack.
template<typename FloatType>
class AnalyticSignalCoherenceFactorProcessor {
public:
	AnalyticSignalCoherenceFactorProcessor() : cf_(0) { }
	AnalyticSignalCoherenceFactorProcessor(ConstParameterMapPtr pm)
			: cf_(AnalyticSignalCoherenceFactor<FloatType>::get(pm)) {
	}
	~AnalyticSignalCoherenceFactorProcessor() { }
	AnalyticSignalCoherenceFactorProcessor(const AnalyticSignalCoherenceFactorProcessor& o)
			: cf_(o.cf_ ? o.cf_->clone() : 0) {
	}
	AnalyticSignalCoherenceFactorProcessor& operator=(const AnalyticSignalCoherenceFactorProcessor& o) {
		if (&o != this) {
			cf_.reset(o.cf_ ? o.cf_->clone() : 0);
		}
		return *this;
	}
	FloatType calculate(const std::complex<FloatType>* data, unsigned int size) {
		if (cf_) {
			return cf_->calculate(data, size);
		} else {
			THROW_EXCEPTION(InvalidStateException, "The object has not been initialized.");
		}
	}
	bool enabled() { return cf_.get() != nullptr; }
	const AnalyticSignalCoherenceFactor<FloatType>& implementation() const { return *cf_; }
private:
	boost::scoped_ptr<AnalyticSignalCoherenceFactor<FloatType>> cf_;
};

//=============================================================================

template<typename FloatType>
class SignCoherenceFactor : public CoherenceFactor<FloatType> {
public:
	SignCoherenceFactor(FloatType p) : p_(p) { }
	virtual ~SignCoherenceFactor() { }

	virtual CoherenceFactor<FloatType>* clone() const {
		return new SignCoherenceFactor(*this);
	}
	virtual FloatType calculate(const FloatType* data, unsigned int size);
private:
	SignCoherenceFactor(const SignCoherenceFactor& o)
		: CoherenceFactor<FloatType>()
		, p_(o.p_) { }
	SignCoherenceFactor& operator=(const SignCoherenceFactor& o) {
		if (&o != this) {
			p_ = o.p_;
		}
	}

	FloatType p_;
};

template<typename FloatType>
FloatType
SignCoherenceFactor<FloatType>::calculate(const FloatType* data, unsigned int size)
{
	unsigned int n = 0;
	int signSum = 0;
	const FloatType* end = data + size;
	while (data != end) {
		if (*data++ >= 0) {
			++signSum;
		} else {
			--signSum;
		}
		++n;
	}

	const FloatType k = static_cast<FloatType>(signSum) / n;
	const FloatType m = 1 - std::sqrt(1 - k * k);
	return std::pow(std::abs(m), p_);
}

//=============================================================================

template<typename FloatType>
class PhaseCoherenceFactor : public AnalyticSignalCoherenceFactor<FloatType> {
public:
	PhaseCoherenceFactor(FloatType gamma)
		: gamma_(gamma)
		, sigma0_(PI / std::sqrt(3.0))
		, factor_(gamma / sigma0_) { }
	virtual ~PhaseCoherenceFactor() { }

	virtual AnalyticSignalCoherenceFactor<FloatType>* clone() const {
		return new PhaseCoherenceFactor(*this);
	}
	virtual FloatType calculate(const std::complex<FloatType>* data, unsigned int size);
	virtual void getConstants(std::vector<FloatType>& paramList) const;
private:
	PhaseCoherenceFactor(const PhaseCoherenceFactor& o)
		: AnalyticSignalCoherenceFactor<FloatType>()
		, gamma_( o.gamma_)
		, sigma0_(o.sigma0_)
		, factor_(o.factor_)
		, phi_(   o.phi_)
		, phiAux_(o.phiAux_) { }
	PhaseCoherenceFactor& operator=(const PhaseCoherenceFactor& o) {
		if (&o != this) {
			gamma_  = o.gamma_;
			sigma0_ = o.sigma0_;
			factor_ = o.factor_;
			phi_    = o.phi_;
			phiAux_ = o.phiAux_;
		}
	}

	const FloatType gamma_;
	const FloatType sigma0_;
	const FloatType factor_;
	std::vector<FloatType> phi_;
	std::vector<FloatType> phiAux_;
};

template<typename FloatType>
FloatType
PhaseCoherenceFactor<FloatType>::calculate(const std::complex<FloatType>* data, unsigned int size)
{
	if (phi_.size() != size) {
		phi_.resize(size);
		phiAux_.resize(size);
	}

	for (unsigned int i = 0; i < size; ++i) {
		const std::complex<FloatType> c = data[i];
		phi_[i] = std::atan2(c.imag(), c.real());
	}
	for (unsigned int i = 0; i < size; ++i) {
		phiAux_[i] = phi_[i] + ((phi_[i] < 0) ? FloatType(PI) : -FloatType(PI));
	}

	const FloatType sf = std::min(
				Statistics::standardDeviation(&phi_[0]   , size),
				Statistics::standardDeviation(&phiAux_[0], size));
	return std::max<FloatType>(0, 1 - factor_ * sf);
}

template<typename FloatType>
void
PhaseCoherenceFactor<FloatType>::getConstants(std::vector<FloatType>& list) const
{
	list.clear();
	list.push_back(gamma_);
	list.push_back(sigma0_);
	list.push_back(factor_);
}

//=============================================================================

template<typename FloatType>
class PRNGPhaseCoherenceFactor : public AnalyticSignalCoherenceFactor<FloatType> {
public:
	PRNGPhaseCoherenceFactor(FloatType gamma)
		: prng_(1)
		, gamma_(gamma)
		, sigma0_(PI / std::sqrt(3.0))
		, factor_(gamma / sigma0_) { }
	virtual ~PRNGPhaseCoherenceFactor() { }

	virtual AnalyticSignalCoherenceFactor<FloatType>* clone() const {
		return new PRNGPhaseCoherenceFactor(*this);
	}
	virtual FloatType calculate(const std::complex<FloatType>* data, unsigned int size);
	virtual void getConstants(std::vector<FloatType>& paramList) const;
private:
	PRNGPhaseCoherenceFactor(const PRNGPhaseCoherenceFactor& o)
		: AnalyticSignalCoherenceFactor<FloatType>()
		, prng_(  o.prng_)
		, gamma_( o.gamma_)
		, sigma0_(o.sigma0_)
		, factor_(o.factor_)
		, phi_(   o.phi_)
		, phiAux_(o.phiAux_) { }
	PRNGPhaseCoherenceFactor& operator=(const PRNGPhaseCoherenceFactor& o) {
		if (&o != this) {
			prng_   = o.prng_;
			gamma_  = o.gamma_;
			sigma0_ = o.sigma0_;
			factor_ = o.factor_;
			phi_    = o.phi_;
			phiAux_ = o.phiAux_;
		}
	}

	MinstdPseudorandomNumberGenerator prng_;
	const FloatType gamma_;
	const FloatType sigma0_;
	const FloatType factor_;
	std::vector<FloatType> phi_;
	std::vector<FloatType> phiAux_;
};

template<typename FloatType>
FloatType
PRNGPhaseCoherenceFactor<FloatType>::calculate(const std::complex<FloatType>* data, unsigned int size)
{
	if (phi_.size() != size) {
		phi_.resize(size);
		phiAux_.resize(size);
	}

	for (unsigned int i = 0; i < size; ++i) {
		const std::complex<FloatType> c = data[i];
		if (c == FloatType(0)) {
			phi_[i] = (2.0 * prng_.get() - 1.0) * PI;
		} else {
			phi_[i] = std::atan2(c.imag(), c.real());
		}
		if (phi_[i] < 0) {
			phiAux_[i] = phi_[i] + FloatType(PI);
		} else {
			phiAux_[i] = phi_[i] - FloatType(PI);
		}
	}

	const FloatType sf = std::min(
				Statistics::standardDeviation(&phi_[0]   , size),
				Statistics::standardDeviation(&phiAux_[0], size));
	return std::max<FloatType>(0, 1 - factor_ * sf);
}

template<typename FloatType>
void
PRNGPhaseCoherenceFactor<FloatType>::getConstants(std::vector<FloatType>& list) const
{
	list.clear();
	list.push_back(gamma_);
	list.push_back(sigma0_);
	list.push_back(factor_);
}







//template<typename FloatType>
//class PCF2 {
//public:
//	PCF2(FloatType p) : p_(p) {}
//	~PCF2() {}

//	FloatType calculate(const std::complex<FloatType>* data, std::size_t size);
//private:
//	const FloatType p_;
//};

//template<typename FloatType>
//FloatType
//PCF2<FloatType>::calculate(const std::complex<FloatType>* data, std::size_t size)
//{
//	std::complex<FloatType> vectorialSum = 0.0;

//	for (std::size_t i = 0; i < size; ++i) {
//		std::complex<FloatType> c = data[i];
//		c /= std::abs(c); //FIXME division by zero
//		vectorialSum += c;
//	}

//	return std::pow(std::abs(vectorialSum) / size, p_);
//}



//template<typename FloatType>
//class PCF3 {
//public:
//	PCF3() {}
//	~PCF3() {}

//	FloatType calculate(const std::complex<FloatType>* data, std::size_t size);
//};

//template<typename FloatType>
//FloatType
//PCF3<FloatType>::calculate(const std::complex<FloatType>* data, std::size_t size)
//{
//	std::complex<FloatType> vectorialSum = 0.0;

//	for (std::size_t i = 0; i < size; ++i) {
//		std::complex<FloatType> c = data[i];
//		c /= std::abs(c);
//		vectorialSum += c;
//	}

//	const FloatType trA = 0.4;//TODO: use parameters
//	const FloatType trB = 0.7;

//	const FloatType f1 = std::abs(vectorialSum) / size;

//	if (f1 <= trA) {
//		return 0.0;
//	} else if (f1 >= trB) {
//		return 1.0;
//	} else {
//		const FloatType x = (f1 - trA) / (trB - trA);//TODO: check trB > trA
//		const FloatType x2 = x * x;
//		return x2 * (FloatType(-2) * x + FloatType(3));
//	}
//}


} // namespace Lab

#endif // COHERENCEFACTOR_H
