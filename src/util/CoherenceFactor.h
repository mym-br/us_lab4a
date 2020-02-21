/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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
#include <memory>
#include <string>

#include "PseudorandomNumberGenerator.h"
#include "ParameterMap.h"
#include "Statistics.h"
#include "Util.h"

// SignCoherenceFactor
// PhaseCoherenceFactor
//
// Implementations of the coherence factors in:
// Camacho, J.
// Parrilla, M.
// Fritsch, C.
// Phase Coherence Imaging.
// IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control,
// vol. 56, no. 5, pp. 958-974, 2009.
// DOI: 10.1109/TUFFC.2009.1128
//
// PRNGPhaseCoherenceFactor
//
// Similar to PhaseCoherenceFactor, but sets a pseudorandom phase when
// the signal is zero.
// Used when the signal has no noise, like in simulated data.

namespace Lab {

template<typename TFloat> class SignCoherenceFactor;
template<typename TFloat> class PhaseCoherenceFactor;
template<typename TFloat> class PRNGPhaseCoherenceFactor;

template<typename TFloat>
class CoherenceFactor {
public:
	CoherenceFactor() = default;
	virtual ~CoherenceFactor() = default;

	virtual std::unique_ptr<CoherenceFactor<TFloat>> clone() const = 0;
	virtual TFloat calculate(const TFloat* data, unsigned int size) = 0;

	static std::unique_ptr<CoherenceFactor<TFloat>> get(const ParameterMap& pm);
private:
	CoherenceFactor(const CoherenceFactor&) = delete;
	CoherenceFactor& operator=(const CoherenceFactor&) = delete;
	CoherenceFactor(CoherenceFactor&&) = delete;
	CoherenceFactor& operator=(CoherenceFactor&&) = delete;
};

template<typename TFloat>
std::unique_ptr<CoherenceFactor<TFloat>>
CoherenceFactor<TFloat>::get(const ParameterMap& pm)
{
	const auto coherenceFactorMethod = pm.value<std::string>("coherence_factor_method");
	if (coherenceFactorMethod == "none") {
		return nullptr;
	} else if (coherenceFactorMethod == "sign_coherence_factor") {
		const auto p = pm.value<TFloat>("sign_coherence_factor_p", 0.0, 100.0);
		return std::make_unique<SignCoherenceFactor<TFloat>>(p);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid coherence factor method: " << coherenceFactorMethod << '.');
	}
}

//=============================================================================

// Copy constructible and assignable.
// May be allocated on the stack.
template<typename TFloat>
class CoherenceFactorProcessor {
public:
	CoherenceFactorProcessor() = default;
	CoherenceFactorProcessor(const ParameterMap& pm)
			: cf_(CoherenceFactor<TFloat>::get(pm)) {
	}
	CoherenceFactorProcessor(const CoherenceFactorProcessor& o)
			: cf_(o.cf_ ? o.cf_->clone() : nullptr) {
	}
	CoherenceFactorProcessor& operator=(const CoherenceFactorProcessor& o) {
		if (&o != this) {
			if (o.cf_) {
				cf_ = o.cf_->clone();
			} else {
				cf_.reset();
			}
		}
		return *this;
	}
	~CoherenceFactorProcessor() = default;

	TFloat calculate(const TFloat* data, unsigned int size) {
		if (cf_) {
			return cf_->calculate(data, size);
		} else {
			THROW_EXCEPTION(InvalidStateException, "The object has not been initialized.");
		}
	}
	bool enabled() const { return cf_.get() != nullptr; }
private:
	CoherenceFactorProcessor(CoherenceFactorProcessor&&) = delete;
	CoherenceFactorProcessor& operator=(CoherenceFactorProcessor&&) = delete;

	std::unique_ptr<CoherenceFactor<TFloat>> cf_;
};

//=============================================================================

template<typename TFloat>
class AnalyticSignalCoherenceFactor {
public:
	AnalyticSignalCoherenceFactor() = default;
	virtual ~AnalyticSignalCoherenceFactor() = default;

	virtual std::unique_ptr<AnalyticSignalCoherenceFactor<TFloat>> clone() const = 0;
	virtual TFloat calculate(const std::complex<TFloat>* data, unsigned int size) = 0;
	virtual void getConstants(std::vector<TFloat>& list) const = 0;

	static std::unique_ptr<AnalyticSignalCoherenceFactor<TFloat>> get(const ParameterMap& pm);
private:
	AnalyticSignalCoherenceFactor(const AnalyticSignalCoherenceFactor&) = delete;
	AnalyticSignalCoherenceFactor& operator=(const AnalyticSignalCoherenceFactor&) = delete;
	AnalyticSignalCoherenceFactor(AnalyticSignalCoherenceFactor&&) = delete;
	AnalyticSignalCoherenceFactor& operator=(AnalyticSignalCoherenceFactor&&) = delete;
};

template<typename TFloat>
std::unique_ptr<AnalyticSignalCoherenceFactor<TFloat>>
AnalyticSignalCoherenceFactor<TFloat>::get(const ParameterMap& pm)
{
	const auto coherenceFactorMethod = pm.value<std::string>("analytic_signal_coherence_factor_method");
	if (coherenceFactorMethod == "none") {
		return nullptr;
	} else if (coherenceFactorMethod == "phase_coherence_factor") {
		const auto gamma = pm.value<TFloat>("phase_coherence_factor_gamma", 0.0, 100.0);
		return std::make_unique<PhaseCoherenceFactor<TFloat>>(gamma);
	} else if (coherenceFactorMethod == "prng_phase_coherence_factor") {
		const auto gamma = pm.value<TFloat>("phase_coherence_factor_gamma", 0.0, 100.0);
		return std::make_unique<PRNGPhaseCoherenceFactor<TFloat>>(gamma);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid coherence factor method: " << coherenceFactorMethod << '.');
	}
}

//=============================================================================

// Copy constructible and assignable.
// May be allocated on the stack.
template<typename TFloat>
class AnalyticSignalCoherenceFactorProcessor {
public:
	AnalyticSignalCoherenceFactorProcessor() = default;
	AnalyticSignalCoherenceFactorProcessor(const ParameterMap& pm)
			: cf_(AnalyticSignalCoherenceFactor<TFloat>::get(pm)) {
	}
	AnalyticSignalCoherenceFactorProcessor(const AnalyticSignalCoherenceFactorProcessor& o)
			: cf_(o.cf_ ? o.cf_->clone() : nullptr) {
	}
	AnalyticSignalCoherenceFactorProcessor& operator=(const AnalyticSignalCoherenceFactorProcessor& o) {
		if (&o != this) {
			if (o.cf_) {
				cf_ = o.cf_->clone();
			} else {
				cf_.reset();
			}
		}
		return *this;
	}
	~AnalyticSignalCoherenceFactorProcessor() = default;

	TFloat calculate(const std::complex<TFloat>* data, unsigned int size) {
		if (cf_) {
			return cf_->calculate(data, size);
		} else {
			THROW_EXCEPTION(InvalidStateException, "The object has not been initialized.");
		}
	}
	bool enabled() const { return cf_.get() != nullptr; }
	const AnalyticSignalCoherenceFactor<TFloat>& implementation() const { return *cf_; }
private:
	AnalyticSignalCoherenceFactorProcessor(AnalyticSignalCoherenceFactorProcessor&&) = delete;
	AnalyticSignalCoherenceFactorProcessor& operator=(AnalyticSignalCoherenceFactorProcessor&&) = delete;

	std::unique_ptr<AnalyticSignalCoherenceFactor<TFloat>> cf_;
};

//=============================================================================

template<typename TFloat>
class SignCoherenceFactor : public CoherenceFactor<TFloat> {
public:
	SignCoherenceFactor(TFloat p) : p_(p) { }
	virtual ~SignCoherenceFactor() = default;

	virtual std::unique_ptr<CoherenceFactor<TFloat>> clone() const {
		return std::make_unique<SignCoherenceFactor>(p_);
	}
	virtual TFloat calculate(const TFloat* data, unsigned int size);
private:
	SignCoherenceFactor(const SignCoherenceFactor&) = delete;
	SignCoherenceFactor& operator=(const SignCoherenceFactor&) = delete;
	SignCoherenceFactor(SignCoherenceFactor&&) = delete;
	SignCoherenceFactor& operator=(SignCoherenceFactor&&) = delete;

	const TFloat p_;
};

template<typename TFloat>
TFloat
SignCoherenceFactor<TFloat>::calculate(const TFloat* data, unsigned int size)
{
	unsigned int n = 0;
	int signSum = 0;
	const TFloat* end = data + size;
	while (data != end) {
		if (*data++ >= 0) {
			++signSum;
		} else {
			--signSum;
		}
		++n;
	}

	const TFloat k = static_cast<TFloat>(signSum) / n;
	const TFloat m = 1 - std::sqrt(1 - k * k);
	return std::pow(std::abs(m), p_);
}

//=============================================================================

template<typename TFloat>
class PhaseCoherenceFactor : public AnalyticSignalCoherenceFactor<TFloat> {
public:
	PhaseCoherenceFactor(TFloat gamma)
		: gamma_(gamma)
		, sigma0_(pi / std::sqrt(3.0))
		, factor_(gamma / sigma0_) { }
	virtual ~PhaseCoherenceFactor() = default;

	virtual std::unique_ptr<AnalyticSignalCoherenceFactor<TFloat>> clone() const {
		return std::make_unique<PhaseCoherenceFactor>(gamma_);
	}
	virtual TFloat calculate(const std::complex<TFloat>* data, unsigned int size);
	virtual void getConstants(std::vector<TFloat>& list) const;
private:
	PhaseCoherenceFactor(const PhaseCoherenceFactor&) = delete;
	PhaseCoherenceFactor& operator=(const PhaseCoherenceFactor&) = delete;
	PhaseCoherenceFactor(PhaseCoherenceFactor&&) = delete;
	PhaseCoherenceFactor& operator=(PhaseCoherenceFactor&&) = delete;

	const TFloat gamma_;
	const TFloat sigma0_;
	const TFloat factor_;
	std::vector<TFloat> phi_;
	std::vector<TFloat> phiAux_;
};

template<typename TFloat>
TFloat
PhaseCoherenceFactor<TFloat>::calculate(const std::complex<TFloat>* data, unsigned int size)
{
	if (phi_.size() != size) {
		phi_.resize(size);
		phiAux_.resize(size);
	}

	for (unsigned int i = 0; i < size; ++i) {
		const std::complex<TFloat> c = data[i];
		phi_[i] = std::atan2(c.imag(), c.real());
		phiAux_[i] = phi_[i] + ((phi_[i] < 0) ? TFloat(pi) : -TFloat(pi));
	}

	const TFloat sf = std::min(
				Statistics::standardDeviation(&phi_[0]   , size),
				Statistics::standardDeviation(&phiAux_[0], size));
	return std::max<TFloat>(0, 1 - factor_ * sf);
}

template<typename TFloat>
void
PhaseCoherenceFactor<TFloat>::getConstants(std::vector<TFloat>& list) const
{
	list.clear();
	list.push_back(gamma_);
	list.push_back(sigma0_);
	list.push_back(factor_);
}

//=============================================================================

template<typename TFloat>
class PRNGPhaseCoherenceFactor : public AnalyticSignalCoherenceFactor<TFloat> {
public:
	PRNGPhaseCoherenceFactor(TFloat gamma)
		: gamma_(gamma)
		, sigma0_(pi / std::sqrt(3.0))
		, factor_(gamma / sigma0_) { }
	virtual ~PRNGPhaseCoherenceFactor() = default;

	virtual std::unique_ptr<AnalyticSignalCoherenceFactor<TFloat>> clone() const {
		return std::make_unique<PRNGPhaseCoherenceFactor>(gamma_);
	}
	virtual TFloat calculate(const std::complex<TFloat>* data, unsigned int size);
	virtual void getConstants(std::vector<TFloat>& list) const;
private:
	PRNGPhaseCoherenceFactor(const PRNGPhaseCoherenceFactor&) = delete;
	PRNGPhaseCoherenceFactor& operator=(const PRNGPhaseCoherenceFactor&) = delete;
	PRNGPhaseCoherenceFactor(PRNGPhaseCoherenceFactor&&) = delete;
	PRNGPhaseCoherenceFactor& operator=(PRNGPhaseCoherenceFactor&&) = delete;

	PseudorandomNumberGenerator prng_;
	const TFloat gamma_;
	const TFloat sigma0_;
	const TFloat factor_;
	std::vector<TFloat> phi_;
	std::vector<TFloat> phiAux_;
};

template<typename TFloat>
TFloat
PRNGPhaseCoherenceFactor<TFloat>::calculate(const std::complex<TFloat>* data, unsigned int size)
{
	if (phi_.size() != size) {
		phi_.resize(size);
		phiAux_.resize(size);
	}

	for (unsigned int i = 0; i < size; ++i) {
		const std::complex<TFloat> c = data[i];
		if (c == TFloat(0)) {
			phi_[i] = (2.0 * prng_.get() - 1.0) * pi;
		} else {
			phi_[i] = std::atan2(c.imag(), c.real());
		}
		phiAux_[i] = phi_[i] + ((phi_[i] < 0) ? TFloat(pi) : -TFloat(pi));
	}

	const TFloat sf = std::min(
				Statistics::standardDeviation(&phi_[0]   , size),
				Statistics::standardDeviation(&phiAux_[0], size));
	return std::max<TFloat>(0, 1 - factor_ * sf);
}

template<typename TFloat>
void
PRNGPhaseCoherenceFactor<TFloat>::getConstants(std::vector<TFloat>& list) const
{
	list.clear();
	list.push_back(gamma_);
	list.push_back(sigma0_);
	list.push_back(factor_);
}

} // namespace Lab

#endif // COHERENCEFACTOR_H
