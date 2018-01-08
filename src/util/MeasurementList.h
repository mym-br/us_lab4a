#ifndef MEASUREMENTLIST_H
#define MEASUREMENTLIST_H

#include <cstddef> /* std::size_t */
#include <vector>

#include "Exception.h"
#include "Statistics.h"



namespace Lab {

template<typename FloatType>
class MeasurementList {
public:
	MeasurementList() : position_(0), list_(0) { }

	void reset(std::size_t n) {
		position_ = 0;
		list_.assign(n, 0);
	}
	void put(FloatType measurement) {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ == list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is full.");
		}
		list_[position_++] = measurement;
	}

	const std::vector<FloatType>& list() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}
		return list_;
	}

	FloatType arithmeticMean() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		return Statistics::arithmeticMean(&list_[0], list_.size());
	}

	FloatType standardDeviation() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		return Statistics::standardDeviation(&list_[0], list_.size());
	}

	FloatType minimum() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		FloatType minValue = list_[0];
		for (std::size_t i = 1; i < list_.size(); ++i) {
			if (list_[i] < minValue) {
				minValue = list_[i];
			}
		}
		return minValue;
	}

	FloatType maximum() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		FloatType maxValue = list_[0];
		for (std::size_t i = 1; i < list_.size(); ++i) {
			if (list_[i] > maxValue) {
				maxValue = list_[i];
			}
		}
		return maxValue;
	}

private:
	std::size_t position_;
	std::vector<FloatType> list_;
};

} // namespace Lab

#endif // MEASUREMENTLIST_H
