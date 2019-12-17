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

#ifndef MEASUREMENTLIST_H
#define MEASUREMENTLIST_H

#include <cstddef> /* std::size_t */
#include <vector>

#include "Exception.h"
#include "Statistics.h"



namespace Lab {

template<typename TFloat>
class MeasurementList {
public:
	MeasurementList() : position_(), list_() { }

	void reset(std::size_t n) {
		position_ = 0;
		list_.assign(n, 0);
	}
	void put(TFloat measurement) {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ == list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is full.");
		}
		list_[position_++] = measurement;
	}

	const std::vector<TFloat>& list() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}
		return list_;
	}

	TFloat arithmeticMean() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		return Statistics::arithmeticMean(&list_[0], list_.size());
	}

	TFloat standardDeviation() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		return Statistics::standardDeviation(&list_[0], list_.size());
	}

	TFloat minimum() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		TFloat minValue = list_[0];
		for (std::size_t i = 1; i < list_.size(); ++i) {
			if (list_[i] < minValue) {
				minValue = list_[i];
			}
		}
		return minValue;
	}

	TFloat maximum() const {
		if (list_.empty()) {
			THROW_EXCEPTION(InvalidStateException, "The list size is zero.");
		}
		if (position_ != list_.size()) {
			THROW_EXCEPTION(InvalidStateException, "The list is not full.");
		}

		TFloat maxValue = list_[0];
		for (std::size_t i = 1; i < list_.size(); ++i) {
			if (list_[i] > maxValue) {
				maxValue = list_[i];
			}
		}
		return maxValue;
	}

private:
	std::size_t position_;
	std::vector<TFloat> list_;
};

} // namespace Lab

#endif // MEASUREMENTLIST_H
