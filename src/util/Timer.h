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

#ifndef TIMER_H_
#define TIMER_H_

#include <time.h>



namespace Lab {

class Timer {
public:
	Timer() : valid_() {
		reset();
	}
	void reset() {
		int rv = clock_gettime(CLOCK_MONOTONIC, &start_);
		if (rv == 0) {
			valid_ = true;
		} else {
			valid_ = false;
		}
	}
	double getTime() {
		if (!valid_) return -1.0;

		struct timespec ts;
		int rv = clock_gettime(CLOCK_MONOTONIC, &ts);
		if (rv == 0) {
			return static_cast<double>(ts.tv_sec - start_.tv_sec) + static_cast<double>(ts.tv_nsec - start_.tv_nsec) * 1.0e-9;
		} else {
			return -2.0;
		}
	}
private:
	bool valid_;
	struct timespec start_;
};

} // namespace Lab

#endif /*TIMER_H_*/
