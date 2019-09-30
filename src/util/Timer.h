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

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>



namespace Lab {

class Timer {
public:
	Timer() {
		reset();
	}
	void reset() {
		start_ = std::chrono::steady_clock::now();
	}
	double getTime() {
		std::chrono::duration<double> d = std::chrono::steady_clock::now() - start_;
		return d.count();
	}
private:
	std::chrono::time_point<std::chrono::steady_clock> start_;
};

} // namespace Lab

#endif /*TIMER_H_*/
