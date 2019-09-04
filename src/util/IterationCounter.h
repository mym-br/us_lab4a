/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

#ifndef ITERATIONCOUNTER_H
#define ITERATIONCOUNTER_H

#include <atomic>

#include "Timer.h"

namespace Lab {

struct IterationCounter {
	static std::atomic_uint count;
	static unsigned int total;
	static Timer timer;

	static void reset(unsigned int numberOfIterations) {
		count = 0;
		total = numberOfIterations;
		timer.reset();
	}
	static void add(unsigned int n) {
		count.fetch_add(n);
	}
};

} // namespace Lab

#endif // ITERATIONCOUNTER_H
