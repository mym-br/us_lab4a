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

#define TIMER_POSIX_USE_CLOCK_GETTIME 1
#define TIMER_POSIX_USE_MONOTONIC_CLOCK 1

#ifdef _WIN32
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
# include <winbase.h>
# include "Exception.h"
#else
# ifdef TIMER_POSIX_USE_CLOCK_GETTIME
#  include <time.h>
#  ifdef TIMER_POSIX_USE_MONOTONIC_CLOCK
#   define TIMER_POSIX_CLOCK_GETTIME_CLOCK_ID CLOCK_MONOTONIC
#  else
#   define TIMER_POSIX_CLOCK_GETTIME_CLOCK_ID CLOCK_REALTIME
#  endif
# else
#  include <sys/time.h>
# endif
#endif



namespace Lab {

#ifdef _WIN32

struct Timer {
	LARGE_INTEGER ticksPerSecond;
	LARGE_INTEGER start;
	Timer() {
		if (!QueryPerformanceFrequency(&ticksPerSecond)) {
			throw Exception("Error in QueryPerformanceFrequency.");
		}
		if (!QueryPerformanceCounter(&start)) {
			throw Exception("Error in QueryPerformanceCounter.");
		}
	}
	double getTime() {
		LARGE_INTEGER cputime, end;
		QueryPerformanceCounter(&end);
		cputime.QuadPart = end.QuadPart - start.QuadPart;
		return static_cast<double>(cputime.QuadPart) / ticksPerSecond.QuadPart;
	}
};

#else /* _WIN32 */
# ifdef TIMER_POSIX_USE_CLOCK_GETTIME

class Timer {
public:
	Timer() : valid_() {
		reset();
	}
	void reset() {
		int rv = clock_gettime(TIMER_POSIX_CLOCK_GETTIME_CLOCK_ID, &start_);
		if (rv == 0) {
			valid_ = true;
		} else {
			valid_ = false;
		}
	}
	double getTime() {
		if (!valid_) return -1.0;

		struct timespec ts;
		int rv = clock_gettime(TIMER_POSIX_CLOCK_GETTIME_CLOCK_ID, &ts);
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

# else

class Timer {
public:
	Timer() : valid_() {
		reset();
	}
	void reset() {
		int rv = gettimeofday(&start_, nullptr); // its man page says it is obsolete
		if (rv == 0) {
			valid_ = true;
		} else {
			valid_ = false;
		}
	}
	double getTime() {
		if (!valid_) return -1.0;

		struct timeval tv;
		int rv = gettimeofday(&tv, nullptr); // its man page says it is obsolete
		if (rv == 0) {
			return static_cast<double>(tv.tv_sec - start_.tv_sec) + static_cast<double>(tv.tv_usec - start_.tv_usec) * 1.0e-6;
		} else {
			return -2.0;
		}
	}
private:
	bool valid_;
	struct timeval start_;
};

# endif
#endif /* _WIN32 */

} // namespace Lab

#endif /*TIMER_H_*/
