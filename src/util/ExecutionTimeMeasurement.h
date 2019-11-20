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

#ifndef EXECUTIONTIMEMEASUREMENT_H
#define EXECUTIONTIMEMEASUREMENT_H

//#define EXECUTION_TIME_MEASUREMENT_ACTIVE 1
#define EXECUTION_TIME_MEASUREMENT_ITERATIONS 10

#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
# include "Log.h"
# include "MeasurementList.h"
# include "Timer.h"
// Show information from a MeasurementList object.
# define EXECUTION_TIME_MEASUREMENT_LOG_TIMES(NAME, OBJ) \
	LOG_INFO << std::setprecision(15) << "TIME " NAME " mean=" << OBJ.arithmeticMean() << \
	" std=" << OBJ.standardDeviation() << \
	" min=" << OBJ.minimum() << \
	" max=" << OBJ.maximum() << \
	" std/mean=" << OBJ.standardDeviation() / OBJ.arithmeticMean()

// Show information from a MeasurementList object.
# define EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N(NAME, OBJ, N) \
	LOG_INFO << std::setprecision(15) << "TIME " NAME " mean=" << (OBJ.arithmeticMean() * (N)) << \
	" std=" << (OBJ.standardDeviation() * (N)) << \
	" min=" << (OBJ.minimum() * (N)) << \
	" max=" << (OBJ.maximum() * (N)) << \
	" std/mean=" << OBJ.standardDeviation() / OBJ.arithmeticMean()
#endif

#endif // EXECUTIONTIMEMEASUREMENT_H
