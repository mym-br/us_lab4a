/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#include <iomanip> /* setprecision */



#define EXECUTION_TIME_MEASUREMENT_ITERATIONS 10U

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
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



#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
# define BEGIN_EXECUTION_TIME_MEASUREMENT \
	{ MeasurementList<double> execTimeMeasML; \
	for (unsigned int nExecTimeMeas = 0; nExecTimeMeas < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++nExecTimeMeas) { \
	if (nExecTimeMeas <= 1U) { /* nExecTimeMeas = 0: initial reset, nExecTimeMeas = 1: ignores the first iteration */ \
		execTimeMeasML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS); \
	} \
	Timer execTimeMeasTimer;
# define END_EXECUTION_TIME_MEASUREMENT \
	execTimeMeasML.put(execTimeMeasTimer.getTime()); } \
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("total:", execTimeMeasML); }

# define BEGIN_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL(OBJ) \
	{ MeasurementList<double> execTimeMeasML; \
	for (unsigned int nExecTimeMeas = 0; nExecTimeMeas < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++nExecTimeMeas) { \
	if (nExecTimeMeas <= 1U) { /* nExecTimeMeas = 0: initial reset, nExecTimeMeas = 1: ignores the first iteration */ \
		execTimeMeasML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS); \
		(OBJ)->execTimeMeasReset(); \
	} \
	Timer execTimeMeasTimer;
# define END_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL(OBJ) \
	execTimeMeasML.put(execTimeMeasTimer.getTime()); } \
	(OBJ)->execTimeMeasShowResults(); \
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("total:", execTimeMeasML); }

# define BEGIN_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL_X_N(OBJ, N) \
	{ MeasurementList<double> execTimeMeasML; \
	for (unsigned int nExecTimeMeas = 0; nExecTimeMeas < EXECUTION_TIME_MEASUREMENT_ITERATIONS + 1; ++nExecTimeMeas) { \
	if (nExecTimeMeas <= 1U) { /* nExecTimeMeas = 0: initial reset, nExecTimeMeas = 1: ignores the first iteration */ \
		execTimeMeasML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS); \
		(OBJ)->execTimeMeasReset(N); \
	} \
	Timer execTimeMeasTimer;
# define END_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL_X_N(OBJ, N) \
	execTimeMeasML.put(execTimeMeasTimer.getTime()); } \
	(OBJ)->execTimeMeasShowResults(N); \
	EXECUTION_TIME_MEASUREMENT_LOG_TIMES("total:", execTimeMeasML); }
#else
# define BEGIN_EXECUTION_TIME_MEASUREMENT
# define END_EXECUTION_TIME_MEASUREMENT
# define BEGIN_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL(OBJ)
# define END_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL(OBJ)
# define BEGIN_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL_X_N(OBJ, N)
# define END_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL_X_N(OBJ, N)
#endif

#endif // EXECUTIONTIMEMEASUREMENT_H
