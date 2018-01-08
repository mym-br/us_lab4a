#ifndef EXECUTIONTIMEMEASUREMENT_H
#define EXECUTIONTIMEMEASUREMENT_H

//#define EXECUTION_TIME_MEASUREMENT_ACTIVE 1
#define EXECUTION_TIME_MEASUREMENT_ITERATIONS 10

#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
# include "Log.h"
# include "MeasurementList.h"
# include "Timer.h"
// Shows information from a MeasurementList object.
# define EXECUTION_TIME_MEASUREMENT_LOG_TIMES(NAME, OBJ) \
	LOG_INFO << std::setprecision(15) << "TIME " NAME " mean=" << OBJ.arithmeticMean() << \
	" std=" << OBJ.standardDeviation() << \
	" min=" << OBJ.minimum() << \
	" max=" << OBJ.maximum() << \
	" std/mean=" << OBJ.standardDeviation() / OBJ.arithmeticMean()

// Shows information from a MeasurementList object.
# define EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N(NAME, OBJ, N) \
	LOG_INFO << std::setprecision(15) << "TIME " NAME " mean=" << (OBJ.arithmeticMean() * (N)) << \
	" std=" << (OBJ.standardDeviation() * (N)) << \
	" min=" << (OBJ.minimum() * (N)) << \
	" max=" << (OBJ.maximum() * (N)) << \
	" std/mean=" << OBJ.standardDeviation() / OBJ.arithmeticMean()
#endif

#endif // EXECUTIONTIMEMEASUREMENT_H
