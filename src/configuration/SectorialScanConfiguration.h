#ifndef SECTORIALSCANCONFIGURATION_H
#define SECTORIALSCANCONFIGURATION_H

#include "ParameterMap.h"



namespace Lab {

template<typename FloatType>
struct SectorialScanConfiguration {
	unsigned int numElements;
	unsigned int numElementsMux;
	unsigned int baseElement;         // first element: 0
	unsigned int signalMode;          // 0: Raw / 1: Envelope
	FloatType pitch;                  // m
	FloatType centerFrequency;        // Hz
	FloatType gain;                   // dB
	FloatType propagationSpeed;       // m/s
	FloatType acquisitionDelay;       // s
	FloatType samplingFrequency;      // Hz
	FloatType focalEmissionDistance;  // m
	FloatType focalReceptionDistance; // m
	FloatType valueScale;

	// Sectorial scan grid.
	FloatType rangeStart; // m
	FloatType rangeEnd;   // m
	FloatType startAngle; // degree
	FloatType endAngle;   // degree
	FloatType angleStep;  // degree

	bool enableFocusing;

	void load(ConstParameterMapPtr pm);
};



template<typename FloatType>
void
SectorialScanConfiguration<FloatType>::load(ConstParameterMapPtr pm)
{
	numElementsMux         = pm->value<unsigned int>("num_elements_mux"        ,       8,      1024);
	numElements            = pm->value<unsigned int>("num_elements"            ,       8, numElementsMux);
	baseElement            = pm->value<unsigned int>("base_element"            ,       0, numElementsMux - numElements);
	signalMode             = pm->value<unsigned int>("signal_mode"             ,       0,         1);
	pitch                  = pm->value<FloatType>(   "pitch"                   , 0.01e-3,   10.0e-3);
	centerFrequency        = pm->value<FloatType>(   "center_frequency"        ,   100.0,   100.0e6);
	gain                   = pm->value<FloatType>(   "gain"                    ,     0.0,     100.0);
	propagationSpeed       = pm->value<FloatType>(   "propagation_speed"       ,   100.0,   10000.0);
	acquisitionDelay       = pm->value<FloatType>(   "acquisition_delay"       ,     0.0,       1.0);
	samplingFrequency      = pm->value<FloatType>(   "sampling_frequency"      ,   100.0,   200.0e6);
	focalEmissionDistance  = pm->value<FloatType>(   "focus_emission_distance" ,  1.0e-3, 1000.0e-3);
	focalReceptionDistance = pm->value<FloatType>(   "focus_reception_distance",  1.0e-3, 1000.0e-3);
	valueScale             = pm->value<FloatType>(   "value_scale"             ,     0.0,    1.0e10);

	// Sectorial scan grid.
	rangeStart             = pm->value<FloatType>(   "range_start"             ,  1.0e-3, 99.0);
	rangeEnd               = pm->value<FloatType>(   "range_end"               , rangeStart + 1.0e-2, 100.0);
	startAngle             = pm->value<FloatType>(   "start_angle"             ,   -90.0, 89.0);
	endAngle               = pm->value<FloatType>(   "end_angle"               , startAngle + 1.0, 90.0);
	angleStep              = pm->value<FloatType>(   "angle_step"              ,  1.0e-3, 10.0);

	enableFocusing         = pm->value<bool>(        "enable_focusing");

	if (numElementsMux & 0x01) {
		THROW_EXCEPTION(InvalidParameterException, "The value of num_elements_mux is not even.");
	}
	if (numElements & 0x01) {
		THROW_EXCEPTION(InvalidParameterException, "The value of num_elements is not even.");
	}
}

} // namespace Lab

#endif // SECTORIALSCANCONFIGURATION_H
