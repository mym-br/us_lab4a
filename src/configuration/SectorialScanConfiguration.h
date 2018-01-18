#ifndef SECTORIALSCANCONFIGURATION_H
#define SECTORIALSCANCONFIGURATION_H



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
};

} // namespace Lab

#endif // SECTORIALSCANCONFIGURATION_H
