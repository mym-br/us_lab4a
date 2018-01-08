#ifndef WAVEFORM_H_
#define WAVEFORM_H_

#include <cmath>

#include "Util.h"



namespace Lab {
namespace Waveform {

template<typename FloatType> void createPulseA(FloatType centerFrequency /* Hz */, FloatType samplingFreq /* Hz */, std::vector<FloatType>& v);



template<typename FloatType>
void
createPulseA(FloatType centerFrequency, FloatType samplingFreq, std::vector<FloatType>& v)
{
	const FloatType end = 4.25 / centerFrequency;
	const FloatType period = 1.0 / samplingFreq;
	const FloatType twoPi = 2.0 * PI;
	const unsigned int numPoints = static_cast<unsigned int>(end / period) + 1;
	v.resize(numPoints);
	for (unsigned int i = 0; i < numPoints; ++i) {
		const FloatType t = period * i;
		v[i] = t * t * t * std::exp(-3.5 * centerFrequency * t) * std::cos(twoPi * centerFrequency * t);
	}
}

} // namespace Waveform
} // namespace Lab

#endif /* WAVEFORM_H_ */
