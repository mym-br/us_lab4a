#ifndef LINEARINTERPOLATOR_H_
#define LINEARINTERPOLATOR_H_

#include <cstddef> /* std::size_t */



namespace Lab {

template<typename FloatType, unsigned int UpsamplingFactor>
class LinearInterpolator {
public:
	LinearInterpolator();
	~LinearInterpolator();

	void interpolate(const FloatType* input, std::size_t inputLength, FloatType* output);
private:
	LinearInterpolator(const LinearInterpolator&);
	LinearInterpolator& operator=(const LinearInterpolator&);
};

template<typename FloatType, unsigned int UpsamplingFactor>
LinearInterpolator<FloatType, UpsamplingFactor>::LinearInterpolator()
{
}

template<typename FloatType, unsigned int UpsamplingFactor>
LinearInterpolator<FloatType, UpsamplingFactor>::~LinearInterpolator()
{
}

// The argument "ouput" must point to an array of size inputLength * upsamplingFactor_.
template<typename FloatType, unsigned int UpsamplingFactor>
void
LinearInterpolator<FloatType, UpsamplingFactor>::interpolate(const FloatType* input, std::size_t inputLength, FloatType* output)
{
	for (std::size_t i = 0; i < inputLength - 1; ++i) {
		const FloatType step = (input[i + 1] - input[i]) / UpsamplingFactor;
		FloatType* destPtr = output + i * UpsamplingFactor;
		FloatType value = input[i];
		for (unsigned int j = 0; j < UpsamplingFactor; ++j) {
			*destPtr++ = value;
			value += step;
		}
	}

	// The last interval.
	const std::size_t i = inputLength - 1;
	const FloatType step = -input[i] / UpsamplingFactor;
	FloatType* destPtr = output + i * UpsamplingFactor;
	FloatType value = input[i];
	for (unsigned int j = 0; j < UpsamplingFactor; ++j) {
		*destPtr++ = value;
		value += step;
	}
}

} // namespace Lab

#endif /* LINEARINTERPOLATOR_H_ */
