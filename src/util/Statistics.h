#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cmath>
#include <cstddef> /* std::size_t */



namespace Lab {
namespace Statistics {

template<typename FloatType> FloatType standardDeviation(const FloatType* data, std::size_t size);
template<typename FloatType> FloatType arithmeticMean(const FloatType* data, std::size_t size);



template<typename FloatType>
FloatType
standardDeviation(const FloatType* data, std::size_t size)
{
	FloatType sum = 0.0;
	const FloatType* end = data + size;
	const FloatType mean = arithmeticMean(data, size);
	while (data != end) {
		const FloatType e = *data++ - mean;
		sum += e * e;
	}
	return std::sqrt(sum / size);
}

template<typename FloatType>
FloatType
arithmeticMean(const FloatType* data, std::size_t size)
{
	FloatType sum = 0.0;
	const FloatType* end = data + size;
	while (data != end) {
		sum += *data++;
	}
	return sum / size;
}

} // namespace Statistics
} // namespace Lab

#endif /* STATISTICS_H_ */
