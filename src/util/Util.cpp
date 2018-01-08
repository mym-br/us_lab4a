#include "Util.h"



namespace Lab {
namespace Util {

template<>
float
minValue()
{
	return -std::numeric_limits<float>::max();
}

template<>
double
minValue()
{
	return -std::numeric_limits<double>::max();
}

template<>
long double
minValue()
{
	return -std::numeric_limits<long double>::max();
}

} // namespace Util
} // namespace Lab
