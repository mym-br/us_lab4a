#include "MinstdPseudorandomNumberGenerator.h"

#include <cmath>

#include "Exception.h"



namespace Lab {

const double MinstdPseudorandomNumberGenerator::a = 16807.0;
const double MinstdPseudorandomNumberGenerator::m = 2147483647.0;
const double MinstdPseudorandomNumberGenerator::invM = 1.0 / MinstdPseudorandomNumberGenerator::m;

MinstdPseudorandomNumberGenerator::MinstdPseudorandomNumberGenerator(long seed) : x_(seed)
{
	if (seed < 1 || seed >= m) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid seed: " << seed << " (1 <= seed < " << static_cast<long>(m) << ").");
	}
}

MinstdPseudorandomNumberGenerator::~MinstdPseudorandomNumberGenerator()
{
}

double
MinstdPseudorandomNumberGenerator::get()
{
	const double ax = a * x_;
	x_ = ax - m * std::floor(ax * invM);
	return x_ * invM;
}

} // namespace Lab
