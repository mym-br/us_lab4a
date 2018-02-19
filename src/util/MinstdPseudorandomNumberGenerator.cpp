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
