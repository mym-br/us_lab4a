/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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

#ifndef OCL_COHERENCE_FACTOR_H
#define OCL_COHERENCE_FACTOR_H

#include <string>

#include "OCLStatistics.h"



namespace Lab {
namespace OCLCoherenceFactor {

inline
std::string
code() {
	return OCLStatistics::code() + R"CLC(

MFloat
calcPCF(MFloat* re, MFloat* im, MFloat factor)
{
	MFloat phi[COHERENCE_FACTOR_N];
	MFloat phiAux[COHERENCE_FACTOR_N];

	for (unsigned int i = 0; i < COHERENCE_FACTOR_N; ++i) {
		phi[i] = atan2(im[i], re[i]);
	}
	for (unsigned int i = 0; i < COHERENCE_FACTOR_N; ++i) {
		phiAux[i] = phi[i] - copysign((MFloat) M_PI, phi[i]);
	}

	const MFloat sf = fmin(standardDeviation(phi), standardDeviation(phiAux));
	return fmax((MFloat) 0, (MFloat) 1 - factor * sf);
}

)CLC";
}

} // namespace OCLCoherenceFactor
} // namespace Lab

#endif // OCL_COHERENCE_FACTOR_H
