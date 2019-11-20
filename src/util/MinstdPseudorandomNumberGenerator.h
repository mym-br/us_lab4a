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

#ifndef MINSTDPSEUDORANDOMNUMBERGENERATOR_H_
#define MINSTDPSEUDORANDOMNUMBERGENERATOR_H_



namespace Lab {

class MinstdPseudorandomNumberGenerator {
public:
	MinstdPseudorandomNumberGenerator(long seed);
	~MinstdPseudorandomNumberGenerator();

	// Return one value of the pseudorandom sequence in ]0.0,1.0[.
	double get();
private:
	static constexpr double a = 16807.0;
	static constexpr double m = 2147483647.0;
	static constexpr double invM = 1.0 / m;

	double x_;
};

} // namespace Lab

#endif /* MINSTDPSEUDORANDOMNUMBERGENERATOR_H_ */
