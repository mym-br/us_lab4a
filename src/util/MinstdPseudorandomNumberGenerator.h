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

	// Returns one value of the pseudorandom sequence ]0.0,1.0[.
	double get();
private:
	//MinstdPseudorandomNumberGenerator(const MinstdPseudorandomNumberGenerator&);
	//MinstdPseudorandomNumberGenerator& operator=(const MinstdPseudorandomNumberGenerator&);

	static const double a;
	static const double m;
	static const double invM;
	double x_;
};

} // namespace Lab

#endif /* MINSTDPSEUDORANDOMNUMBERGENERATOR_H_ */
