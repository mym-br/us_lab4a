/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

#ifndef PSEUDORANDOMNUMBERGENERATOR_H_
#define PSEUDORANDOMNUMBERGENERATOR_H_

#include <random>



namespace Lab {

class PseudorandomNumberGenerator {
public:
	PseudorandomNumberGenerator();
	~PseudorandomNumberGenerator();

	// Return one value of the pseudorandom sequence in [0.0,1.0).
	double get();
private:
	std::mt19937 engine_;
	std::uniform_real_distribution<double> dist_;
};

} // namespace Lab

#endif /* PSEUDORANDOMNUMBERGENERATOR_H_ */
