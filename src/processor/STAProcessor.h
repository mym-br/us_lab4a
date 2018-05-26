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

#ifndef STAPROCESSOR_H_
#define STAPROCESSOR_H_

#include "Matrix2.h"
#include "XZValueFactor.h"



namespace Lab {

template<typename FloatType>
class STAProcessor {
public:
	STAProcessor() {}
	virtual ~STAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType>>& gridData) = 0;
};

} // namespace Lab

#endif /* STAPROCESSOR_H_ */
