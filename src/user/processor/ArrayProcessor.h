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

#ifndef ARRAYPROCESSOR_H_
#define ARRAYPROCESSOR_H_

#include "Matrix.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename FloatType>
class ArrayProcessor {
public:
	ArrayProcessor() = default;
	virtual ~ArrayProcessor() = default;

	virtual void prepare(unsigned int /*baseElement*/) {}
	virtual void process(Matrix<XYZValueFactor<FloatType>>& gridData) = 0;
private:
	ArrayProcessor(const ArrayProcessor&) = delete;
	ArrayProcessor& operator=(const ArrayProcessor&) = delete;
	ArrayProcessor(ArrayProcessor&&) = delete;
	ArrayProcessor& operator=(ArrayProcessor&&) = delete;
};

} // namespace Lab

#endif /* ARRAYPROCESSOR_H_ */
