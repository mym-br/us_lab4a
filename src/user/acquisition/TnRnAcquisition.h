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

#ifndef TNRNACQUISITION_H
#define TNRNACQUISITION_H

#include "Matrix.h"

#include <vector>



namespace Lab {

template<typename FloatType>
class TnRnAcquisition {
public:
	typedef Matrix<FloatType> AcquisitionDataType;

	TnRnAcquisition() = default;
	virtual ~TnRnAcquisition() = default;

	virtual void prepare(unsigned int /*baseElement*/, const std::vector<FloatType>& /*txDelays*/) {}
	virtual void execute(AcquisitionDataType& acqData) = 0;
private:
	TnRnAcquisition(const TnRnAcquisition&) = delete;
	TnRnAcquisition& operator=(const TnRnAcquisition&) = delete;
	TnRnAcquisition(TnRnAcquisition&&) = delete;
	TnRnAcquisition& operator=(TnRnAcquisition&&) = delete;
};

} // namespace Lab

#endif // TNRNACQUISITION_H
