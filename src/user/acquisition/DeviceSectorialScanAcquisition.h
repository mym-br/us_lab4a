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

#ifndef DEVICESECTORIALSCANACQUISITION_H
#define DEVICESECTORIALSCANACQUISITION_H

#include "Matrix.h"
#include "XZValue.h"



namespace Lab {

// The image formation is executed on the acquisition device.
template<typename TFloat>
class DeviceSectorialScanAcquisition {
public:
	typedef Matrix<XZValue<TFloat>> AcquisitionDataType;

	DeviceSectorialScanAcquisition() = default;
	virtual ~DeviceSectorialScanAcquisition() = default;

	virtual void execute(AcquisitionDataType& acqData) = 0;
private:
	DeviceSectorialScanAcquisition(const DeviceSectorialScanAcquisition&) = delete;
	DeviceSectorialScanAcquisition& operator=(const DeviceSectorialScanAcquisition&) = delete;
	DeviceSectorialScanAcquisition(DeviceSectorialScanAcquisition&&) = delete;
	DeviceSectorialScanAcquisition& operator=(DeviceSectorialScanAcquisition&&) = delete;
};

} // namespace Lab

#endif // DEVICESECTORIALSCANACQUISITION_H
