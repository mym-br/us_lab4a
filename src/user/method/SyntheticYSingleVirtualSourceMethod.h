/***************************************************************************
 *  Copyright 2018, 2019 Marcelo Y. Matuda                                 *
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

#ifndef SYNTHETIC_Y_SINGLE_VIRTUAL_SOURCE_METHOD_H
#define SYNTHETIC_Y_SINGLE_VIRTUAL_SOURCE_METHOD_H

#include "Method.h"



namespace Lab {

class Project;

template<typename TFloat>
class SyntheticYSingleVirtualSourceMethod : public Method {
public:
	SyntheticYSingleVirtualSourceMethod(Project& project);
	virtual ~SyntheticYSingleVirtualSourceMethod() = default;

	virtual void execute();
private:
	static constexpr const char* yFile    = "/y";
	static constexpr const char* yDataset = "y";

	SyntheticYSingleVirtualSourceMethod(const SyntheticYSingleVirtualSourceMethod&) = delete;
	SyntheticYSingleVirtualSourceMethod& operator=(const SyntheticYSingleVirtualSourceMethod&) = delete;
	SyntheticYSingleVirtualSourceMethod(SyntheticYSingleVirtualSourceMethod&&) = delete;
	SyntheticYSingleVirtualSourceMethod& operator=(SyntheticYSingleVirtualSourceMethod&&) = delete;

	Project& project_;
};

} // namespace Lab

#endif // SYNTHETIC_Y_SINGLE_VIRTUAL_SOURCE_METHOD_H
