/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#include "Method.h"

#include "DeviceSectorialScanMethod.h"
#include "Exception.h"
#include "MultiLayerImageMethod.h"
#include "NetworkSyncSingleVirtualSourceMethod.h"
#include "NetworkSyncSTAMethod.h"
#include "Project.h"
#include "ShowImageMethod.h"
#include "SimRectangularFlatSourceMethod.h"
#include "SingleAcquisitionMethod.h"
#include "SingleVirtualSourceMethod.h"
#include "STA3DMethod.h"
#include "STAMethod.h"
#include "SyntheticYSingleVirtualSourceMethod.h"
#include "T1R1SAFT3DMethod.h"
#include "TestMethod.h"
#include "VTKFileMultiImageMethod.h"

namespace Lab {

MethodNameMap Method::nameMap_;

MethodNameMap::MethodNameMap()
{
#define METHOD_ITEM(A, B) map_[#A] = MethodType::A;
	METHOD_TABLE
#undef METHOD_ITEM
}

MethodNameMap::~MethodNameMap()
{
}

MethodType
MethodNameMap::findByName(const std::string& name)
{
	Map::const_iterator iter = map_.find(name);
	if (iter == map_.end()) {
		THROW_EXCEPTION(InvalidParameterException, "Could not find a method with name \"" << name << "\".");
	}
	return iter->second;
}

void
Method::execute()
{
}

Method*
Method::get(Project& project)
{
	switch (project.method()) {
#define METHOD_ITEM(A, B) case MethodType::A: return new B(project);
	METHOD_TABLE
#undef METHOD_ITEM
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project.method()) << '.');
	}
}

} // namespace Lab
