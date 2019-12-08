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

#include "method_includes.h"

namespace Lab {

MethodNameMap Method::nameMap_;

MethodNameMap::MethodNameMap()
{
#define METHOD_ITEM(A, B) map_[#B] = MethodEnum::B;
	METHOD_TABLE
#undef METHOD_ITEM
}

MethodEnum
MethodNameMap::findByName(const std::string& name)
{
	Map::const_iterator iter = map_.find(name);
	if (iter == map_.end()) {
		THROW_EXCEPTION(InvalidParameterException, "Could not find a method with name \"" << name << "\".");
	}
	return iter->second;
}

std::unique_ptr<Method>
Method::get(Project& project)
{
	switch (project.method()) {
#define METHOD_ITEM(A, B) case MethodEnum::B: return std::make_unique<A>(project);
	METHOD_TABLE
#undef METHOD_ITEM
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project.method()) << '.');
	}
}

} // namespace Lab
