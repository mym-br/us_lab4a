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

#ifndef METHOD_H_
#define METHOD_H_

#include <string>
#include <unordered_map>

#include "method_table.h"

namespace Lab {

enum class MethodType {
	invalid
#define METHOD_ITEM(A, B) ,A
	METHOD_TABLE
#undef METHOD_ITEM
};

class Project;

class MethodNameMap {
public:
	MethodNameMap();
	~MethodNameMap();

	MethodType findByName(const std::string& name);
private:
	typedef std::unordered_map<std::string, MethodType> Map;

	Map map_;
};

class Method {
public:
	Method() {}
	virtual ~Method() {}

	virtual void execute();

	static MethodType findByName(const std::string& name) { return nameMap_.findByName(name); }
	static Method* get(Project& project);
private:
	static MethodNameMap nameMap_;
};

} // namespace Lab

#endif /* METHOD_H_ */
