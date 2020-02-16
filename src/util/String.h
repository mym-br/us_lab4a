/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#ifndef LAB_STRING_H
#define LAB_STRING_H

#include <sstream>
#include <string>



namespace Lab {
namespace String {

struct End {
};

class Begin {
public:
	template<typename T>
	Begin& operator<<(const T& v) {
		out_ << v;
		return *this;
	}

	std::string operator<<(const End& /*end*/) {
		return out_.str();
	}
private:
	std::ostringstream out_;
};

} // namespace String
} // namespace Lab

#endif // LAB_STRING_H
