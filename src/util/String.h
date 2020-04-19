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

#include <cstddef> /* std::size_t */
#include <sstream>
#include <string>
#include <string_view>



namespace Lab {
namespace String {

constexpr const char* spaces = " \t";

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

inline
std::string
trim(const std::string_view& s)
{
	const std::size_t pos1 = s.find_first_not_of(spaces);
	if (pos1 == s.npos) return std::string();

	const std::size_t pos2 = s.find_last_not_of(spaces);
	return std::string(s.data() + pos1, pos2 - pos1 + 1U);
}

inline
std::string
trim(const std::string& s)
{
	const std::size_t pos1 = s.find_first_not_of(spaces);
	if (pos1 == s.npos) return std::string();

	const std::size_t pos2 = s.find_last_not_of(spaces);
	return std::string(s.data() + pos1, pos2 - pos1 + 1U);
}

inline
std::string_view
trimToView(const std::string_view& s)
{
	const std::size_t pos1 = s.find_first_not_of(spaces);
	if (pos1 == s.npos) return std::string_view();

	const std::size_t pos2 = s.find_last_not_of(spaces);
	return std::string_view(s.data() + pos1, pos2 - pos1 + 1U);
}

inline
std::string_view
trimToView(const std::string& s)
{
	const std::size_t pos1 = s.find_first_not_of(spaces);
	if (pos1 == s.npos) return std::string_view();

	const std::size_t pos2 = s.find_last_not_of(spaces);
	return std::string_view(s.data() + pos1, pos2 - pos1 + 1U);
}

inline
bool
hasSpace(const std::string_view& s)
{
	return s.find_first_of(spaces) != s.npos;
}

inline
bool
hasSpace(const std::string& s)
{
	return s.find_first_of(spaces) != s.npos;
}

} // namespace String
} // namespace Lab

#endif // LAB_STRING_H
