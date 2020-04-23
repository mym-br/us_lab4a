/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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

#ifndef TEMPLATE_UTIL_H_
#define TEMPLATE_UTIL_H_

#include <string>
#include <typeinfo>
#include <type_traits>



namespace Lab {

// Primary template.
template<typename, typename =void>
struct has_y_member : std::false_type { };
// Partial template specialization.
template<typename T>
struct has_y_member<T, std::void_t<decltype(T::y)>> : std::true_type { };



namespace TemplateUtil {

template<typename T, std::enable_if_t<std::is_same_v<T, float>, int> = 0>
std::string
typeName()
{
	return "float";
}

template<typename T, std::enable_if_t<std::is_same_v<T, double>, int> = 0>
std::string
typeName()
{
	return "double";
}

} // namespace TemplateUtil
} // namespace Lab

#endif /*TEMPLATE_UTIL_H_*/
