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

#include "ParameterMap.h"

#include <limits>



namespace Lab {

ParameterMap::ParameterMap(const std::string& filePath)
		: filePath_(filePath)
		, reader_(filePath)
{
}

bool
ParameterMap::contains(const char* key) const
{
	const auto& map = reader_.map();
	return map.find(key) != map.end();
}

template<>
short
ParameterMap::convertValue<short>(const std::string& s)
{
	std::size_t idx;
	int v;
	try {
		v = std::stoi(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid short integer: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid short integer: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid short integer: " << s << " (garbage at the end).");
	}
	if (v < std::numeric_limits<short>::min() || v > std::numeric_limits<short>::max()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid short integer: " << s << " (out of range).");
	}
	return static_cast<short>(v);
}

template<>
unsigned short
ParameterMap::convertValue<unsigned short>(const std::string& s)
{
	std::size_t idx;
	int v;
	try {
		v = std::stoi(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned short integer: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned short integer: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned short integer: " << s << " (garbage at the end).");
	}
	if (v < 0 || v > std::numeric_limits<short>::max()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned short integer: " << s << " (out of range).");
	}
	return static_cast<unsigned short>(v);
}

template<>
int
ParameterMap::convertValue<int>(const std::string& s)
{
	std::size_t idx;
	int v;
	try {
		v = std::stoi(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid integer: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid integer: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid integer: " << s << " (garbage at the end).");
	}
	return v;
}

// Range: [0, std::numeric_limits<int>::max()]
template<>
unsigned int
ParameterMap::convertValue<unsigned int>(const std::string& s)
{
	std::size_t idx;
	int v;
	try {
		v = std::stoi(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned integer: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned integer: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned integer: " << s << " (garbage at the end).");
	}
	if (v < 0) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned integer: " << s << " (negative).");
	}
	return static_cast<unsigned int>(v);
}

template<>
long
ParameterMap::convertValue<long>(const std::string& s)
{
	std::size_t idx;
	long v;
	try {
		v = std::stol(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid long integer: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid long integer: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid long integer: " << s << " (garbage at the end).");
	}
	return v;
}

// Range: [0, std::numeric_limits<long>::max()]
template<>
unsigned long
ParameterMap::convertValue<unsigned long>(const std::string& s)
{
	std::size_t idx;
	long v;
	try {
		v = std::stol(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned long integer: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned long integer: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned long integer: " << s << " (garbage at the end).");
	}
	if (v < 0) {
		THROW_EXCEPTION(InvalidValueException, "Invalid unsigned long integer: " << s << " (negative).");
	}
	return static_cast<unsigned long>(v);
}

template<>
float
ParameterMap::convertValue<float>(const std::string& s)
{
	std::size_t idx;
	float v;
	try {
		v = std::stof(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid float: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid float: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid float: " << s << " (garbage at the end).");
	}
	return v;
}

template<>
double
ParameterMap::convertValue<double>(const std::string& s)
{
	std::size_t idx;
	double v;
	try {
		v = std::stod(s, &idx);
	} catch (std::invalid_argument&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid double: " << s << " (wrong format).");
	} catch (std::out_of_range&) {
		THROW_EXCEPTION(InvalidValueException, "Invalid double: " << s << " (out of range).");
	}
	if (idx != s.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid double: " << s << " (garbage at the end).");
	}
	return v;
}

template<>
bool
ParameterMap::convertValue<bool>(const std::string& s)
{
	if (s == "true") return true;
	if (s == "false") return false;
	THROW_EXCEPTION(InvalidValueException, "Invalid bool: " << s << '.');
}

template<>
std::string
ParameterMap::convertValue<std::string>(const std::string& s)
{
	return s;
}

} // namespace Lab
