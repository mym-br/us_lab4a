/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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

#include <iomanip>
#include <sstream>



namespace Lab {
namespace FileUtil {

std::string
path(const std::string& baseDir,
	const std::string& str0, unsigned int index0)
{
	std::ostringstream out;
	out << baseDir << std::setfill('0')
		<< str0 << std::setw(4) << index0;
	return out.str();
}

std::string
path(const std::string& baseDir,
	const std::string& str0, unsigned int index0,
	const std::string& str1, unsigned int index1)
{
	std::ostringstream out;
	out << baseDir << std::setfill('0')
		<< str0 << std::setw(4) << index0
		<< str1 << std::setw(4) << index1;
	return out.str();
}

std::string
staSignalPath(const std::string& baseDir, unsigned int baseElement, unsigned int txElem)
{
	std::ostringstream out;
	out << baseDir << std::setfill('0')
		<< "/signal-base" << std::setw(4) << baseElement
		<< "-tx"          << std::setw(4) << txElem;
	return out.str();
}

std::string
signalPath(const std::string& baseDir, unsigned int baseElement, unsigned int txElem, unsigned int rxElem)
{
	std::ostringstream out;
	out << baseDir << std::setfill('0')
		<< "/signal-base" << std::setw(4) << baseElement
		<< "-tx"          << std::setw(4) << txElem
		<< "-rx"          << std::setw(4) << rxElem;
	return out.str();
}

} // namespace FileUtil
} // namespace Lab
