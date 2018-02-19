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

#ifndef FILEUTIL_H
#define FILEUTIL_H

#include <fstream>
#include <iterator> /* istreambuf_iterator */
#include <string>

#include "Exception.h"



namespace Lab {
namespace FileUtil {

std::string
loadASCIIFileToString(const char* fileName)
{
	std::ifstream in(fileName);
	if (!in) THROW_EXCEPTION(IOException, "File not found: " << fileName << '.');

	std::string s(
		(std::istreambuf_iterator<char>(in)),
		(std::istreambuf_iterator<char>()));
	if (!in) THROW_EXCEPTION(IOException, "An error occurred while reading the file " << fileName << '.');
	return s;
}

} // namespace FileUtil
} // namespace Lab

#endif // FILEUTIL_H
