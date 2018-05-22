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

#ifndef FILEUTIL_H
#define FILEUTIL_H

#include <string>



namespace Lab {
namespace FileUtil {

std::string path(const std::string& baseDir,
			const std::string& str0, unsigned int index0);
std::string path(const std::string& baseDir,
			const std::string& str0, unsigned int index0,
			const std::string& str1, unsigned int index1);
std::string staSignalPath(const std::string& baseDir, unsigned int baseElement, unsigned int txElem);
std::string signalPath(const std::string& baseDir, unsigned int baseElement, unsigned int txElem, unsigned int rxElem);

} // namespace FileUtil
} // namespace Lab

#endif // FILEUTIL_H
