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

#include "KeyValueFileReader.h"

#include <cstddef> /* std::size_t */
#include <fstream>

#include "Exception.h"
#include "String.h"



namespace Lab {

KeyValueFileReader::KeyValueFileReader(const std::string& filePath)
{
	std::ifstream in(filePath);
	if (!in) {
		THROW_EXCEPTION(IOException, "Could not open the file " << filePath << '.');
	}

	unsigned int lineNumber = 0;
	std::string line;
	while (std::getline(in, line)) {
		++lineNumber;

		std::string_view trimmedLine = String::trimToView(line);
		if (trimmedLine.empty()) continue;
		if (trimmedLine[0] == '#') continue; // comment

		std::size_t eqPos = trimmedLine.find('=');
		if (eqPos == trimmedLine.npos) {
			THROW_EXCEPTION(InvalidValueException, "Missing '=' separator at line " << lineNumber
					<< " of file " << filePath << '.');
		}
		if (eqPos == 0) {
			THROW_EXCEPTION(InvalidValueException, "Missing key at line " << lineNumber
					<< " of file " << filePath << '.');
		}
		if (eqPos == trimmedLine.size() - 1U) {
			THROW_EXCEPTION(InvalidValueException, "Missing value at line " << lineNumber
					<< " of file " << filePath << '.');
		}

		std::string key   = String::trim(trimmedLine.substr(0, eqPos));
		std::string value = String::trim(trimmedLine.substr(eqPos + 1));
		if (String::hasSpace(key)) {
			THROW_EXCEPTION(InvalidValueException, "Key has space at line " << lineNumber
					<< " of file " << filePath << '.');
		}

		if (map_.find(key) != map_.end()) {
			THROW_EXCEPTION(InvalidValueException, "Duplicate key '" << key << "' at line " << lineNumber
					<< " of file " << filePath << '.');
		}
		map_[key] = value;
	}
}

} // namespace Lab
