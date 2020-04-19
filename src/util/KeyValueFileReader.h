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

#ifndef KEYVALUEFILEREADER_H_
#define KEYVALUEFILEREADER_H_

#include <string>
#include <unordered_map>



namespace Lab {

class KeyValueFileReader {
public:
	explicit KeyValueFileReader(const std::string& filePath);
	~KeyValueFileReader() = default;

	const std::unordered_map<std::string, std::string>& map() const { return map_; }
private:
	KeyValueFileReader(const KeyValueFileReader&) = delete;
	KeyValueFileReader& operator=(const KeyValueFileReader&) = delete;
	KeyValueFileReader(KeyValueFileReader&&) = delete;
	KeyValueFileReader& operator=(KeyValueFileReader&&) = delete;

	std::unordered_map<std::string, std::string> map_;
};

} // namespace Lab

#endif /* KEYVALUEFILEREADER_H_ */
