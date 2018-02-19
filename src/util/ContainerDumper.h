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

#ifndef CONTAINERDUMPER_H_
#define CONTAINERDUMPER_H_

#include <fstream>



namespace Lab {

class ContainerDumper {
public:
	template<typename InputIterator> static void save(const char* fileName, InputIterator first, InputIterator last);
private:
	ContainerDumper();
};

template<typename InputIterator>
void
ContainerDumper::save(const char* fileName, InputIterator first, InputIterator last)
{
	std::ofstream out(fileName);
	while (first != last) {
		out << *first++ << '\n';
	}
}

} // namespace Lab

#endif /* CONTAINERDUMPER_H_ */
