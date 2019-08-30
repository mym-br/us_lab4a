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

#ifndef ARRAY_H_
#define ARRAY_H_

#include <cstddef> /* std::size_t */

namespace Lab {

template<typename T, std::size_t N>
class Array {
public:
	typedef T value_type;
	typedef std::size_t size_type;

	Array() {}
	Array(const T& value)
	{
		fill(value);
	}
	~Array() {}

	T& operator[](size_type i)
	{
		return data_[i];
	}

	const T& operator[](size_type i) const
	{
		return data_[i];
	}

	static size_type size() { return N; }

	void fill(const T& value)
	{
		for (size_type i = 0; i < N; ++i) {
			data_[i] = value;
		}
	}
private:
	Array(const Array&);
	Array& operator=(const Array&);

	T data_[N];
};

} // namespace Lab

#endif /* ARRAY_H_ */
