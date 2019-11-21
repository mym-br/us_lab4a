/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

#ifndef WAVEFRONTOBJFILEWRITER_H
#define WAVEFRONTOBJFILEWRITER_H

#include <fstream>
#include <vector>

#include "Exception.h"
#include "XYZ.h"

namespace Lab {

template<typename FloatType>
class WavefrontObjFileWriter {
public:
	WavefrontObjFileWriter(const char* filePath);

	void addPoint(FloatType x, FloatType y, FloatType z);
	void adjustPointIndex(int& i);
	void addTri(int p0, int p1, int p2); // first index: 0
	void addQuad(int p0, int p1, int p2, int p3); // first index: 0
	void write();
private:
	std::ofstream out_;
	std::vector<XYZ<FloatType>> pointList_;
	std::vector<std::vector<int>> faceList_;
};

template<typename FloatType>
WavefrontObjFileWriter<FloatType>::WavefrontObjFileWriter(const char* filePath)
		: out_(filePath)
{
	out_ << "g group0\n\n";
}

template<typename FloatType>
void WavefrontObjFileWriter<FloatType>::addPoint(FloatType x, FloatType y, FloatType z)
{
	pointList_.push_back(XYZ<FloatType>{x, y, z});
}

template<typename FloatType>
void
WavefrontObjFileWriter<FloatType>::adjustPointIndex(int& i)
{
	if (i < 0) {
		i = static_cast<int>(pointList_.size()) + i;
		if (i < 0) THROW_EXCEPTION(InvalidValueException, "Invalid point index (negative after adjustment).");
	} else if (i >= static_cast<int>(pointList_.size())) {
		THROW_EXCEPTION(InvalidValueException, "Invalid point index.");
	}
}

template<typename FloatType>
void
WavefrontObjFileWriter<FloatType>::addTri(int p0, int p1, int p2)
{
	adjustPointIndex(p0);
	adjustPointIndex(p1);
	adjustPointIndex(p2);
	faceList_.push_back(std::vector<int>{p0, p1, p2});
}

template<typename FloatType>
void
WavefrontObjFileWriter<FloatType>::addQuad(int p0, int p1, int p2, int p3)
{
	adjustPointIndex(p0);
	adjustPointIndex(p1);
	adjustPointIndex(p2);
	adjustPointIndex(p3);
	faceList_.push_back(std::vector<int>{p0, p1, p2, p3});
}

template<typename FloatType>
void
WavefrontObjFileWriter<FloatType>::write()
{
	for (const auto& p : pointList_) {
		out_ << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
	}
	out_ << '\n';

	for (const auto& f : faceList_) {
		out_ << 'f';
		for (const auto& v : f) {
			out_ << ' ' << v + 1;
		}
		out_ << '\n';
	}
	out_ << '\n';
}

} // namespace Lab

#endif // WAVEFRONTOBJFILEWRITER_H
