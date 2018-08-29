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

#ifndef OGL_H
#define OGL_H

namespace Lab {

struct OGLPos3D {
	float x;
	float y;
	float z;
	OGLPos3D(float ix, float iy, float iz) : x(ix), y(iy), z(iz) {}
};

struct OGLColor {
	float red;
	float green;
	float blue;
	OGLColor(float r, float g, float b) : red(r), green(g), blue(b) {}
};

struct OGLColorA {
	float red;
	float green;
	float blue;
	float alpha;
	OGLColorA(float r, float g, float b, float a) : red(r), green(g), blue(b), alpha(a) {}
};

struct OGLPoint3D {
	OGLPos3D pos;
	OGLColor color;
	OGLPoint3D() : pos{0.0, 0.0, 0.0}, color{0.0, 0.0, 0.0} {}
	OGLPoint3D(float x, float y, float z, float r, float g, float b)
		: pos{x, y, z}
		, color{r, g, b} {}
};

struct OGLPoint3DA {
	OGLPos3D pos;
	OGLColorA color;
	OGLPoint3DA() : pos{0.0, 0.0, 0.0}, color{0.0, 0.0, 0.0, 0.0} {}
	OGLPoint3DA(float x, float y, float z, float r, float g, float b, float a)
		: pos{x, y, z}
		, color{r, g, b, a} {}
};

} // namespace Lab

#endif // OGL_H
