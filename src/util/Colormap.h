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
#ifndef COLORMAP_H
#define COLORMAP_H

#include "Util.h"

#define COLORMAP_TABLE \
COLORMAP_ITEM(GRADIENT_GRAY            , "Gray"            , Colormap::GrayScale                    ) \
COLORMAP_ITEM(GRADIENT_INVERTED_GRAY   , "Inverted gray"   , Colormap::Inverted<Colormap::GrayScale>) \
COLORMAP_ITEM(GRADIENT_VIRIDIS         , "Viridis"         , Colormap::Viridis                      ) \
COLORMAP_ITEM(GRADIENT_INVERTED_VIRIDIS, "Inverted Viridis", Colormap::Inverted<Colormap::Viridis>  ) \
COLORMAP_ITEM(GRADIENT_CIVIDIS         , "Cividis"         , Colormap::Cividis                      ) \
COLORMAP_ITEM(GRADIENT_INVERTED_CIVIDIS, "Inverted Cividis", Colormap::Inverted<Colormap::Cividis>  ) \
COLORMAP_ITEM(GRADIENT_PLASMA          , "Plasma"          , Colormap::Plasma                       ) \
COLORMAP_ITEM(GRADIENT_INVERTED_PLASMA , "Inverted Plasma" , Colormap::Inverted<Colormap::Plasma>   ) \
COLORMAP_ITEM(GRADIENT_INFERNO         , "Inferno"         , Colormap::Inferno                      ) \
COLORMAP_ITEM(GRADIENT_INVERTED_INFERNO, "Inverted Inferno", Colormap::Inverted<Colormap::Inferno>  ) \
COLORMAP_ITEM(GRADIENT_MAGMA           , "Magma"           , Colormap::Magma                        ) \
COLORMAP_ITEM(GRADIENT_INVERTED_MAGMA  , "Inverted Magma"  , Colormap::Inverted<Colormap::Magma>    ) \
COLORMAP_ITEM(GRADIENT_RED_WHITE_BLUE  , "Red-white-blue"  , Colormap::RedWhiteBlue                 )

namespace Lab {
namespace Colormap {

enum class Id {
#define COLORMAP_ITEM(A, B, C) A,
	COLORMAP_TABLE
#undef COLORMAP_ITEM
	DEFAULT // must be the last
};

extern const char* nameList[];

template<typename Colormap>
struct TableColormap {
	template<typename Point>
	static void setColor(float value, Point& point) {
		Util::clip(value, 0.0f, 1.0f);
		const float pos = value * 255.0f;
		const unsigned int basePos = static_cast<unsigned int>(pos);
		const unsigned int nextPos = basePos + 1U;
		if (nextPos == 256U) {
			float* color = Colormap::table[basePos];
			point.color.red   = color[0];
			point.color.green = color[1];
			point.color.blue  = color[2];
		} else {
			const float coef = pos - basePos;
			float* color0 = Colormap::table[basePos];
			float* color1 = Colormap::table[nextPos];
			point.color.red   = color0[0] + (color1[0] - color0[0]) * coef;
			point.color.green = color0[1] + (color1[1] - color0[1]) * coef;
			point.color.blue  = color0[2] + (color1[2] - color0[2]) * coef;
		}
	}
};

struct GrayScale {
	template<typename Point>
	static void setColor(float value, Point& point) {
		Util::clip(value, 0.0f, 1.0f);
		point.color.red   = value;
		point.color.green = value;
		point.color.blue  = value;
	}
};

struct Viridis : TableColormap<Viridis> {
	static float table[256][3];
};

struct Plasma : TableColormap<Plasma> {
	static float table[256][3];
};

struct Inferno : TableColormap<Inferno> {
	static float table[256][3];
};

struct Magma : TableColormap<Magma> {
	static float table[256][3];
};

struct Cividis : TableColormap<Cividis> {
	static float table[256][3];
};

struct RedWhiteBlue {
	template<typename Point>
	static void setColor(float value, Point& point) {
		Util::clip(value, 0.0f, 1.0f);
		if (value >= 0.5f) {
			const float v = 2.0f - 2.0f * value;
			point.color.red   = v;
			point.color.green = v;
			point.color.blue  = 1.0f;
		} else {
			const float v = 2.0f * value;
			point.color.red   = 1.0f;
			point.color.green = v;
			point.color.blue  = v;
		}
	}
};

template<typename Colormap>
struct Inverted {
	template<typename Point>
	static void setColor(float value, Point& point) {
		Colormap::setColor(1.0f - value, point);
	}
};

} // namespace Colormap
} // namespace Lab

#endif // COLORMAP_H
