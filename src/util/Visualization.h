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
#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#define VISUALIZATION_VALUE_TABLE \
VISUALIZATION_VALUE_ITEM(RAW_LINEAR      , "Raw - linear"      ) \
VISUALIZATION_VALUE_ITEM(RECTIFIED_LINEAR, "Rectified - linear") \
VISUALIZATION_VALUE_ITEM(RECTIFIED_LOG   , "Rectified - dB"    ) \
VISUALIZATION_VALUE_ITEM(ENVELOPE_LINEAR , "Envelope - linear" ) \
VISUALIZATION_VALUE_ITEM(ENVELOPE_LOG    , "Envelope - dB"     ) \

namespace Lab {
namespace Visualization {

enum class Value {
#define VISUALIZATION_VALUE_ITEM(A, B) A,
	VISUALIZATION_VALUE_TABLE
#undef VISUALIZATION_VALUE_ITEM
	DEFAULT // must be the last
};

extern const char* valueNameList[];

} // Visualization
} // Lab

#endif // VISUALIZATION_H
