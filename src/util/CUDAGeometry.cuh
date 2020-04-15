/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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

#ifndef CUDA_GEOMETRY_CUH
#define CUDA_GEOMETRY_CUH

namespace Lab {

__device__
inline
float
distance2D(float x0, float y0, float x1, float y1)
{
	const float dx = x1 - x0;
	const float dy = y1 - y0;
	return sqrtf(dx * dx + dy * dy);
}

__device__
inline
float
distance2DY0(float x0, float x1, float y1)
{
	const float dx = x1 - x0;
	return sqrtf(dx * dx + y1 * y1);
}

} // namespace Lab

#endif // CUDA_GEOMETRY_CUH
