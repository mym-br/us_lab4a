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

#ifndef OCL_REDUCE_H
#define OCL_REDUCE_H

#include <string>



namespace Lab {
namespace OCLReduce {

inline
std::string
code() {
	return R"CLC(

// Call with global size: n * 1024, group size: 1024.
__kernel
void
groupReduceMinMaxKernel(
		__global unsigned int* in,
		__global unsigned int* outMin,
		__global unsigned int* outMax,
		unsigned int n,
		__local unsigned int* groupMinData,
		__local unsigned int* groupMaxData) {

	unsigned int minValue = 0xffffffff;
	unsigned int maxValue = 0;
	const unsigned int totalNumThreads = get_global_size(0);
	for (unsigned int i = get_global_id(0); i < n; i += totalNumThreads) {
		minValue = min(minValue, in[i]);
		maxValue = max(maxValue, in[i]);
	}
	groupMinData[get_local_id(0)] = minValue;
	groupMaxData[get_local_id(0)] = maxValue;

	barrier(CLK_LOCAL_MEM_FENCE);

	if (get_local_id(0) == 0) {
		minValue = groupMinData[0];
		maxValue = groupMaxData[0];
		for (int i = 1; i < get_local_size(0); ++i) {
			minValue = min(minValue, groupMinData[i]);
			maxValue = max(maxValue, groupMaxData[i]);
		}
		outMin[get_group_id(0)] = minValue; // one value for each group
		outMax[get_group_id(0)] = maxValue;
	}
}

// Call with global size: 1, group size: 1.
__kernel
void
reduceMinMaxKernel(
		__global unsigned int* groupMinVal,
		__global unsigned int* groupMaxVal,
		unsigned int numGroups) {

	if (get_global_id(0) != 0) return; // single thread

	unsigned int minValue = 0xffffffff;
	unsigned int maxValue = 0;
	for (int i = 0; i < numGroups; ++i) {
		minValue = min(minValue, groupMinVal[i]);
		maxValue = max(maxValue, groupMaxVal[i]);
	}

	groupMinVal[0] = minValue;
	groupMaxVal[0] = maxValue;
}

)CLC";
}

} // namespace OCLReduce
} // namespace Lab

#endif // OCL_REDUCE_H
