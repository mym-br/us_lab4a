/***************************************************************************
 *  Copyright 2024 Marcelo Y. Matuda                                       *
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

#include "ThreadUtil.h"

#include <algorithm> /* min */

#include <tbb/info.h>



namespace Lab {
namespace ThreadUtil {

std::unique_ptr<tbb::global_control>
createNumThreadsControl(unsigned int numThreads)
{
	if (numThreads == 0) return std::unique_ptr<tbb::global_control>();

	const int defConc = tbb::info::default_concurrency();
	const unsigned int boundedNumThreads = (defConc <= 0) ? 1U : std::min(static_cast<unsigned int>(defConc), numThreads);
	return std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, boundedNumThreads);
}

} // namespace ThreadUtil
} // namespace Lab
