/***************************************************************************
 *  Copyright 2018, 2019 Marcelo Y. Matuda                                 *
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

#ifndef NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_H
#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_H

#include <string>
#include <vector>

#include "Method.h"



namespace Lab {

class Project;
template<typename TFloat> class TnRnAcquisition;

template<typename TFloat>
class NetworkSyncSingleVirtualSourceMethod : public Method {
public:
	NetworkSyncSingleVirtualSourceMethod(Project& project);
	virtual ~NetworkSyncSingleVirtualSourceMethod() = default;

	virtual void execute();
private:
	static constexpr unsigned int maxSteps = 10000;
	static constexpr const char* timeFile    = "/time";
	static constexpr const char* timeDataset = "time";
	static constexpr const char* yFile    = "/y";
	static constexpr const char* yDataset = "y";

	NetworkSyncSingleVirtualSourceMethod(const NetworkSyncSingleVirtualSourceMethod&) = delete;
	NetworkSyncSingleVirtualSourceMethod& operator=(const NetworkSyncSingleVirtualSourceMethod&) = delete;
	NetworkSyncSingleVirtualSourceMethod(NetworkSyncSingleVirtualSourceMethod&&) = delete;
	NetworkSyncSingleVirtualSourceMethod& operator=(NetworkSyncSingleVirtualSourceMethod&&) = delete;

	void saveSignals(TnRnAcquisition<TFloat>& acq,
				unsigned int baseElement, const std::vector<TFloat>& txDelays,
				const std::string& dataDir);

	Project& project_;
};

} // namespace Lab

#endif // NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_H
