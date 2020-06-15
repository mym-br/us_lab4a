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

#ifndef NETWORKSYNCSTAMETHOD_H
#define NETWORKSYNCSTAMETHOD_H

#include <string>

#include "Method.h"



namespace Lab {

class Project;
template<typename TFloat> class STAAcquisition;
template<typename TFloat> struct STAConfiguration;

template<typename TFloat>
class NetworkSyncSTAMethod : public Method {
public:
	NetworkSyncSTAMethod(Project& project);
	virtual ~NetworkSyncSTAMethod() = default;

	virtual void execute();
private:
	static constexpr unsigned int maxSteps = 10000;
	static constexpr const char* timeFile    = "/time";
	static constexpr const char* timeDataset = "time";
	static constexpr const char* yFile    = "/y";
	static constexpr const char* yDataset = "y";

	NetworkSyncSTAMethod(const NetworkSyncSTAMethod&) = delete;
	NetworkSyncSTAMethod& operator=(const NetworkSyncSTAMethod&) = delete;
	NetworkSyncSTAMethod(NetworkSyncSTAMethod&&) = delete;
	NetworkSyncSTAMethod& operator=(NetworkSyncSTAMethod&&) = delete;

	void saveSignals(const STAConfiguration<TFloat>& config, STAAcquisition<TFloat>& acq,
				unsigned int baseElement, const std::string& dataDir);

	Project& project_;
};

} // namespace Lab

#endif // NETWORKSYNCSTAMETHOD_H
