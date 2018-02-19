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

#include "Method.h"

#include "Exception.h"
#include "Project.h"
#include "SectorialScanMethod.h"
#include "ShowImageMethod.h"
#include "SingleAcquisitionMethod.h"
#include "STAMethod.h"
#include "TestMethod.h"

#define ADD_MAP_ITEM(A) map_[#A] = MethodType::A



namespace Lab {

MethodNameMap Method::nameMap_;

MethodNameMap::MethodNameMap()
{
	ADD_MAP_ITEM(single_acquisition);
	ADD_MAP_ITEM(sectorial_scan_sp_network);
	ADD_MAP_ITEM(sectorial_scan_sp_network_continuous);
	ADD_MAP_ITEM(sectorial_scan_sp_network_trigger);
	ADD_MAP_ITEM(sectorial_scan_sp_saved);
	ADD_MAP_ITEM(sta_sectorial_simple_simulated);
	ADD_MAP_ITEM(sta_sectorial_simple_saved);
	ADD_MAP_ITEM(sta_sectorial_simple_simulated);
	ADD_MAP_ITEM(sta_sectorial_dp_network);
	ADD_MAP_ITEM(sta_sectorial_dp_saved);
	ADD_MAP_ITEM(sta_sectorial_vectorial_dp_saved);
	ADD_MAP_ITEM(sta_sectorial_sp_saved);
	ADD_MAP_ITEM(sta_sectorial_vectorial_sp_saved);
	ADD_MAP_ITEM(sta_save_signals);
	ADD_MAP_ITEM(show_image);
	ADD_MAP_ITEM(test);
}

MethodNameMap::~MethodNameMap()
{
}

MethodType
MethodNameMap::findByName(const std::string& name)
{
	Map::const_iterator iter = map_.find(name);
	if (iter == map_.end()) {
		THROW_EXCEPTION(InvalidParameterException, "Could not find a method with name \"" << name << "\".");
	}
	return iter->second;
}

Method*
Method::get(Project& project)
{
	switch (project.method()) {
	case MethodType::single_acquisition:
		return new SingleAcquisitionMethod(project);
	case MethodType::sectorial_scan_sp_network:            // falls through
	case MethodType::sectorial_scan_sp_network_continuous: // falls through
	case MethodType::sectorial_scan_sp_network_trigger:    // falls through
	case MethodType::sectorial_scan_sp_saved:
		return new SectorialScanMethod<float>(project);
	case MethodType::sta_sectorial_simple_simulated:   // falls through
	case MethodType::sta_sectorial_simple_saved:       // falls through
	case MethodType::sta_sectorial_simulated:          // falls through
	case MethodType::sta_sectorial_dp_network:         // falls through
	case MethodType::sta_sectorial_dp_saved:           // falls through
	case MethodType::sta_sectorial_vectorial_dp_saved: // falls through
	case MethodType::sta_save_signals:
		return new STAMethod<double>(project);
	case MethodType::sta_sectorial_sp_saved:           // falls through
	case MethodType::sta_sectorial_vectorial_sp_saved:
		return new STAMethod<float>(project);
	case MethodType::show_image:
		return new ShowImageMethod(project);
	case MethodType::test:
		return new TestMethod(project);
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project.method()) << '.');
	}
}

} // namespace Lab
