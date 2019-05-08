/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#include "DeviceSectorialScanMethod.h"
#include "Exception.h"
#include "MultiLayerImageMethod.h"
#include "NetworkSyncSingleVirtualSourceMethod.h"
#include "NetworkSyncSTAMethod.h"
#include "Project.h"
#include "ShowImageMethod.h"
#include "SimRectangularFlatSourceMethod.h"
#include "SingleAcquisitionMethod.h"
#include "SingleVirtualSourceMethod.h"
#include "STA3DMethod.h"
#include "STAMethod.h"
#include "SyntheticYSingleVirtualSourceMethod.h"
#include "T1R1SAFT3DMethod.h"
#include "TestMethod.h"
#include "VTKFileMultiImageMethod.h"

#define ADD_MAP_ITEM(A) map_[#A] = MethodType::A



namespace Lab {

MethodNameMap Method::nameMap_;

MethodNameMap::MethodNameMap()
{
	ADD_MAP_ITEM(single_acquisition);
	ADD_MAP_ITEM(device_sectorial_scan_sp_network);
	ADD_MAP_ITEM(device_sectorial_scan_sp_network_continuous);
	ADD_MAP_ITEM(device_sectorial_scan_sp_network_trigger);
	ADD_MAP_ITEM(device_sectorial_scan_sp_saved);
	ADD_MAP_ITEM(sta_simple_saved);
	ADD_MAP_ITEM(sta_simple_simulated);
	ADD_MAP_ITEM(sta_simulated);
	ADD_MAP_ITEM(sta_dp_network);
	ADD_MAP_ITEM(sta_vectorial_dp_network);
	ADD_MAP_ITEM(sta_dp_saved);
	ADD_MAP_ITEM(sta_vectorial_dp_saved);
	ADD_MAP_ITEM(sta_sp_saved);
	ADD_MAP_ITEM(sta_vectorial_sp_saved);
	ADD_MAP_ITEM(sta_save_signals);
	ADD_MAP_ITEM(sta_network_sync);
	ADD_MAP_ITEM(sta_network_sync_save_signals);
	ADD_MAP_ITEM(sim_acoustic_field_array_of_rectangular_flat_sources_transient);
	ADD_MAP_ITEM(sim_acoustic_field_rectangular_flat_source_transient);
	ADD_MAP_ITEM(sim_impulse_response_array_of_rectangular_flat_sources);
	ADD_MAP_ITEM(sim_impulse_response_rectangular_flat_source);
	ADD_MAP_ITEM(sim_propagation_array_of_rectangular_flat_sources_transient);
	ADD_MAP_ITEM(sim_propagation_rectangular_flat_source_transient);
	ADD_MAP_ITEM(sim_radiation_pattern_array_of_rectangular_flat_sources_transient);
	ADD_MAP_ITEM(sim_radiation_pattern_rectangular_flat_source_transient);
	ADD_MAP_ITEM(single_virtual_source_3d_network_save_signals);
	ADD_MAP_ITEM(single_virtual_source_3d_network_save_signal_sequence);
	ADD_MAP_ITEM(single_virtual_source_3d_simulated_save_signals);
	ADD_MAP_ITEM(single_virtual_source_3d_vectorial_dp_network);
	ADD_MAP_ITEM(single_virtual_source_3d_vectorial_dp_saved);
	ADD_MAP_ITEM(single_virtual_source_3d_vectorial_dp_saved_sequence);
	ADD_MAP_ITEM(single_virtual_source_3d_vectorial_simulated);
	ADD_MAP_ITEM(single_virtual_source_3d_vectorial_sp_network_continuous);
	ADD_MAP_ITEM(single_virtual_source_network_sync_imaging);
	ADD_MAP_ITEM(single_virtual_source_network_sync_save_signals);
	ADD_MAP_ITEM(single_virtual_source_network_sync_synth_y_imaging);
	ADD_MAP_ITEM(sta_3d_simulated_save_signals);
	ADD_MAP_ITEM(sta_3d_simulated_seq_y_save_signals);
	ADD_MAP_ITEM(sta_3d_vectorial_simulated);
	ADD_MAP_ITEM(t1r1saft_3d_simulated_save_signals);
	ADD_MAP_ITEM(t1r1saft_3d_simulated_seq_y_save_signals);
	ADD_MAP_ITEM(t1r1saft_3d_vectorial_simulated);
	ADD_MAP_ITEM(show_image);
	ADD_MAP_ITEM(multi_layer_image);
	ADD_MAP_ITEM(multi_image_vtk_file);
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
	case MethodType::device_sectorial_scan_sp_network:            // falls through
	case MethodType::device_sectorial_scan_sp_network_continuous: // falls through
	case MethodType::device_sectorial_scan_sp_network_trigger:    // falls through
	case MethodType::device_sectorial_scan_sp_saved:
		return new DeviceSectorialScanMethod<float>(project);
	case MethodType::sta_simple_simulated:                // falls through
	case MethodType::sta_simple_saved:                    // falls through
	case MethodType::sta_simulated:                       // falls through
	case MethodType::sta_dp_network:                      // falls through
	case MethodType::sta_vectorial_dp_network:            // falls through
	case MethodType::sta_dp_saved:                        // falls through
	case MethodType::sta_vectorial_dp_saved:              // falls through
	case MethodType::sta_save_signals:
		return new STAMethod<double>(project);
	case MethodType::sta_sp_saved:           // falls through
	case MethodType::sta_vectorial_sp_saved:
		return new STAMethod<float>(project);
	case MethodType::sta_3d_simulated_save_signals:       // falls through
	case MethodType::sta_3d_simulated_seq_y_save_signals: // falls through
	case MethodType::sta_3d_vectorial_simulated:
		return new STA3DMethod<double>(project);
	case MethodType::t1r1saft_3d_simulated_save_signals:       // falls through
	case MethodType::t1r1saft_3d_simulated_seq_y_save_signals: // falls through
	case MethodType::t1r1saft_3d_vectorial_simulated:
		return new T1R1SAFT3DMethod<double>(project);
	case MethodType::sta_network_sync:              // falls through
	case MethodType::sta_network_sync_save_signals:
		return new NetworkSyncSTAMethod<double>(project);
	case MethodType::sim_acoustic_field_array_of_rectangular_flat_sources_transient:       // falls through
	case MethodType::sim_acoustic_field_rectangular_flat_source_transient:                 // falls through
	case MethodType::sim_impulse_response_array_of_rectangular_flat_sources:               // falls through
	case MethodType::sim_impulse_response_rectangular_flat_source:                         // falls through
	case MethodType::sim_propagation_array_of_rectangular_flat_sources_transient:          // falls through
	case MethodType::sim_propagation_rectangular_flat_source_transient:                    // falls through
	case MethodType::sim_radiation_pattern_array_of_rectangular_flat_sources_transient:    // falls through
	case MethodType::sim_radiation_pattern_rectangular_flat_source_transient:
		return new SimRectangularFlatSourceMethod<double>(project);
	case MethodType::single_virtual_source_3d_network_save_signals:         // falls through
	case MethodType::single_virtual_source_3d_network_save_signal_sequence: // falls through
	case MethodType::single_virtual_source_3d_simulated_save_signals:       // falls through
	case MethodType::single_virtual_source_3d_vectorial_dp_network:         // falls through
	case MethodType::single_virtual_source_3d_vectorial_dp_saved:           // falls through
	case MethodType::single_virtual_source_3d_vectorial_dp_saved_sequence:  // falls through
	case MethodType::single_virtual_source_3d_vectorial_simulated:
		return new SingleVirtualSourceMethod<double>(project);
	case MethodType::single_virtual_source_3d_vectorial_sp_network_continuous:
		return new SingleVirtualSourceMethod<float>(project);
	case MethodType::single_virtual_source_network_sync_imaging:      // falls through
	case MethodType::single_virtual_source_network_sync_save_signals:
		return new NetworkSyncSingleVirtualSourceMethod<double>(project);
	case MethodType::single_virtual_source_network_sync_synth_y_imaging:
		return new SyntheticYSingleVirtualSourceMethod<double>(project);
	case MethodType::show_image:
		return new ShowImageMethod(project);
	case MethodType::multi_layer_image:
		return new MultiLayerImageMethod(project);
	case MethodType::multi_image_vtk_file:
		return new VTKFileMultiImageMethod(project);
	case MethodType::test:
		return new TestMethod(project);
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project.method()) << '.');
	}
}

} // namespace Lab
