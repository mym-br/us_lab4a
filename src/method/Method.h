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

#ifndef METHOD_H_
#define METHOD_H_

#include <string>
#include <unordered_map>



namespace Lab {

enum class MethodType {
	invalid,
	single_acquisition,
	device_sectorial_scan_sp_network,
	device_sectorial_scan_sp_network_continuous,
	device_sectorial_scan_sp_network_trigger,
	device_sectorial_scan_sp_saved,
	sta_simple_simulated,
	sta_simple_saved,
	sta_simulated,
	sta_dp_network,
	sta_vectorial_dp_network,
	sta_dp_saved,
	sta_vectorial_dp_saved,
	sta_sp_saved,
	sta_vectorial_sp_saved,
	sta_save_signals,
	sta_network_sync,
	sta_network_sync_save_signals,
	sim_acoustic_field_array_of_rectangular_flat_sources_transient,
	sim_acoustic_field_rectangular_flat_source_transient,
	sim_impulse_response_array_of_rectangular_flat_sources,
	sim_impulse_response_rectangular_flat_source,
	sim_propagation_array_of_rectangular_flat_sources_transient,
	sim_propagation_rectangular_flat_source_transient,
	sim_radiation_pattern_array_of_rectangular_flat_sources_transient,
	sim_radiation_pattern_rectangular_flat_source_transient,
	single_virtual_source_3d_network_save_signals,
	single_virtual_source_3d_network_save_signal_sequence,
	single_virtual_source_3d_simulated_save_signals,
	single_virtual_source_3d_vectorial_dp_network,
	single_virtual_source_3d_vectorial_dp_saved,
	single_virtual_source_3d_vectorial_dp_saved_sequence,
	single_virtual_source_3d_vectorial_simulated,
	single_virtual_source_3d_vectorial_sp_network_continuous,
	single_virtual_source_network_sync_imaging,
	single_virtual_source_network_sync_save_signals,
	single_virtual_source_network_sync_synth_y_imaging,
	sta_3d_simulated_save_signals,
	sta_3d_simulated_seq_y_save_signals,
	sta_3d_vectorial_simulated,
	t1r1saft_3d_simulated_save_signals,
	t1r1saft_3d_simulated_seq_y_save_signals,
	t1r1saft_3d_vectorial_simulated,
	show_image,
	multi_layer_image,
	multi_image_vtk_file,
	test
};

class Project;

class MethodNameMap {
public:
	MethodNameMap();
	~MethodNameMap();

	MethodType findByName(const std::string& name);
private:
	typedef std::unordered_map<std::string, MethodType> Map;

	Map map_;
};

class Method {
public:
	Method() {}
	virtual ~Method() {}

	virtual void execute();

	static MethodType findByName(const std::string& name) { return nameMap_.findByName(name); }
	static Method* get(Project& project);
private:
	static MethodNameMap nameMap_;
};

} // namespace Lab

#endif /* METHOD_H_ */
