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

#ifndef METHOD_H_
#define METHOD_H_

#include <string>

#include <boost/unordered_map.hpp>



namespace Lab {

enum class MethodType {
	invalid,
	single_acquisition,
	sectorial_scan_sp_network,
	sectorial_scan_sp_network_continuous,
	sectorial_scan_sp_network_trigger,
	sectorial_scan_sp_saved,
	sta_sectorial_simple_simulated,
	sta_sectorial_simple_saved,
	sta_sectorial_simulated,
	sta_sectorial_dp_network,
	sta_sectorial_dp_saved,
	sta_sectorial_vectorial_dp_saved,
	sta_sectorial_sp_saved,
	sta_sectorial_vectorial_sp_saved,
	sta_save_signals,
	sim_acoustic_beam_array_of_rectangular_flat_sources_transient,
	sim_acoustic_beam_rectangular_flat_source_transient,
	sim_acoustic_field_array_of_rectangular_flat_sources_transient,
	sim_acoustic_field_rectangular_flat_source_transient,
	sim_impulse_response_array_of_rectangular_flat_sources,
	sim_impulse_response_rectangular_flat_source,
	show_image,
	test
};

class Project;

class MethodNameMap {
public:
	MethodNameMap();
	~MethodNameMap();

	MethodType findByName(const std::string& name);
private:
	typedef boost::unordered_map<std::string, MethodType> Map;

	Map map_;
};

class Method {
public:
	Method() {}
	virtual ~Method() {}

	virtual void execute() = 0;

	static MethodType findByName(const std::string& name) { return nameMap_.findByName(name); }
	static Method* get(Project& project);
private:
	static MethodNameMap nameMap_;
};

} // namespace Lab

#endif /* METHOD_H_ */
