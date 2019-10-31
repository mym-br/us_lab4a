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

#ifndef METHOD_TABLE_H
#define METHOD_TABLE_H

#define METHOD_TABLE \
METHOD_ITEM(DeviceSectorialScanMethod<float>, device_sectorial_scan_sp_network) \
METHOD_ITEM(DeviceSectorialScanMethod<float>, device_sectorial_scan_sp_network_continuous) \
METHOD_ITEM(DeviceSectorialScanMethod<float>, device_sectorial_scan_sp_network_trigger) \
METHOD_ITEM(DeviceSectorialScanMethod<float>, device_sectorial_scan_sp_saved) \
METHOD_ITEM(MultiLayerImageMethod, multi_layer_image) \
METHOD_ITEM(NetworkSyncSingleVirtualSourceMethod<double>, single_virtual_source_network_sync_imaging) \
METHOD_ITEM(NetworkSyncSingleVirtualSourceMethod<double>, single_virtual_source_network_sync_save_signals) \
METHOD_ITEM(NetworkSyncSTAMethod<double>, sta_network_sync) \
METHOD_ITEM(NetworkSyncSTAMethod<double>, sta_network_sync_save_signals) \
METHOD_ITEM(ShowImageMethod, show_image) \
METHOD_ITEM(SimCircularSourceMethod<double>, sim_acoustic_field_circular_source_transient) \
METHOD_ITEM(SimCircularSourceMethod<double>, sim_impulse_response_circular_source) \
METHOD_ITEM(SimCircularSourceMethod<double>, sim_propagation_circular_source_transient) \
METHOD_ITEM(SimCircularSourceMethod<double>, sim_radiation_pattern_circular_source_transient) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_acoustic_field_array_of_rectangular_flat_sources_transient) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_acoustic_field_rectangular_flat_source_transient) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_impulse_response_array_of_rectangular_flat_sources) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_impulse_response_rectangular_flat_source) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_propagation_array_of_rectangular_flat_sources_transient) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_propagation_rectangular_flat_source_transient) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_radiation_pattern_array_of_rectangular_flat_sources_transient) \
METHOD_ITEM(SimRectangularFlatSourceMethod<double>, sim_radiation_pattern_rectangular_flat_source_transient) \
METHOD_ITEM(SingleAcquisitionMethod, single_acquisition) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_network_save_signals) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_network_save_signal_sequence) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_simulated_save_signals) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_vectorial_dp_network) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_vectorial_dp_saved) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_vectorial_dp_saved_sequence) \
METHOD_ITEM(SingleVirtualSourceMethod<double>, single_virtual_source_3d_vectorial_simulated) \
METHOD_ITEM(SingleVirtualSourceMethod<float>, single_virtual_source_3d_vectorial_sp_network_continuous) \
METHOD_ITEM(STA3DMethod<double>, sta_3d_simulated_save_signals) \
METHOD_ITEM(STA3DMethod<double>, sta_3d_simulated_seq_y_save_signals) \
METHOD_ITEM(STA3DMethod<double>, sta_3d_vectorial_simulated) \
METHOD_ITEM(STAMethod<double>, sta_dp_network) \
METHOD_ITEM(STAMethod<double>, sta_dp_saved) \
METHOD_ITEM(STAMethod<double>, sta_save_signals) \
METHOD_ITEM(STAMethod<double>, sta_simple_saved) \
METHOD_ITEM(STAMethod<double>, sta_simple_simulated) \
METHOD_ITEM(STAMethod<double>, sta_simulated) \
METHOD_ITEM(STAMethod<double>, sta_vectorial_dp_network) \
METHOD_ITEM(STAMethod<double>, sta_vectorial_dp_saved) \
METHOD_ITEM(STAMethod<float>, sta_sp_saved) \
METHOD_ITEM(STAMethod<float>, sta_vectorial_sp_saved) \
METHOD_ITEM(SyntheticYSingleVirtualSourceMethod<double>, single_virtual_source_network_sync_synth_y_imaging) \
METHOD_ITEM(T1R1SAFT3DMethod<double>, t1r1saft_3d_simulated_save_signals) \
METHOD_ITEM(T1R1SAFT3DMethod<double>, t1r1saft_3d_simulated_seq_y_save_signals) \
METHOD_ITEM(T1R1SAFT3DMethod<double>, t1r1saft_3d_vectorial_simulated) \
METHOD_ITEM(TestMethod, test) \
METHOD_ITEM(VTKFileMultiImageMethod, multi_image_vtk_file)

#endif // METHOD_TABLE_H
