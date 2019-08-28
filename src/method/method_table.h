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
  METHOD_ITEM(single_acquisition, SingleAcquisitionMethod) \
  METHOD_ITEM(device_sectorial_scan_sp_network           , DeviceSectorialScanMethod<float>) \
  METHOD_ITEM(device_sectorial_scan_sp_network_continuous, DeviceSectorialScanMethod<float>) \
  METHOD_ITEM(device_sectorial_scan_sp_network_trigger   , DeviceSectorialScanMethod<float>) \
  METHOD_ITEM(device_sectorial_scan_sp_saved             , DeviceSectorialScanMethod<float>) \
  METHOD_ITEM(sta_simple_simulated    , STAMethod<double>) \
  METHOD_ITEM(sta_simple_saved        , STAMethod<double>) \
  METHOD_ITEM(sta_simulated           , STAMethod<double>) \
  METHOD_ITEM(sta_dp_network          , STAMethod<double>) \
  METHOD_ITEM(sta_vectorial_dp_network, STAMethod<double>) \
  METHOD_ITEM(sta_dp_saved            , STAMethod<double>) \
  METHOD_ITEM(sta_vectorial_dp_saved  , STAMethod<double>) \
  METHOD_ITEM(sta_sp_saved            , STAMethod<float>) \
  METHOD_ITEM(sta_vectorial_sp_saved  , STAMethod<float>) \
  METHOD_ITEM(sta_save_signals        , STAMethod<double>) \
  METHOD_ITEM(sta_network_sync             , NetworkSyncSTAMethod<double>) \
  METHOD_ITEM(sta_network_sync_save_signals, NetworkSyncSTAMethod<double>) \
  METHOD_ITEM(sim_acoustic_field_array_of_rectangular_flat_sources_transient   , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_acoustic_field_rectangular_flat_source_transient             , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_impulse_response_array_of_rectangular_flat_sources           , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_impulse_response_rectangular_flat_source                     , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_propagation_array_of_rectangular_flat_sources_transient      , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_propagation_rectangular_flat_source_transient                , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_radiation_pattern_array_of_rectangular_flat_sources_transient, SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(sim_radiation_pattern_rectangular_flat_source_transient          , SimRectangularFlatSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_network_save_signals           , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_network_save_signal_sequence   , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_simulated_save_signals         , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_vectorial_dp_network           , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_vectorial_dp_saved             , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_vectorial_dp_saved_sequence    , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_vectorial_simulated            , SingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_3d_vectorial_sp_network_continuous, SingleVirtualSourceMethod<float>) \
  METHOD_ITEM(single_virtual_source_network_sync_imaging     , NetworkSyncSingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_network_sync_save_signals, NetworkSyncSingleVirtualSourceMethod<double>) \
  METHOD_ITEM(single_virtual_source_network_sync_synth_y_imaging, SyntheticYSingleVirtualSourceMethod<double>) \
  METHOD_ITEM(sta_3d_simulated_save_signals      , STA3DMethod<double>) \
  METHOD_ITEM(sta_3d_simulated_seq_y_save_signals, STA3DMethod<double>) \
  METHOD_ITEM(sta_3d_vectorial_simulated         , STA3DMethod<double>) \
  METHOD_ITEM(t1r1saft_3d_simulated_save_signals      , T1R1SAFT3DMethod<double>) \
  METHOD_ITEM(t1r1saft_3d_simulated_seq_y_save_signals, T1R1SAFT3DMethod<double>) \
  METHOD_ITEM(t1r1saft_3d_vectorial_simulated         , T1R1SAFT3DMethod<double>) \
  METHOD_ITEM(show_image, ShowImageMethod) \
  METHOD_ITEM(multi_layer_image, MultiLayerImageMethod) \
  METHOD_ITEM(multi_image_vtk_file, VTKFileMultiImageMethod) \
  METHOD_ITEM(test, TestMethod)

#endif // METHOD_TABLE_H
