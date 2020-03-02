/***************************************************************************
 *  Copyright 2019, 2020 Marcelo Y. Matuda                                 *
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

#ifdef USE_OPENCL
# define OCL_METHOD_TABLE \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_ocl_sp)
#else
# define OCL_METHOD_TABLE
#endif

#define METHOD_TABLE \
OCL_METHOD_TABLE \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_acquisition_simulated_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_acquisition_network_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_vectorial_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_cc_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_ccbf_pulse_echo_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_ccbf_pitch_catch_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_arc_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_point_detection_tangent_curve_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_circle_fitting_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_two_medium_imaging_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_two_medium_imaging_vectorial_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_two_medium_imaging_combined_sta_dp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_speed_1_measurement) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<double>, cylinder_detection_and_fermat_distance_measurement) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_acquisition_simulated_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_acquisition_network_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_vectorial_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_cc_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_ccbf_pulse_echo_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_ccbf_pitch_catch_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_arc_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_point_detection_tangent_curve_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_circle_fitting_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_two_medium_imaging_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_two_medium_imaging_vectorial_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_two_medium_imaging_combined_cyl_wave_sp) \
METHOD_ITEM(CylinderDetectionAndFermatMethod<float>, cylinder_detection_and_fermat_two_medium_imaging_combined_sta_sp) \
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
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_acoustic_field_array_of_rectangular_sources_transient) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_acoustic_field_rectangular_source_transient) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_impulse_response_array_of_rectangular_sources) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_impulse_response_rectangular_source) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_propagation_array_of_rectangular_sources_transient) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_propagation_rectangular_source_transient) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_radiation_pattern_array_of_rectangular_sources_transient) \
METHOD_ITEM(SimRectangularSourceMethod<double>, sim_radiation_pattern_rectangular_source_transient) \
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
