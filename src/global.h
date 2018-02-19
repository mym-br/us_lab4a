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

#ifndef GLOBAL_H
#define GLOBAL_H

#define NETWORK_AQUISITION_CONFIG_FILE "config-network_acquisition.txt"

#define SETTINGS_KEY_PROJECT_DIR             "main/projectDir"
#define SETTINGS_KEY_LOGLEVEL_COMBOBOX_INDEX "main/logLevelComboBoxIndex"
#define SETTINGS_KEY_SELECTED_TASK           "main/selectedTask"

#define LOG_FILE_NAME "log-us_lab4a.txt"

#define FLOAT_SCIENTIFIC_NOTATION_NUM_DIGITS_AFTER_DECIMAL_POINT 7



namespace Lab {
namespace Figure {

enum Colormap {
	COLORMAP_GRAY,
	COLORMAP_INVERTED_GRAY,
	COLORMAP_VIRIDIS,
	COLORMAP_INVERTED_VIRIDIS,
	COLORMAP_PLASMA,
	COLORMAP_INVERTED_PLASMA,
	COLORMAP_INFERNO,
	COLORMAP_INVERTED_INFERNO,
	COLORMAP_MAGMA,
	COLORMAP_INVERTED_MAGMA,
	COLORMAP_RED_WHITE_BLUE,
	COLORMAP_DEFAULT // must be the last
};
enum Visualization {
	VISUALIZATION_RAW_LINEAR,
	VISUALIZATION_RECTIFIED_LINEAR,
	VISUALIZATION_RECTIFIED_LOG,
	VISUALIZATION_ENVELOPE_LINEAR,
	VISUALIZATION_ENVELOPE_LOG,
	VISUALIZATION_DEFAULT // must be the last
};

} // namespace Figure
} // namespace Lab

#endif // GLOBAL_H
