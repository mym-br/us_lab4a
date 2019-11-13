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

constexpr const char* NETWORK_AQUISITION_CONFIG_FILE = "config-network_acquisition.txt";

constexpr const char* SETTINGS_KEY_PROJECT_DIR             = "main/projectDir";
constexpr const char* SETTINGS_KEY_LOGLEVEL_COMBOBOX_INDEX = "main/logLevelComboBoxIndex";
constexpr const char* SETTINGS_KEY_SELECTED_TASK           = "main/selectedTask";
constexpr const char* SETTINGS_KEY_SELECTED_EXP            = "main/selectedExp";

constexpr const char* LOG_FILE_NAME = "log-us_lab4a.txt";

constexpr unsigned int SCIENTIFIC_NOTATION_NUM_DIGITS_AFTER_DECIMAL_POINT = 7;



namespace Lab {
namespace Figure {

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
