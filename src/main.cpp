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

#include "USLab4a.h"

#include <cfenv> /* fesetround */
#include <fftw3.h>
#include <iostream>
#include <locale>
#include <thread>

#include <QApplication>
#include <QCoreApplication>
#include <QIcon>
#include <QLocale>

#include "Log.h"
#include "lzf_filter.h"

namespace {

constexpr const char* FFTW_WISDOM_FILE_NAME_SP = "fftw_wisdom_sp-us_lab4a.txt";
constexpr const char* FFTW_WISDOM_FILE_NAME_DP = "fftw_wisdom_dp-us_lab4a.txt";

}

int
main(int argc, char* argv[])
{
	std::fesetround(FE_TONEAREST);

	if (register_lzf() < 0) {
		std::cerr << "Could not register the LZF filter (HDF5)." << std::endl;
		return EXIT_FAILURE;
	}

	fftwf_import_wisdom_from_filename(FFTW_WISDOM_FILE_NAME_SP);
	fftw_import_wisdom_from_filename(FFTW_WISDOM_FILE_NAME_DP);

	std::cout << "std::thread::hardware_concurrency(): "
			<< std::thread::hardware_concurrency() << std::endl;

	Lab::Log::setLevel(Lab::Log::LEVEL_DEBUG);

	QApplication a(argc, argv);
	a.setWindowIcon(QIcon{":/img/window_icon.png"});

	// For QSettings.
	QCoreApplication::setOrganizationName("LUS EP USP");
	QCoreApplication::setOrganizationDomain("lus.poli.usp.br");
	QCoreApplication::setApplicationName("USLab4a");

	// Force "C" locale.
	QLocale::setDefault(QLocale::c());
	std::locale::global(std::locale::classic());

	Lab::USLab4a w;
	int returnValue = a.exec();

	fftw_export_wisdom_to_filename(FFTW_WISDOM_FILE_NAME_DP);
	fftwf_export_wisdom_to_filename(FFTW_WISDOM_FILE_NAME_SP);

	return returnValue;
}
