/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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

#include <cfenv> /* fesetround */
#include <fftw3.h>
#include <iostream>
#include <string>
#include <thread>

#include "Log.h"
#include "lzf_filter.h"
#include "Method.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Timer.h"



namespace {

constexpr const char* FFTW_WISDOM_FILE_NAME_SP = "fftw_wisdom_sp-us_lab4a.txt";
constexpr const char* FFTW_WISDOM_FILE_NAME_DP = "fftw_wisdom_dp-us_lab4a.txt";

}

namespace Lab {

struct USLab4a {};

}

void
showUsage(const char* name)
{
	std::cout << "\nUsage: " << name << " project_dir task_file [experiment_dir]\n" << std::endl;
}

void
showLogMessages()
{
	std::string logMessages;
	Lab::Log::transferTo(logMessages);
	std::cout << logMessages << std::endl;
}

int
main(int argc, char* argv[])
{
	if (argc != 3 && argc != 4) {
		showUsage(argv[0]);
		return EXIT_FAILURE;
	}

	std::cout << "std::thread::hardware_concurrency(): "
			<< std::thread::hardware_concurrency() << std::endl;

	std::cout << "================================================================================" << std::endl;
	std::string projectDir = argv[1];
	std::cout << "Project directory:    " << projectDir << std::endl;
	std::string taskFile = argv[2];
	std::cout << "Task file:            " << taskFile << std::endl;
	std::string experimentDir = (argc == 4) ? argv[3] : "";
	std::cout << "Experiment directory: " << experimentDir << std::endl;

	std::fesetround(FE_TONEAREST);

	if (register_lzf() < 0) {
		std::cerr << "Could not register the LZF filter (HDF5)." << std::endl;
		return EXIT_FAILURE;
	}

	fftwf_import_wisdom_from_filename(FFTW_WISDOM_FILE_NAME_SP);
	fftw_import_wisdom_from_filename(FFTW_WISDOM_FILE_NAME_DP);

	Lab::Log::setLogToStdOut();
	Lab::Log::setLevel(Lab::Log::LEVEL_DEBUG);

	Lab::USLab4a dummy;
	Lab::Project project(dummy);

	project.setDirectory(projectDir);
	project.setExpDirectory(experimentDir);

	try {
		project.loadTaskParameters(taskFile);
		const Lab::ParameterMap& pm = project.taskParamMap();
		const auto methodName = pm.value<std::string>("method");
		std::cout << "Method:               " << methodName << std::endl;
		std::cout << "================================================================================" << std::endl;
		Lab::MethodEnum methodCode = Lab::Method::findByName(methodName);
		project.setMethod(methodCode);

		auto method = Lab::Method::get(project);

		Lab::Timer procTimer;
		method->execute();
		LOG_INFO << "########## Processing time = " << procTimer.getTime() << " s";

	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	fftw_export_wisdom_to_filename(FFTW_WISDOM_FILE_NAME_DP);
	fftwf_export_wisdom_to_filename(FFTW_WISDOM_FILE_NAME_SP);

	return EXIT_SUCCESS;
}
