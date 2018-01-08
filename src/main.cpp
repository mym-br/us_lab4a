#include "USLab4a.h"

#include <fftw3.h>
#include <iostream>

#include "tbb/task_scheduler_init.h"

#include <QApplication>
#include <QCoreApplication>

#ifdef _WIN32
# define WIN32_LEAN_AND_MEAN
# include <windows.h> /* timeBeginPeriod, timeEndPeriod */
#endif

#include "Log.h"
#include "Method.h"

#define FFTW_WISDOM_FILE_NAME_SP "fftw_wisdom_sp-us_lab4a.txt"
#define FFTW_WISDOM_FILE_NAME_DP "fftw_wisdom_dp-us_lab4a.txt"



int
main(int argc, char* argv[])
{
#ifdef _WIN32
	timeBeginPeriod(1); // set the minimum timer resolution to 1 ms
#endif

	fftwf_import_wisdom_from_filename(FFTW_WISDOM_FILE_NAME_SP);
	fftw_import_wisdom_from_filename(FFTW_WISDOM_FILE_NAME_DP);

	//tbb::task_scheduler_init init(2);
	std::cout << "tbb::task_scheduler_init::default_num_threads() = " << tbb::task_scheduler_init::default_num_threads() << std::endl;

	//QApplication::setGraphicsSystem("native"); // "raster" / "native" / "opengl" (experimental in Qt-4.8.4)

	Lab::Log::setLevel(Lab::Log::LEVEL_DEBUG);

	QApplication a(argc, argv);

	// For QSettings.
	QCoreApplication::setOrganizationName("LUS EP USP");
	QCoreApplication::setOrganizationDomain("lus.poli.usp.br");
	QCoreApplication::setApplicationName("USLab4a");

	Lab::Method::fillNameMap();

	Lab::USLab4a w;
	//w.show();
	int returnValue = a.exec();

	fftw_export_wisdom_to_filename(FFTW_WISDOM_FILE_NAME_DP);
	fftwf_export_wisdom_to_filename(FFTW_WISDOM_FILE_NAME_SP);

#ifdef _WIN32
	timeEndPeriod(1);
#endif

	return returnValue;
}
