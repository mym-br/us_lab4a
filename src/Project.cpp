/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#include "Project.h"

#include <filesystem>
#include <fstream>
#include <string>

#include <QProcess>
#include <QStringList>

#include "KeyValueFileReader.h"
#include "USLab4a.h"

namespace fs = std::filesystem;

namespace Lab {

Matrix<XYZValue<float>>* const Project::emptyGridData = nullptr;
std::vector<XYZ<float>>* const Project::emptyPointList = nullptr;

Project::Project(USLab4a& mainClass)
		: method_(MethodEnum::invalid)
		, mainClass_(mainClass)
		, useGUI_(true)
{
}

void
Project::loadTaskParameters(const std::string& taskFileName)
{
	std::ostringstream taskFilePath;
	taskFilePath << directory_ << '/' << taskFileName;
	taskParameterMap_ = std::make_unique<const ParameterMap>(taskFilePath.str());
}

ParamMapPtr
Project::getParamMap(const char* fileName) const
{
	std::ostringstream filePath;
	filePath << expDirectory_ << '/' << fileName;
	return std::make_unique<const ParameterMap>(filePath.str());
}

ParamMapPtr
Project::getSubParamMap(const ParameterMap& pm, const char* fileNameKey) const
{
	const auto fileName = pm.value<std::string>(fileNameKey);
	std::ostringstream filePath;
	filePath << expDirectory_ << '/' << fileName;
	return std::make_unique<const ParameterMap>(filePath.str());
}

ParamMapPtr
Project::getSubParamMap(const char* fileNameKey) const
{
	const auto fileName = taskParamMap().value<std::string>(fileNameKey);
	std::ostringstream filePath;
	filePath << expDirectory_ << '/' << fileName;
	return std::make_unique<const ParameterMap>(filePath.str());
}

void
Project::handleShowFigure2DRequest()
{
	const std::lock_guard<std::mutex> locker(figure2DData_.mutex);
	if (!figure2DData_.showFigureRequested) return;
	mainClass_.showFigure2D(
			figure2DData_.figureId,
			figure2DData_.figureName,
			figure2DData_.xList,
			figure2DData_.yList,
			figure2DData_.markPoints);
	figure2DData_.showFigureRequested = false;
	figure2DData_.requestHandledCondition.notify_all();
}

void
Project::handleShowFigure3DRequest()
{
	const std::lock_guard<std::mutex> locker(figure3DData_.mutex);
	if (!figure3DData_.showFigureRequested) return;
	mainClass_.showFigure3D(
			figure3DData_.figureId,
			figure3DData_.figureName,
			figure3DData_.newGridData ? &figure3DData_.gridData : nullptr,
			figure3DData_.newPointList ? &figure3DData_.pointList : nullptr,
			figure3DData_.visualization,
			figure3DData_.colormap,
			figure3DData_.valueScale);
	figure3DData_.showFigureRequested = false;
	figure3DData_.requestHandledCondition.notify_all();
}

void
Project::handleShowMultiLayer3DRequest()
{
	const std::lock_guard<std::mutex> locker(multiLayer3DData_.mutex);
	if (!multiLayer3DData_.showFigureRequested) return;
	mainClass_.showMultiLayer3D(
			multiLayer3DData_.figureId,
			multiLayer3DData_.figureName,
			multiLayer3DData_.pointArray,
			multiLayer3DData_.indexArray);
	multiLayer3DData_.showFigureRequested = false;
	multiLayer3DData_.requestHandledCondition.notify_all();
}

void
Project::executeProgram(std::string& programPath, std::vector<std::string>& programArgs)
{
	QStringList args;
	for (std::vector<std::string>::const_iterator iter = programArgs.begin(); iter != programArgs.end(); ++iter) {
		args << QString(iter->c_str());
	}

	auto proc = std::make_unique<QProcess>();
	proc->setWorkingDirectory(expDirectory_.c_str());

	LOG_DEBUG << "Executing program: " << programPath << " with arguments: " << args.join(" ").toStdString();
	proc->start(programPath.c_str(), args);
	bool ok = proc->waitForFinished(-1 /* no time limit */);
	if (!ok) {
		THROW_EXCEPTION(ExternalProgramExecutionException, "Could not execute the program: " << programPath <<
				" with arguments: " << args.join(" ").toStdString());
	}
	LOG_DEBUG << "Program standard output:\n" << QString(proc->readAllStandardOutput()).toStdString();

	int exitCode = proc->exitCode();
	if (exitCode != 0) {
		THROW_EXCEPTION(ExternalProgramExecutionException, "An error ocurred while executing the program: " << programPath <<
				" with arguments: " << args.join(" ").toStdString() <<
				"\nExit code: " << exitCode <<
				"\nError message:\n" << QString(proc->readAllStandardError()).toStdString());
	}
}

void
Project::requestProcessingCancellation()
{
	const std::lock_guard<std::mutex> locker(control_.mutex);
	control_.processingCancellationRequested = true;
	control_.triggerCondition.notify_all();
}

bool
Project::processingCancellationRequested()
{
	const std::lock_guard<std::mutex> locker(control_.mutex);
	if (control_.processingCancellationRequested) {
		control_.processingCancellationRequested = false;
		return true;
	} else {
		return false;
	}
}

void
Project::resetTrigger()
{
	const std::lock_guard<std::mutex> locker(control_.mutex);
	control_.trigger = false;
	control_.triggerCount = 0;
}

void
Project::trigger()
{
	const std::lock_guard<std::mutex> locker(control_.mutex);
	if (control_.trigger) {
		THROW_EXCEPTION(InvalidStateException, "Duplicated trigger.");
	}
	control_.trigger = true;
	control_.triggerCondition.notify_all();
}

bool
Project::waitForTrigger(std::size_t* triggerCount)
{
	std::unique_lock<std::mutex> locker(control_.mutex);

	if (control_.processingCancellationRequested) {
		control_.processingCancellationRequested = false;
		control_.trigger = false;
		control_.triggerCount = 0;
		return false;
	}

	while (!control_.trigger) {
		LOG_INFO << ">>>>> Waiting for trigger >>>>>";
		control_.triggerCondition.wait(locker);

		if (control_.processingCancellationRequested) {
			control_.processingCancellationRequested = false;
			control_.trigger = false;
			control_.triggerCount = 0;
			return false;
		}
	}
	control_.trigger = false;

	if (triggerCount != nullptr) {
		*triggerCount = control_.triggerCount;
	}
	++control_.triggerCount;
	return true;
}

void
Project::createDirectory(const std::string &path, bool mustNotExist) const
{
	std::ostringstream fullPath;
	fullPath << expDirectory_ << '/' << path;
	fs::path dir(fullPath.str());

	if (mustNotExist && fs::exists(dir)) {
		THROW_EXCEPTION(InvalidDirectoryException, "The directory/file " << dir << " already exists.");
	}

	fs::create_directories(dir);
}

bool
Project::directoryExists(const std::string& path) const
{
	std::ostringstream fullPath;
	fullPath << expDirectory_ << '/' << path;
	fs::path dir(fullPath.str());
	return fs::exists(dir) && fs::is_directory(dir);
}

} // namespace Lab
