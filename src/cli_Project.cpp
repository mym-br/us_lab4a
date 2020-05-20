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

#include "Project.h"

#include <filesystem>
#include <fstream>
#include <string>

#include "KeyValueFileReader.h"

#include "cli_USLab4a.h"

namespace fs = std::filesystem;

namespace Lab {

Matrix<XYZValue<float>>* const Project::emptyGridData = nullptr;
std::vector<XYZ<float>>* const Project::emptyPointList = nullptr;

Project::Project(USLab4a& mainClass)
		: method_(MethodEnum::invalid)
		, mainClass_(mainClass)
		, useGUI_()
{
	mainClass_.showMessage("[cli_Project.cpp] Project");
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
}

void
Project::handleShowFigure3DRequest()
{
}

void
Project::handleShowMultiLayer3DRequest()
{
}

void
Project::executeProgram(std::string& /*programPath*/, std::vector<std::string>& /*programArgs*/)
{
	THROW_EXCEPTION(InvalidCallException, "External program execution is disabled.");
}

void
Project::requestProcessingCancellation()
{
}

bool
Project::processingCancellationRequested()
{
	return true;
}

void
Project::resetTrigger()
{
}

void
Project::trigger()
{
}

bool
Project::waitForTrigger(std::size_t* /*triggerCount*/)
{
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
