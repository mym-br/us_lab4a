#include "Project.h"

#include <string>

#include <boost/make_shared.hpp>

#include <QDir>
#include <QFile>
#include <QMutexLocker>
#include <QProcess>
#include <QStringList>

#include "KeyValueFileReader.h"
#include "Log.h"
#include "USLab4a.h"



namespace Lab {

Matrix2<XZValue<float>>* Project::emptyGridData;
std::vector<XZ<float>>*  Project::emptyPointList;

Project::Project(USLab4a& mainWindow)
		: method_(MethodType::invalid)
		, mainWindow_(mainWindow)
{
}

Project::~Project()
{
}

void
Project::loadTaskParameters(const std::string& taskFileName)
{
	QFile taskFile(directory_ + '/' + taskFileName.c_str());
	if (taskFile.exists()) {
//		ui_.taskFileLineEdit->setText(filePath);
//		QFileInfo fileInfo(taskFile);

//		project_.setTaskFilePath(filePath);
		taskParameterMap_ = boost::make_shared<const ParameterMap>(taskFile.fileName().toStdString().c_str());
	} else {
		THROW_EXCEPTION(InvalidFileException, "The file \"" << taskFile.fileName().toStdString() << "\" does not exist.");
	}
}

ConstParameterMapPtr
Project::loadParameterMap(const char* fileName) const
{
	QString filePath = directory_ + '/' + fileName;
	return boost::make_shared<const ParameterMap>(filePath);
}

ConstParameterMapPtr
Project::loadChildParameterMap(ConstParameterMapPtr pm, const char* fileNameKey) const
{
	std::string fileName = pm->value<std::string>(fileNameKey);
	QString filePath = directory_ + '/' + fileName.c_str();
	return boost::make_shared<const ParameterMap>(filePath);
}

void
Project::handleShowFigure2DRequest()
{
	QMutexLocker locker(&figure2DData_.mutex);
	if (!figure2DData_.showFigureRequested) return;
	mainWindow_.showFigure2D(
			figure2DData_.figureId,
			figure2DData_.figureName,
			figure2DData_.xList,
			figure2DData_.yList,
			figure2DData_.markPoints);
	figure2DData_.showFigureRequested = false;
	figure2DData_.requestHandledCondition.wakeAll();
}

void
Project::handleShowFigure3DRequest()
{
	QMutexLocker locker(&figure3DData_.mutex);
	if (!figure3DData_.showFigureRequested) return;
	mainWindow_.showFigure3D(
			figure3DData_.figureId,
			figure3DData_.figureName,
			figure3DData_.newGridData ? &figure3DData_.gridData : 0,
			figure3DData_.newPointList ? &figure3DData_.pointList : 0,
			figure3DData_.visualization,
			figure3DData_.colormap,
			figure3DData_.valueScale);
	figure3DData_.showFigureRequested = false;
	figure3DData_.requestHandledCondition.wakeAll();
}

void
Project::executeProgram(std::string& programPath, std::vector<std::string>& programArgs)
{
	QStringList args;
	for (std::vector<std::string>::const_iterator iter = programArgs.begin(); iter != programArgs.end(); ++iter) {
		args << QString(iter->c_str());
	}

	boost::scoped_ptr<QProcess> proc(new QProcess);
	proc->setWorkingDirectory(directory_);

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
	QMutexLocker locker(&control_.mutex);
	control_.processingCancellationRequested = true;
	control_.triggerCondition.wakeAll();
}

bool
Project::processingCancellationRequested()
{
	QMutexLocker locker(&control_.mutex);
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
	QMutexLocker locker(&control_.mutex);
	control_.trigger = false;
	control_.triggerCount = 0;
}

void
Project::trigger()
{
	QMutexLocker locker(&control_.mutex);
	if (control_.trigger) {
		THROW_EXCEPTION(InvalidStateException, "Duplicated trigger.");
	}
	control_.trigger = true;
	control_.triggerCondition.wakeAll();
}

bool
Project::waitForTrigger(std::size_t* triggerCount)
{
	QMutexLocker locker(&control_.mutex);
	if (control_.processingCancellationRequested) {
		control_.processingCancellationRequested = false;
		control_.trigger = false;
		control_.triggerCount = 0;
		return false;
	}

	while (!control_.trigger) {
		LOG_INFO << ">>>>> Waiting for trigger >>>>>";
		control_.triggerCondition.wait(&control_.mutex);

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
Project::createDirectory(const std::string &path, bool mustNotExist)
{
	QDir dir;
	QString fullDir = directory_;
	fullDir.append('/');
	fullDir.append(path.c_str());
	if (mustNotExist && dir.exists(fullDir)) {
		THROW_EXCEPTION(InvalidDirectoryException, "The directory/file " << fullDir.toStdString() << " already exists.");
	}

	if (!dir.mkpath(fullDir)) {
		THROW_EXCEPTION(InvalidDirectoryException, "Could not create the directory " << fullDir.toStdString() << '.');
	}
}

} // namespace Lab
