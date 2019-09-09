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

#include <cstring>
#include <exception>
#include <memory>
#include <sstream>
#include <string>

#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QList>
#include <QListWidgetItem>
#include <QScrollBar>
#include <QSettings>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QTime>

#include "Controller.h"
#include "Exception.h"
#include "Figure2DWindow.h"
#include "Figure3DWindow.h"
#include "IterationCounter.h"
#include "Log.h"
#include "LogSyntaxHighlighter.h"
#include "Method.h"
#include "MultiLayer3DWindow.h"
#include "ParameterMap.h"
#include "Util.h"

#define FIGURE_WINDOWS_TIMER_PERIOD_MS 10
#define INFO_TIMER_PERIOD_MS 250
#define MAX_LOG_BLOCK_COUNT 500
#define FLOATING_POINT_VALUE_DECIMALS 6
#define TASK_FILE_PREFIX "task-"
#define TASK_FILE_SUFFIX ".txt"

namespace Lab {

USLab4a::USLab4a(QWidget* parent)
		: QMainWindow(parent)
		, controller_(std::make_unique<Controller>(project_))
		, project_(*this)
		, figureWindowsTimer_(this)
		, infoTimer_(this)
		, logFile_(LOG_FILE_NAME)
		, nextScriptEntry_()
{
	ui_.setupUi(this);

	connect(controller_.get(), SIGNAL(processingFinished()), this, SLOT(handleControllerFinishedProcessing()));

	ui_.logPlainTextEdit->setMaximumBlockCount(MAX_LOG_BLOCK_COUNT);

	figureWindowsTimer_.start(FIGURE_WINDOWS_TIMER_PERIOD_MS);
	connect(&figureWindowsTimer_, SIGNAL(timeout()), this, SLOT(updateFigureWindows()));

	infoTimer_.start(INFO_TIMER_PERIOD_MS);
	connect(&infoTimer_, SIGNAL(timeout()), this, SLOT(updateInfo()));

	new LogSyntaxHighlighter(ui_.logPlainTextEdit->document());

	QSettings settings;

	// Initial project directory.
	QVariant projectDir = settings.value(SETTINGS_KEY_PROJECT_DIR);
	if (!projectDir.isNull()) {
		QString path = projectDir.toString();
		ui_.projectDirLineEdit->setText(path);
		project_.setDirectory(path.toStdString());
	} else {
		QString currentPath = QDir::currentPath();
		ui_.projectDirLineEdit->setText(currentPath);
		project_.setDirectory(currentPath.toStdString());
	}
	try {
		fillTaskListWidget();

		QVariant selectedTask = settings.value(SETTINGS_KEY_SELECTED_TASK);
		if (!selectedTask.isNull()) {
			QList<QListWidgetItem*> taskItems = ui_.taskListWidget->findItems(selectedTask.toString(), Qt::MatchFixedString | Qt::MatchCaseSensitive);
			if (taskItems.size() == 1) {
				taskItems[0]->setSelected(true);
			}
		}
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::USLab4a] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::USLab4a] Caught an unknown exception.";
	}

	ui_.logLevelComboBox->addItem(tr("Error")  , Log::LEVEL_ERROR);
	ui_.logLevelComboBox->addItem(tr("Warning"), Log::LEVEL_WARNING);
	ui_.logLevelComboBox->addItem(tr("Info")   , Log::LEVEL_INFO);
	ui_.logLevelComboBox->addItem(tr("Debug")  , Log::LEVEL_DEBUG);
	// Initial log level.
	QVariant logLevelComboBoxIndex = settings.value(SETTINGS_KEY_LOGLEVEL_COMBOBOX_INDEX);
	if (!logLevelComboBoxIndex.isNull()) {
		ui_.logLevelComboBox->setCurrentIndex(logLevelComboBoxIndex.toInt());
		Log::setLevel(static_cast<Log::Level>(logLevelComboBoxIndex.toInt()));
	} else {
		ui_.logLevelComboBox->setCurrentIndex(0);
		Log::setLevel(static_cast<Log::Level>(0));
	}

	bool fileRetValue = logFile_.open(QIODevice::WriteOnly /*| QIODevice::Truncate*/);
	if (!fileRetValue) {
		THROW_EXCEPTION(IOException, "Could not open the log file.");
	}

	show();
}

USLab4a::~USLab4a()
{
}



void
USLab4a::fillTaskListWidget()
{
	//LOG_DEBUG << "USLab4a::fillTaskListWidget";

	ui_.taskListWidget->clear();

	QDir dir(project_.directory().c_str());
	if (dir.exists()) {
		QStringList nameFilters;
		nameFilters << TASK_FILE_PREFIX "*" TASK_FILE_SUFFIX;

		const std::size_t prefixLen = strlen(TASK_FILE_PREFIX);
		const std::size_t suffixLen = strlen(TASK_FILE_SUFFIX);
		const QStringList fileList = dir.entryList(nameFilters, QDir::Files | QDir::Readable, QDir::Name);
		if (!fileList.empty()) {
			for (QStringList::const_iterator iter = fileList.begin(); iter != fileList.end(); ++iter) {
				QString fileName = *iter;
				ui_.taskListWidget->addItem(fileName.mid(prefixLen, fileName.size() - prefixLen - suffixLen));
			}
		}
	} else {
		THROW_EXCEPTION(InvalidDirectoryException, "The directory " << dir.absolutePath().toStdString() << " does not exist.");
	}
}

void
USLab4a::processScriptEntry()
{
	if (nextScriptEntry_ < static_cast<unsigned int>(scriptEntryList_.size())) {

		QString projectDir = scriptEntryList_[nextScriptEntry_].project;
		QFileInfo projectInfo(projectDir);
		if (!projectInfo.exists() || !projectInfo.isDir()) {
			LOG_ERROR << "[USLab4a::processScriptEntry] Project path " << projectDir.toStdString() << " doesn't exist or is not a directory.";
			resetScriptData();
			return;
		}
		ui_.projectDirLineEdit->setText(projectDir);
		project_.setDirectory(projectDir.toStdString());

		QString taskFile = scriptEntryList_[nextScriptEntry_].task;
		LOG_INFO << "[SCRIPT] PROJECT " << projectDir.toStdString();
		LOG_INFO << "         TASK    " << taskFile.toStdString();

		try {
			project_.loadTaskParameters(taskFile.toStdString());

			ConstParameterMapPtr pm = project_.taskParameterMap();
			if (!pm) {
				THROW_EXCEPTION(InvalidStateException, "The task parameter map has not been initialized.");
			}

			std::string methodName = pm->value<std::string>("method");
			ui_.methodLineEdit->setText(methodName.c_str());

			MethodEnum method = Method::findByName(methodName);
			project_.setMethod(method);
		} catch (std::exception& e) {
			LOG_ERROR << "[USLab4a::processScriptEntry] Caught exception: " << e.what() << '.';
			resetScriptData();
			return;
		} catch (...) {
			LOG_ERROR << "[USLab4a::processScriptEntry] Caught an unknown exception.";
			resetScriptData();
			return;
		}

		closeAllFigures();

		ui_.projectDirLineEdit->setDisabled(true);
		ui_.selectProjectDirButton->setDisabled(true);
		ui_.scanProjectDirButton->setDisabled(true);
		ui_.taskListWidget->setDisabled(true);
		ui_.enableTaskButton->setDisabled(true);

		controller_->enableProcessing();

		++nextScriptEntry_;
	} else {
		resetScriptData();
		LOG_INFO << "[SCRIPT] Processing finished.";
	}
}

void
USLab4a::resetScriptData()
{
	scriptEntryList_.clear();
	nextScriptEntry_ = 0;
}

void
USLab4a::handleControllerFinishedProcessing()
{
	ui_.projectDirLineEdit->setDisabled(false);
	ui_.selectProjectDirButton->setDisabled(false);
	ui_.scanProjectDirButton->setDisabled(false);
	ui_.taskListWidget->setDisabled(false);
	ui_.enableTaskButton->setText(tr("Enable"));
	ui_.enableTaskButton->setDisabled(false);

	IterationCounter::reset(0);
	ui_.progressBar->setValue(0);
	ui_.progressBar->setEnabled(false);

	if (!scriptEntryList_.empty()) {
		processScriptEntry();
	}
}

void
USLab4a::on_openScriptAction_triggered()
{
	QString filePath = QFileDialog::getOpenFileName(this, tr("Select script file"), project_.directory().c_str(), tr("Text files (*.txt)"));
	if (!filePath.isEmpty()) {
		qDebug("on_openScriptAction_triggered filePath=%s", filePath.toStdString().c_str());

		QFileInfo fileInfo(filePath);
		QString scriptDir = fileInfo.canonicalPath();
		qDebug("on_openScriptAction_triggered scriptDir=%s", scriptDir.toStdString().c_str());

		QFile file(filePath);
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
			LOG_ERROR << "[USLab4a::on_openScriptAction_triggered] The file " << filePath.toStdString() << " could not be opened.";
			return;
		}

		QTextStream in(&file);
		in.setCodec("UTF-8");
		int lineNumber = 0;
		while (!in.atEnd()) {
			++lineNumber;

			QString line = in.readLine().trimmed();

			if (line.startsWith('#')) continue; // comment
			if (line.isEmpty()) continue;

			QStringList fieldList = line.split(' ', QString::SkipEmptyParts);
			if (fieldList.size() == 1) {
				LOG_ERROR << "[USLab4a::on_openScriptAction_triggered] Missing task at line " << lineNumber << " of file " << filePath.toStdString() << '.';
				resetScriptData();
				return;
			}
			if (fieldList.size() != 2) {
				LOG_ERROR << "[USLab4a::on_openScriptAction_triggered] Invalid syntax at line " << lineNumber << " of file " << filePath.toStdString() << '.';
				resetScriptData();
				return;
			}

			scriptEntryList_ << ScriptEntry{scriptDir + '/' + fieldList[0], fieldList[1]};
		}

		ui_.taskListWidget->clear();

		processScriptEntry();
	}
}

void
USLab4a::on_exitAction_triggered()
{
	qApp->closeAllWindows();
}

void
USLab4a::on_closeAllFiguresAction_triggered()
{
	closeAllFigures();
}

void
USLab4a::on_logLevelComboBox_activated(int index)
{
	Log::setLevel(static_cast<Log::Level>(ui_.logLevelComboBox->itemData(index).toInt()));
}

void
USLab4a::on_enableTaskButton_clicked()
{
	if (controller_->processingEnabled()) {
		ui_.enableTaskButton->setDisabled(true);
		project_.requestProcessingCancellation();
		return;
	}

	QList<QListWidgetItem*> selectedTasks = ui_.taskListWidget->selectedItems();
	if (selectedTasks.size() != 1) {
		LOG_ERROR << "A task must be selected.";
		return;
	}

	ui_.progressBar->setEnabled(true);
	ui_.progressBar->setValue(0);
	ui_.progressLabel->clear();
	IterationCounter::reset(0);

	try {
		project_.loadTaskParameters((QString(TASK_FILE_PREFIX) + selectedTasks.first()->text() + TASK_FILE_SUFFIX).toStdString());

		ConstParameterMapPtr pm = project_.taskParameterMap();
		if (!pm) {
			THROW_EXCEPTION(InvalidStateException, "The task parameter map has not been initialized.");
		}

		std::string methodName = pm->value<std::string>("method");
		ui_.methodLineEdit->setText(methodName.c_str());

		MethodEnum method = Method::findByName(methodName);
		project_.setMethod(method);
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::on_enableTaskButton_clicked] Caught exception: " << e.what() << '.';
		return;
	} catch (...) {
		LOG_ERROR << "[USLab4a::on_enableTaskButton_clicked] Caught an unknown exception.";
		return;
	}

	closeAllFigures();

	ui_.projectDirLineEdit->setDisabled(true);
	ui_.selectProjectDirButton->setDisabled(true);
	ui_.scanProjectDirButton->setDisabled(true);
	ui_.taskListWidget->setDisabled(true);
	ui_.enableTaskButton->setText(tr("Disable"));

	controller_->enableProcessing();
}

void
USLab4a::on_scanProjectDirButton_clicked()
{
	try {
		fillTaskListWidget();
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::on_scanProjectDirButton_clicked] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::on_scanProjectDirButton_clicked] Caught an unknown exception.";
	}
}

void
USLab4a::on_selectProjectDirButton_clicked()
{
	try {
		QString dirPath = QFileDialog::getExistingDirectory(this,
					tr("Select project directory"),
					project_.directory().c_str(),
					QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
						/*| QFileDialog::DontUseNativeDialog*/);
		if (!dirPath.isEmpty()) {
			ui_.projectDirLineEdit->setText(dirPath);
			project_.setDirectory(dirPath.toStdString());
			fillTaskListWidget();
		}
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::on_selectProjectDirButton_clicked] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::on_selectProjectDirButton_clicked] Caught an unknown exception.";
	}
}

void
USLab4a::on_triggerButton_clicked()
{
	try {
		project_.trigger();
		LOG_INFO << "##### TRIGGER #####";
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::on_triggerButton_clicked] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::on_triggerButton_clicked] Caught an unknown exception.";
	}
}

void
USLab4a::closeEvent(QCloseEvent* event)
{
	if (controller_->processingEnabled()) {
		project_.requestProcessingCancellation();
	}

	controller_->exit();

	QSettings settings;
	settings.setValue(SETTINGS_KEY_PROJECT_DIR            , project_.directory().c_str());
	settings.setValue(SETTINGS_KEY_LOGLEVEL_COMBOBOX_INDEX, ui_.logLevelComboBox->currentIndex());

	QList<QListWidgetItem*> selectedTasks = ui_.taskListWidget->selectedItems();
	if (selectedTasks.size() == 1) {
		settings.setValue(SETTINGS_KEY_SELECTED_TASK, selectedTasks.first()->text());
	} else {
		settings.remove(SETTINGS_KEY_SELECTED_TASK);
	}

	event->ignore();
	qApp->quit();
}

void
USLab4a::updateFigureWindows()
{
	project_.handleShowFigure2DRequest();
	project_.handleShowFigure3DRequest();
	project_.handleShowMultiLayer3DRequest();
}

void
USLab4a::updateInfo()
{
	try {
		// Update log.
		std::string s;
		Log::transferTo(s);
		if (!s.empty()) {
			ui_.logPlainTextEdit->appendPlainText(QString::fromStdString(s));

			//ui_.logPlainTextEdit->ensureCursorVisible();
			QScrollBar* scrollBar = ui_.logPlainTextEdit->verticalScrollBar();
			scrollBar->setValue(scrollBar->maximum());

			qint64 n = logFile_.write("\n");
			if (n == -1) {
				LOG_ERROR << "[USLab4a::updateInfo] An error ocurred while writing to the log file.";
			}
			n = logFile_.write(s.c_str());
			if (n == -1) {
				LOG_ERROR << "[USLab4a::updateInfo] An error ocurred while writing to the log file.";
			}
		}

		// Update progress bar.
		const unsigned int icTotal = IterationCounter::total;
		if (icTotal != 0) {
			const unsigned int icCount = IterationCounter::count;
			const double icTime = IterationCounter::timer.getTime();
			ui_.progressBar->setValue(static_cast<int>(((static_cast<double>(icCount) / icTotal) * 100.0 + 0.5)));
			QTime t0(0, 0);
			QTime t1 = t0.addSecs(static_cast<int>(icTime));
			ui_.progressLabel->setText(t1.toString() + ' ');
		}
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::updateInfo] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::updateInfo] Caught an unknown exception.";
	}
}

void
USLab4a::closeAllFigures()
{
	figure2DWindowList_.clear();
	figure3DWindowList_.clear();
	multiLayer3DWindowList_.clear();
}

void
USLab4a::showFigure2D(
		int id,
		const std::string& figureName,
		const std::vector<double>& xList,
		const std::vector<double>& yList,
		bool markPoints)
{
	try {
		Figure2DWindow& fig = figure2DWindowList_.get(id);
		std::ostringstream title;
		title << "2-" << id << ' ' << figureName;
		fig.setWindowTitle(title.str().c_str());
		fig.updateData(xList, yList, markPoints);
		fig.show();
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::showFigure2D] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::showFigure2D] Caught an unknown exception.";
	}
}

void
USLab4a::showFigure3D(
		int id,
		const std::string& figureName,
		Project::GridDataType* gridData,
		const std::vector<Project::PointType>* pointList,
		Figure::Visualization visualization,
		Figure::Colormap colormap,
		double valueScale)
{
	try {
		Figure3DWindow& fig = figure3DWindowList_.get(id);
		std::ostringstream title;
		title << "3-" << id << ' ' << figureName;
		fig.setWindowTitle(title.str().c_str());
		fig.setVisualization(visualization);
		fig.setColormap(colormap);
		fig.updateData(valueScale, gridData, pointList);
		fig.show();
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::showFigure3D] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::showFigure3D] Caught an unknown exception.";
	}
}

void
USLab4a::showMultiLayer3D(
		int id,
		const std::string& figureName,
		const std::vector<XYZValue<float>>& pointArray,
		const std::vector<unsigned int>& indexArray)
{
	try {
		MultiLayer3DWindow& fig = multiLayer3DWindowList_.get(id);
		std::ostringstream title;
		title << "M3-" << id << ' ' << figureName;
		fig.setWindowTitle(title.str().c_str());
		fig.updateData(pointArray, indexArray);
		fig.show();
	} catch (std::exception& e) {
		LOG_ERROR << "[USLab4a::showMultiLayer3D] Caught exception: " << e.what() << '.';
	} catch (...) {
		LOG_ERROR << "[USLab4a::showMultiLayer3D] Caught an unknown exception.";
	}
}

} // namespace Lab
