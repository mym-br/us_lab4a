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

#ifndef USLAB4A_H
#define USLAB4A_H

#include <memory>
#include <vector>

#include <QFile>
#include <QMainWindow>
#include <QTimer>

#include "Colormap.h"
#include "Controller.h"
#include "FigureWindowList.h"
#include "Project.h"
#include "ui_USLab4a.h"
#include "Visualization.h"



QT_BEGIN_NAMESPACE
class QCloseEvent;
class QString;
QT_END_NAMESPACE

namespace Lab {

class Figure2DWindow;
class Figure3DWindow;
class MultiLayer3DWindow;

struct BatchEntry {
	QString projectDir;
	QString task;
	QString experiment;
};

class USLab4a : public QMainWindow {
	Q_OBJECT
public:
	USLab4a(QWidget* parent=nullptr);
	virtual ~USLab4a();

	void closeAllFigures();
	void showFigure2D(
		int id,
		const std::string& figureName,
		const std::vector<double>& xList,
		const std::vector<double>& yList,
		bool markPoints);
	// *gridData will contain old data after the call.
	void showFigure3D(
		int id,
		const std::string& figureName,
		Project::GridDataType* gridData=nullptr,
		const std::vector<Project::PointType>* pointList=nullptr,
		Visualization::Value visualization=Visualization::Value::DEFAULT,
		Colormap::Id colormap=Colormap::Id::DEFAULT,
		double valueScale=0.0);
	void showMultiLayer3D(
		int id,
		const std::string& figureName,
		const std::vector<XYZValue<float>>& pointArray,
		const std::vector<unsigned int>& indexArray);
private Q_SLOTS:
	void handleControllerFinishedProcessing();
	void on_openBatchFileAction_triggered();
	void on_exitAction_triggered();
	void on_closeAllFiguresAction_triggered();
	void on_aboutAction_triggered();
	void on_logLevelComboBox_activated(int index);
	void on_enableTaskButton_clicked();
	void on_scanProjectDirButton_clicked();
	void on_selectProjectDirButton_clicked();
	void on_triggerButton_clicked();
	void updateFigureWindows();
	void updateInfo();
private:
	USLab4a(const USLab4a&) = delete;
	USLab4a& operator=(const USLab4a&) = delete;

	virtual void closeEvent(QCloseEvent* event);
	void fillTaskAndExpListWidget();
	void processBatchEntry();
	void resetBatchData();

	std::unique_ptr<Controller> controller_;
	FigureWindowList<Figure2DWindow> figure2DWindowList_;
	FigureWindowList<Figure3DWindow> figure3DWindowList_;
	FigureWindowList<MultiLayer3DWindow> multiLayer3DWindowList_;
	Project project_;
	QTimer figureWindowsTimer_;
	QTimer infoTimer_;
	Ui::USLab4aClass ui_;
	QFile logFile_;
	QList<BatchEntry> batchEntryList_;
	unsigned int nextBatchEntry_;
};

} // namespace Lab

#endif // USLAB4A_H
