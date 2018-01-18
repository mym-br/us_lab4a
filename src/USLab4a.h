#ifndef USLAB4A_H
#define USLAB4A_H

#include <vector>

#include <boost/scoped_ptr.hpp>

#include <QFile>
#include <QMainWindow>
#include <QTimer>

#include "Controller.h"
#include "FigureWindowList.h"
#include "global.h"
#include "Project.h"
#include "ui_USLab4a.h"



QT_BEGIN_NAMESPACE
class QCloseEvent;
//class QErrorMessage;
class QString;
QT_END_NAMESPACE

namespace Lab {

class Figure2DWindow;
class Figure3DWindow;

struct ScriptEntry {
	QString project;
	QString task;
};

class USLab4a : public QMainWindow {
	Q_OBJECT
public:
	USLab4a(QWidget* parent = 0);
	virtual ~USLab4a();

	void closeAllFigures();
	void showFigure2D(
		int id,
		const std::string& figureName,
		const std::vector<double>& xList,
		const std::vector<double>& yList,
		bool markPoints);
	void showFigure3D(
		int id,
		const std::string& figureName,
		const Project::GridDataType* gridData = 0,
		const std::vector<Project::PointType>* pointList = 0,
		Figure::Visualization visualization = Figure::VISUALIZATION_DEFAULT,
		Figure::Colormap colormap = Figure::COLORMAP_DEFAULT,
		double valueScale = 0.0);
private slots:
	void handleControllerFinishedProcessing();
	void on_openScriptAction_triggered();
	void on_exitAction_triggered();
	void on_closeAllFiguresAction_triggered();
	void on_logLevelComboBox_activated(int index);
	void on_enableTaskButton_clicked();
	void on_scanProjectDirButton_clicked();
	void on_selectProjectDirButton_clicked();
	void on_triggerButton_clicked();
	void updateFigureWindows();
	void updateLogWidget();
private:
	USLab4a(const USLab4a&);
	USLab4a& operator=(const USLab4a&);

	virtual void closeEvent(QCloseEvent* event);
	void fillTaskListWidget();
	void startAcquisition();
	void stopAcquisition();
	void processScriptEntry();
	void resetScriptData();

	//QErrorMessage* errorMessage_;
	boost::scoped_ptr<Controller> controller_;
	FigureWindowList<Figure2DWindow> figure2DWindowList_;
	FigureWindowList<Figure3DWindow> figure3DWindowList_;
	Project project_;
	QTimer figureWindowsTimer_;
	QTimer logWidgetTimer_;
	Ui::USLab4aClass ui_;
	QFile logFile_;

	QList<ScriptEntry> scriptEntryList_;
	unsigned int nextScriptEntry_;
};

} // namespace Lab

#endif // USLAB4A_H
