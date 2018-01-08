#ifndef FIGURE2DWINDOW_H_
#define FIGURE2DWINDOW_H_

#include <vector>

#include "ui_Figure2DWindow.h"



QT_FORWARD_DECLARE_CLASS(QMouseEvent)
class QCPAxis;

namespace Lab {

class Figure2DWindow : public QWidget {
	Q_OBJECT
public:
	Figure2DWindow(QWidget* parent=0);
	~Figure2DWindow();

	void updateData(const std::vector<double>& xList, const std::vector<double>& yList, bool markPoints);

private slots:
	void on_figure2DWidget_axisClick(QCPAxis* axis, QCPAxis::SelectablePart part, QMouseEvent* event);
	void on_figure2DWidget_selectionChangedByUser();
	void on_figure2DWidget_mouseDoubleClick(QMouseEvent* event);
private:
	Ui::Figure2DWindowClass ui_;
};

} // namespace Lab

#endif /* FIGURE2DWINDOW_H_ */
