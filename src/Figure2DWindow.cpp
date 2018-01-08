#include "Figure2DWindow.h"

#include <algorithm> /* copy */

#include <QList>
#include <QMouseEvent>
#include <QVector>

#include "qcustomplot.h"



namespace Lab {

Figure2DWindow::Figure2DWindow(QWidget* parent)
		: QWidget(parent)
{
	ui_.setupUi(this);

	ui_.figure2DWidget->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes);

	ui_.figure2DWidget->addGraph();
	ui_.figure2DWidget->graph(0)->setPen(QPen(Qt::black));
}

Figure2DWindow::~Figure2DWindow()
{
}

void
Figure2DWindow::updateData(const std::vector<double>& xList, const std::vector<double>& yList, bool markPoints)
{
	QVector<double> x = QVector<double>::fromStdVector(xList);
	QVector<double> y = QVector<double>::fromStdVector(yList);
	ui_.figure2DWidget->graph(0)->setData(x, y);
	ui_.figure2DWidget->graph(0)->rescaleAxes();
	if (markPoints) {
		ui_.figure2DWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssCross);
	}
	ui_.figure2DWidget->replot();
}

void
Figure2DWindow::on_figure2DWidget_axisClick(QCPAxis* axis, QCPAxis::SelectablePart /*part*/, QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton) {
		if (axis->axisType() == QCPAxis::atBottom) {
			ui_.figure2DWidget->axisRect()->setRangeZoom(Qt::Horizontal);
		} else if (axis->axisType() == QCPAxis::atLeft) {
			ui_.figure2DWidget->axisRect()->setRangeZoom(Qt::Vertical);
		}
	}
}

void
Figure2DWindow::on_figure2DWidget_selectionChangedByUser()
{
	QList<QCPAxis*> axes = ui_.figure2DWidget->selectedAxes();
	if (axes.empty()) {
		ui_.figure2DWidget->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
	}
}

void
Figure2DWindow::on_figure2DWidget_mouseDoubleClick(QMouseEvent* /*event*/)
{
	ui_.figure2DWidget->graph(0)->rescaleAxes();
	ui_.figure2DWidget->replot();
}

} // namespace Lab
