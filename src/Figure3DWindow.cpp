#include "Figure3DWindow.h"

#define MAX_MIN_DECIBELS (-1)
#define MIN_MIN_DECIBELS (-100)
#define MIN_DECIBELS_STEP 1



namespace Lab {

Figure3DWindow::Figure3DWindow(QWidget *parent)
		: QWidget(parent)
{
	ui_.setupUi(this);

	visualization_ = ui_.oglFigureWidget->visualization();
	minDecibels_ = static_cast<int>(ui_.oglFigureWidget->minDecibels());
	colormap_ = ui_.oglFigureWidget->colormap();

//	scanWidgetTimer_.start(SCAN_WIDGET_TIMER_PERIOD_MS);
//	connect(&scanWidgetTimer_, SIGNAL(timeout()), this, SLOT(updateScanWidget()));

	ui_.visualizationComboBox->blockSignals(true); // to avoid calling on_*ComboBox_currentIndexChanged
	ui_.visualizationComboBox->addItem(tr("Raw - linear"));
	ui_.visualizationComboBox->addItem(tr("Rectified - linear"));
	ui_.visualizationComboBox->addItem(tr("Rectified - dB"));
	ui_.visualizationComboBox->addItem(tr("Envelope - linear"));
	ui_.visualizationComboBox->addItem(tr("Envelope - dB"));
	ui_.visualizationComboBox->blockSignals(false);
	ui_.visualizationComboBox->setCurrentIndex(static_cast<int>(visualization_));

	ui_.minDecibelsComboBox->blockSignals(true); // to avoid calling on_*ComboBox_currentIndexChanged
	for (int i = MAX_MIN_DECIBELS; i >= MIN_MIN_DECIBELS; i -= MIN_DECIBELS_STEP) {
		ui_.minDecibelsComboBox->addItem(QString::number(i));
	}
	ui_.minDecibelsComboBox->blockSignals(false);
	const int comboIdx = ui_.minDecibelsComboBox->findText(QString::number(minDecibels_));
	if (comboIdx != -1) {
		ui_.minDecibelsComboBox->setCurrentIndex(comboIdx);
	}

	ui_.colormapComboBox->blockSignals(true); // to avoid calling on_*ComboBox_currentIndexChanged
	ui_.colormapComboBox->addItem(tr("Gray"));
	ui_.colormapComboBox->addItem(tr("Inverted gray"));
	ui_.colormapComboBox->addItem(tr("Viridis"));
	ui_.colormapComboBox->addItem(tr("Inverted Viridis"));
	ui_.colormapComboBox->addItem(tr("Plasma"));
	ui_.colormapComboBox->addItem(tr("Inverted Plasma"));
	ui_.colormapComboBox->addItem(tr("Inferno"));
	ui_.colormapComboBox->addItem(tr("Inverted Inferno"));
	ui_.colormapComboBox->addItem(tr("Magma"));
	ui_.colormapComboBox->addItem(tr("Inverted Magma"));
	ui_.colormapComboBox->addItem(tr("Red-white-blue"));
	ui_.colormapComboBox->blockSignals(false);
	ui_.colormapComboBox->setCurrentIndex(static_cast<int>(colormap_));
}

Figure3DWindow::~Figure3DWindow()
{
}

void
Figure3DWindow::updateData(const Matrix2<XZValue<float> >* gridData, const std::vector<XZ<float> >* pointList)
{
	ui_.oglFigureWidget->resetScale();

//	bool newData = false;
	if (gridData != 0) {
		gridData_ = *gridData; // gets a copy
		ui_.oglFigureWidget->updateGridData(gridData_);
//		newData = true;
	}
	if (pointList != 0) {
		pointList_ = *pointList; // gets a copy
		ui_.oglFigureWidget->updatePointList(pointList_);
//		newData = true;
	}
//	if (newData) {
//		ui_.oglFigureWidget->update();
//	}
}

void
Figure3DWindow::on_visualizationComboBox_currentIndexChanged(int index)
{
	visualization_ = static_cast<Figure::Visualization>(index);
	ui_.oglFigureWidget->setVisualization(visualization_);
	ui_.oglFigureWidget->updateGridData(gridData_);
	//	ui_.oglFigureWidget->update();
}

void
Figure3DWindow::on_showInfoCheckBox_stateChanged(int state)
{
	ui_.oglFigureWidget->setShowInfo(state == Qt::Checked);
}

void
Figure3DWindow::on_rotationCheckBox_stateChanged(int state)
{
	ui_.oglFigureWidget->setRotationMode(state == Qt::Checked);
}

void
Figure3DWindow::on_minDecibelsComboBox_currentIndexChanged(const QString& text)
{
	minDecibels_ = text.toInt();
	ui_.oglFigureWidget->setMinDecibels(static_cast<float>(minDecibels_));

	if (visualization_ == Figure::VISUALIZATION_RECTIFIED_LOG ||
			visualization_ == Figure::VISUALIZATION_ENVELOPE_LOG) {
		ui_.oglFigureWidget->updateGridData(gridData_);
//		ui_.oglFigureWidget->update();
	}
}

void
Figure3DWindow::on_colormapComboBox_currentIndexChanged(int index)
{
	colormap_ = static_cast<Figure::Colormap>(index);
	ui_.oglFigureWidget->setColormap(colormap_);
	ui_.oglFigureWidget->updateGridData(gridData_);
//	ui_.oglFigureWidget->update();
}

void
Figure3DWindow::setVisualization(Figure::Visualization visualization)
{
	if (visualization != Figure::VISUALIZATION_DEFAULT) {
		ui_.oglFigureWidget->setVisualization(visualization);
		ui_.visualizationComboBox->setCurrentIndex(static_cast<int>(visualization));
	}
}

void
Figure3DWindow::setColormap(Figure::Colormap colormap)
{
	if (colormap != Figure::COLORMAP_DEFAULT) {
		ui_.oglFigureWidget->setColormap(colormap);
		ui_.colormapComboBox->setCurrentIndex(static_cast<int>(colormap));
	}
}

//void
//FigureWindow::on_valuesInDecibelCheckBox_clicked(bool checked)
//{
//	valuesInDecibel_ = checked;
//	ui_.oglFigureWidget->setDecibel(valuesInDecibel_, minDecibel_);
//	ui_.oglFigureWidget->updateGridData(gridData_);
//	ui_.oglFigureWidget->update();
//}

//void
//FigureWindow::updateScanWidget()
//{
//	if (scanWidgetNeedsUpdate_) {
//		scanWidgetNeedsUpdate_ = false;
//		if (project_.scanData().n1() == 0 || project_.scanData().n2() == 0) {
//			return;
//		}
//
//		if (project_.useEnvelope()) {
//			scanWidget_->updateData(project_.envelope());
//			oglScanWidget_->updateData(project_.envelope(), project_.reflectorList());
//oglScanWidget2_->updateData(project_.envelope(), project_.reflectorList());
//oglScanWidget3_->updateData(project_.envelope(), project_.reflectorList());
//
//
//figureWindowList_[0]->updateData(project_.envelope(), project_.reflectorList());
////figureWindowList_[1]->updateData(project_.envelope(), project_.reflectorList());
//
//
//		} else {
//			scanWidget_->updateData(project_.scanData());
//			oglScanWidget_->updateData(project_.scanData(), project_.reflectorList());
//oglScanWidget2_->updateData(project_.scanData(), project_.reflectorList());
//oglScanWidget3_->updateData(project_.scanData(), project_.reflectorList());
//
//figureWindowList_[0]->updateData(project_.envelope(), project_.reflectorList());
//		}
//
//		LOG_DEBUG << "oglScanWidget_->updateData()";
//	}
//}

} // namespace Lab
