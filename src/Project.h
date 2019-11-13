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

#ifndef PROJECT_H_
#define PROJECT_H_

#include <algorithm> /* copy */
#include <cstddef> /* std::size_t */
#include <iomanip>
#include <sstream>
#include <string>
#include <utility> /* pair */
#include <vector>

#include <QMutex>
#include <QString>
#include <QWaitCondition>

#include "Colormap.h"
#include "Exception.h"
#include "FileUtil.h"
#include "global.h"
#include "HDF5Util.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "ParameterMap.h"
#include "Util.h"
#include "Value.h"
#include "Visualization.h"
#include "XYZ.h"
#include "XYZValue.h"



namespace Lab {

class ParameterMap;
class USLab4a;

// Only the member functions marked as thread-safe can be called by more than
// one thread at a time.
class Project {
public:
	typedef Matrix<XYZValue<float>> GridDataType;
	typedef XYZ<float> PointType;

	Project(USLab4a& mainWindow);
	~Project();

	const std::string directory() const { return directory_.toStdString(); }
	void setDirectory(const std::string& dirPath) { directory_ = dirPath.c_str(); }
	const std::string expDirectory() const { return expDirectory_.toStdString(); }
	void setExpDirectory(const std::string& dirPath) {
		if (dirPath.empty()) {
			expDirectory_ = directory_;
		} else {
			expDirectory_ = directory_ + '/' + dirPath.c_str();
		}
	}

	void loadTaskParameters(const std::string& taskFile);
	ParamMapPtr taskParameterMap() const {
		if (!taskParameterMap_) THROW_EXCEPTION(InvalidStateException, "The task parameter map has not been loaded.");
		return taskParameterMap_;
	}

	MethodEnum method() const { return method_; }
	void setMethod(MethodEnum method) { method_ = method; }

	// The caller becomes the owner of the ParameterMap object.
	ParamMapPtr loadParameterMap(const char* fileName) const;
	ParamMapPtr loadChildParameterMap(ParamMapPtr pm, const char* fileNameKey) const;

	template<typename FloatType> void loadHDF5(const std::string& fileName, const std::string& datasetName, std::vector<FloatType>& container) const;
	template<typename FloatType> void loadHDF5(const std::string& fileName, const std::string& datasetName, Matrix<FloatType>& container) const;
	template<typename FloatType> void saveHDF5(const std::vector<FloatType>& container, const std::string& fileName, const std::string& datasetName) const;
	template<typename FloatType> void saveHDF5(const Matrix<FloatType>& container, const std::string& fileName, const std::string& datasetName) const;
	template<typename T> void saveSignalToHDF5(const std::vector<T>& container, const std::string& outputDir,
					unsigned int acqNumber, unsigned int baseElement,
					unsigned int txElem, unsigned int rxElem);
	template<typename T> void saveSignalsToHDF5(const Matrix<T>& container, const std::string& outputDir,
					unsigned int acqNumber, unsigned int baseElement);
	template<typename T> void saveTxElemSignalsToHDF5(const Matrix<T>& container, const std::string& outputDir,
					unsigned int acqNumber, unsigned int baseElement, unsigned int txElem) const;
	template<typename T> void loadTxElemSignalsFromHDF5(const std::string& inputDir, unsigned int acqNumber,
					unsigned int baseElement, unsigned int txElem, Matrix<T>& container) const;

	// These functions may be called by only one thread.
	template<typename T, typename U> void loadHDF5(const std::string& fileName, const std::string& datasetName, Matrix<T>& container, U copyOp);
	template<typename T, typename U> void saveHDF5(const std::vector<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp);
	template<typename T, typename U> void saveHDF5(const Matrix<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp);
	template<typename T> void saveImageToHDF5(const Matrix<T>& container, const std::string& outputDir,
					const std::string& imageFile="image_value", const std::string& imageDataset="value");
	template<typename T> void saveFactorToHDF5(const Matrix<T>& container, const std::string& outputDir,
					const std::string& file, const std::string& dataset);
	template<typename T> void saveXYZToHDF5(const Matrix<T>& container, const std::string& outputDir,
					const std::string& xFile="image_x", const std::string& xDataset="x",
					const std::string& yFile="image_y", const std::string& yDataset="y",
					const std::string& zFile="image_z", const std::string& zDataset="z");
	template<typename T> void loadImageFromHDF5(const std::string& inputDir,
					const std::string& imageFile, const std::string& imageDataset,
					Matrix<T>& container);
	template<typename T> void loadXYZFromHDF5(const std::string& inputDir,
					const std::string& xFile, const std::string& xDataset,
					const std::string& yFile, const std::string& yDataset,
					const std::string& zFile, const std::string& zDataset,
					Matrix<T>& container);

	// Called by the producer.
	// This function blocks if waitPending = true and there is a pending request.
	// If waitPending = false and there is a pending request, this function does nothing.
	// This function is thread-safe.
	template<typename FloatType> void showFigure2D(int id,
						const char* figureName,
						const std::vector<FloatType>& xList,
						const std::vector<FloatType>& yList,
						bool waitPending=true,
						bool markPoints=false);

	// Called by the producer.
	// This function blocks if waitPending = true and there is a pending request.
	// If waitPending = false and there is a pending request, this function does nothing.
	// This function is thread-safe.
	template<typename T, typename U> void showFigure3D(
						int id,
						const char* figureName,
						const T* gridData,
						const U* pointList,
						bool waitPending=true,
						Visualization::Value visualization=Visualization::VALUE_DEFAULT,
						Colormap::Gradient colormap=Colormap::GRADIENT_DEFAULT,
						double valueScale=0.0);

	// Called by the producer.
	// This function blocks if there is a pending request.
	// If waitPending = false and there is a pending request, this function does nothing.
	// This function is thread-safe.
	template<typename T, typename U> void showMultiLayer3D(
						int id,
						const char* figureName,
						const T& pointArray,
						const U& indexArray);

	// Called by the consumer.
	// These functions are thread-safe.
	void handleShowFigure2DRequest();
	void handleShowFigure3DRequest();
	void handleShowMultiLayer3DRequest();

	void executeProgram(std::string& programPath, std::vector<std::string>& programArgs);

	// These functions are thread-safe.
	void requestProcessingCancellation();
	bool processingCancellationRequested();
	void resetTrigger();
	void trigger();
	bool waitForTrigger(std::size_t* triggerCount=nullptr); // returns false if processing cancellation has been requested

	static Matrix<XYZValue<float>>* const emptyGridData;
	static std::vector<XYZ<float>>* const emptyPointList;

	void createDirectory(const std::string& path, bool mustNotExist) const;
	bool directoryExists(const std::string& path) const;
private:
	struct Figure2DData {
		Figure2DData()
			: showFigureRequested()
			, figureId()
			, figureName("Figure")
			, markPoints()
		{ }
		QMutex mutex;
		QWaitCondition requestHandledCondition;
		bool showFigureRequested;
		int figureId;
		std::string figureName;
		std::vector<double> xList;
		std::vector<double> yList;
		bool markPoints;
	};

	struct Figure3DData {
		Figure3DData()
			: showFigureRequested()
			, newGridData()
			, newPointList()
			, figureId()
			, visualization(Visualization::VALUE_DEFAULT)
			, colormap(Colormap::GRADIENT_DEFAULT)
			, figureName("Figure")
			, valueScale()
		{ }
		QMutex mutex;
		QWaitCondition requestHandledCondition;
		bool showFigureRequested;
		bool newGridData;
		bool newPointList;
		int figureId;
		Visualization::Value visualization;
		Colormap::Gradient colormap;
		std::string figureName;
		GridDataType gridData;
		std::vector<PointType> pointList;
		double valueScale;
	};

	struct MultiLayer3DData {
		MultiLayer3DData()
			: showFigureRequested{}
			, figureId{}
			, figureName{"Figure"}
		{ }
		QMutex mutex;
		QWaitCondition requestHandledCondition;
		bool showFigureRequested;
		int figureId;
		std::string figureName;
		std::vector<XYZValue<float>> pointArray;
		std::vector<unsigned int> indexArray;
	};

	struct Control {
		Control()
			: triggerCount()
			, processingCancellationRequested()
			, trigger()
		{ }
		QMutex mutex;
		QWaitCondition triggerCondition;
		std::size_t triggerCount;
		bool processingCancellationRequested;
		bool trigger;
	};

	Project(const Project&) = delete;
	Project& operator=(const Project&) = delete;

	MethodEnum method_;
	USLab4a& mainWindow_;
	ParamMapPtr taskParameterMap_;
	QString directory_;
	QString expDirectory_;
	Figure2DData figure2DData_;
	Figure3DData figure3DData_;
	MultiLayer3DData multiLayer3DData_;
	Control control_;

	// These containers are here to avoid reallocations.
	std::vector<double> auxHDF5Vector_;
	Matrix<double> auxHDF5Matrix_;
	Matrix<double> aux2HDF5Matrix_;
	Matrix<double> aux3HDF5Matrix_;
};



template<typename FloatType>
void
Project::loadHDF5(const std::string& fileName, const std::string& datasetName, std::vector<FloatType>& container) const
{
	QString filePath = expDirectory_ + '/' + QString::fromStdString(fileName) + HDF5Util::fileSuffix;
	HDF5Util::load2(filePath.toStdString(), datasetName, container);
}

template<typename FloatType>
void
Project::loadHDF5(const std::string& fileName, const std::string& datasetName, Matrix<FloatType>& container) const
{
	QString filePath = expDirectory_ + '/' + QString::fromStdString(fileName) + HDF5Util::fileSuffix;
	HDF5Util::load2(filePath.toStdString(), datasetName, container);
}

template<typename FloatType>
void
Project::saveHDF5(const std::vector<FloatType>& container, const std::string& fileName, const std::string& datasetName) const
{
	QString filePath = expDirectory_ + '/' + QString::fromStdString(fileName) + HDF5Util::fileSuffix;
	HDF5Util::save2(container, filePath.toStdString(), datasetName);
}

template<typename FloatType>
void
Project::saveHDF5(const Matrix<FloatType>& container, const std::string& fileName, const std::string& datasetName) const
{
	QString filePath = expDirectory_ + '/' + QString::fromStdString(fileName) + HDF5Util::fileSuffix;
	HDF5Util::save2(container, filePath.toStdString(), datasetName);
}

template<typename T>
void
Project::saveSignalToHDF5(const std::vector<T>& container, const std::string& outputDir,
				unsigned int acqNumber, unsigned int baseElement,
				unsigned int txElem, unsigned int rxElem)
{
	std::string dirPath = FileUtil::path(outputDir, "/", acqNumber);
	createDirectory(dirPath, false);

	std::string filePath = FileUtil::signalPath(dirPath, baseElement, txElem, rxElem);

	LOG_DEBUG << "Saving " << filePath << "...";
	saveHDF5(container, filePath, "signal");
}

template<typename T>
void
Project::saveSignalsToHDF5(const Matrix<T>& container, const std::string& outputDir,
				unsigned int acqNumber, unsigned int baseElement)
{
	std::string dirPath = FileUtil::path(outputDir, "/", acqNumber);
	createDirectory(dirPath, false);

	std::string filePath = FileUtil::signalsPath(dirPath, baseElement);

	LOG_DEBUG << "Saving " << filePath << "...";
	saveHDF5(container, filePath, "signal");
}

template<typename T>
void
Project::saveTxElemSignalsToHDF5(const Matrix<T>& container, const std::string& outputDir,
					unsigned int acqNumber, unsigned int baseElement, unsigned int txElem) const
{
	std::string dirPath = FileUtil::path(outputDir, "/", acqNumber);
	createDirectory(dirPath, false);

	std::string filePath = FileUtil::txElemSignalsPath(dirPath, baseElement, txElem);

	LOG_DEBUG << "Saving " << filePath << "...";
	saveHDF5(container, filePath, "signal");
}

template<typename T>
void
Project::loadTxElemSignalsFromHDF5(const std::string& inputDir, unsigned int acqNumber, unsigned int baseElement,
				unsigned int txElem, Matrix<T>& container) const
{
	std::string dirPath = FileUtil::path(inputDir, "/", acqNumber);
	createDirectory(dirPath, false);

	std::string filePath = FileUtil::txElemSignalsPath(dirPath, baseElement, txElem);

	LOG_DEBUG << "Loading " << filePath << "...";
	loadHDF5(filePath, "signal", container);
}

template<typename T, typename U>
void
Project::loadHDF5(const std::string& fileName, const std::string& datasetName, Matrix<T>& container, U copyOp)
{
	loadHDF5(fileName, datasetName, auxHDF5Matrix_);
	container.resize(auxHDF5Matrix_.n1(), auxHDF5Matrix_.n2());
	Util::copyUsingOperator(auxHDF5Matrix_.begin(), auxHDF5Matrix_.end(), container.begin(), copyOp);
}

template<typename T, typename U>
void
Project::saveHDF5(const std::vector<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp)
{
	auxHDF5Vector_.resize(container.size());
	Util::copyUsingOperator(container.begin(), container.end(), auxHDF5Vector_.begin(), copyOp);
	saveHDF5(auxHDF5Vector_, fileName, datasetName);
}

template<typename T, typename U>
void
Project::saveHDF5(const Matrix<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp)
{
	auxHDF5Matrix_.resize(container.n1(), container.n2());
	Util::copyUsingOperator(container.begin(), container.end(), auxHDF5Matrix_.begin(), copyOp);
	saveHDF5(auxHDF5Matrix_, fileName, datasetName);
}

template<typename T>
void
Project::saveImageToHDF5(const Matrix<T>& container, const std::string& outputDir,
				const std::string& imageFile, const std::string& imageDataset)
{
	Util::copyValueToSimpleMatrix(container, auxHDF5Matrix_);
	LOG_DEBUG << "Saving " << imageFile << "...";
	saveHDF5(auxHDF5Matrix_, outputDir + '/' + imageFile, imageDataset);
}

template<typename T>
void
Project::saveFactorToHDF5(const Matrix<T>& container, const std::string& outputDir,
				const std::string& file, const std::string& dataset)
{
	Util::copyFactorToSimpleMatrix(container, auxHDF5Matrix_);
	LOG_DEBUG << "Saving " << file << "...";
	saveHDF5(auxHDF5Matrix_, outputDir + '/' + file, dataset);
}

template<typename T>
void
Project::saveXYZToHDF5(const Matrix<T>& container, const std::string& outputDir,
				const std::string& xFile, const std::string& xDataset,
				const std::string& yFile, const std::string& yDataset,
				const std::string& zFile, const std::string& zDataset)
{
	Util::copyXYZToSimpleMatrices(container, auxHDF5Matrix_, aux2HDF5Matrix_, aux3HDF5Matrix_);
	LOG_DEBUG << "Saving " << xFile << "...";
	saveHDF5(auxHDF5Matrix_ , outputDir + '/' + xFile, xDataset);
	LOG_DEBUG << "Saving " << yFile << "...";
	saveHDF5(aux2HDF5Matrix_, outputDir + '/' + yFile, yDataset);
	LOG_DEBUG << "Saving " << zFile << "...";
	saveHDF5(aux3HDF5Matrix_, outputDir + '/' + zFile, zDataset);
}

template<typename T>
void
Project::loadImageFromHDF5(const std::string& inputDir,
				const std::string& imageFile, const std::string& imageDataset,
				Matrix<T>& container)
{
	LOG_DEBUG << "Loading " << imageFile << "...";
	loadHDF5(inputDir + '/' + imageFile, imageDataset, auxHDF5Matrix_);
	container.resize(auxHDF5Matrix_.n1(), auxHDF5Matrix_.n2());
	Util::copyValueFromSimpleMatrix(auxHDF5Matrix_, container);
}

template<typename T>
void
Project::loadXYZFromHDF5(const std::string& inputDir,
				const std::string& xFile, const std::string& xDataset,
				const std::string& yFile, const std::string& yDataset,
				const std::string& zFile, const std::string& zDataset,
				Matrix<T>& container)
{
	LOG_DEBUG << "Loading " << xFile << "...";
	loadHDF5(inputDir + '/' + xFile, xDataset, auxHDF5Matrix_);
	LOG_DEBUG << "Loading " << yFile << "...";
	loadHDF5(inputDir + '/' + yFile, yDataset, aux2HDF5Matrix_);
	LOG_DEBUG << "Loading " << zFile << "...";
	loadHDF5(inputDir + '/' + zFile, zDataset, aux3HDF5Matrix_);
	Util::copyXYZFromSimpleMatrices(auxHDF5Matrix_, aux2HDF5Matrix_, aux3HDF5Matrix_, container);
}

template<typename FloatType>
void
Project::showFigure2D(int id,
		const char* figureName,
		const std::vector<FloatType>& xList,
		const std::vector<FloatType>& yList,
		bool waitPending,
		bool markPoints)
{
	if (xList.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "xList is empty.");
	}
	if (yList.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "yList is empty.");
	}
	if (xList.size() != yList.size()) {
		THROW_EXCEPTION(InvalidParameterException, "xList.size() != yList.size().");
	}
	{
		QMutexLocker locker(&figure2DData_.mutex);

		if (figure2DData_.showFigureRequested) {
			if (waitPending) {
				do {
					figure2DData_.requestHandledCondition.wait(&figure2DData_.mutex);
				} while (figure2DData_.showFigureRequested);
			} else {
				return;
			}
		}

		figure2DData_.showFigureRequested = true;
		figure2DData_.figureId = id;
		figure2DData_.figureName = figureName ? figureName : "Figure";

		figure2DData_.xList.resize(xList.size());
		std::copy(xList.begin(), xList.end(), figure2DData_.xList.begin());

		figure2DData_.yList.resize(yList.size());
		std::copy(yList.begin(), yList.end(), figure2DData_.yList.begin());

		figure2DData_.markPoints = markPoints;
	}
}

template<typename T, typename U>
void
Project::showFigure3D(
		int id,
		const char* figureName,
		const T* gridData,
		const U* pointList,
		bool waitPending,
		Visualization::Value visualization,
		Colormap::Gradient colormap,
		double valueScale)
{
	QMutexLocker locker(&figure3DData_.mutex);

	if (figure3DData_.showFigureRequested) {
		if (waitPending) {
			do {
				figure3DData_.requestHandledCondition.wait(&figure3DData_.mutex);
			} while (figure3DData_.showFigureRequested);
		} else {
			return;
		}
	}

	figure3DData_.showFigureRequested = true;
	figure3DData_.figureId = id;
	figure3DData_.figureName = figureName ? figureName : "Figure";
	if (gridData) {
		figure3DData_.gridData.resize(gridData->n1(), gridData->n2());
		Value::copyXYZValueSequence(gridData->begin(), gridData->end(), figure3DData_.gridData.begin());
		figure3DData_.newGridData = true;
	} else {
		figure3DData_.newGridData = false;
	}
	if (pointList) {
		figure3DData_.pointList.resize(pointList->size());
		Value::copyXYZSequence(pointList->begin(), pointList->end(), figure3DData_.pointList.begin());
		figure3DData_.newPointList = true;
	} else {
		figure3DData_.newPointList = false;
	}

	figure3DData_.visualization = visualization;
	figure3DData_.colormap = colormap;
	figure3DData_.valueScale = valueScale;
}

template<typename T, typename U>
void
Project::showMultiLayer3D(
		int id,
		const char* figureName,
		const T& pointArray,
		const U& indexArray)
{
	QMutexLocker locker(&multiLayer3DData_.mutex);

	if (multiLayer3DData_.showFigureRequested) {
		do {
			multiLayer3DData_.requestHandledCondition.wait(&multiLayer3DData_.mutex);
		} while (multiLayer3DData_.showFigureRequested);
	}

	multiLayer3DData_.showFigureRequested = true;
	multiLayer3DData_.figureId = id;
	multiLayer3DData_.figureName = figureName ? figureName : "Figure";

	multiLayer3DData_.pointArray.resize(pointArray.size());
	Value::copyXYZValueSequence(pointArray.begin(), pointArray.end(), multiLayer3DData_.pointArray.begin());
	multiLayer3DData_.indexArray = indexArray;
}

} // namespace Lab

#endif /* PROJECT_H_ */
