#ifndef PROJECT_H_
#define PROJECT_H_

#include <algorithm> /* copy */
#include <sstream>
#include <string>
#include <utility> /* pair */
#include <vector>

#include <QMutex>
#include <QString>
#include <QWaitCondition>

#include "Exception.h"
#include "global.h"
#include "HDF5Util.h"
#include "Matrix2.h"
#include "Method.h"
#include "ParameterMap.h"
#include "Util.h"
#include "Value.h"
#include "XZ.h"
#include "XZValue.h"



namespace Lab {

class ParameterMap;
class USLab4a;

// Only the member functions marked as thread-safe can be called by more than
// one thread at a time.
class Project {
public:
	typedef Matrix2<XZValue<float> > GridDataType;
	typedef XZ<float> PointType;

	Project(USLab4a& mainWindow);
	~Project();

	const std::string directory() const { return directory_.toStdString(); }
	void setDirectory(const std::string& dirPath) { directory_ = dirPath.c_str(); }

//	const QString& taskFile() const { return taskFile_; }
//	void setTaskFile(const QString& taskFile) { taskFile_ = taskFile; }

	void loadTaskParameters(const std::string& taskFile);
	ConstParameterMapPtr taskParameterMap() const {
		if (!taskParameterMap_) THROW_EXCEPTION(InvalidStateException, "The task parameter map has not been loaded.");
		return taskParameterMap_;
	}

	Method::Type method() const { return method_; }
	void setMethod(Method::Type method) { method_ = method; }

	// The caller becomes the owner of the ParameterMap object.
	ConstParameterMapPtr loadParameterMap(const char* fileName) const;
	ConstParameterMapPtr loadChildParameterMap(ConstParameterMapPtr pm, const char* fileNameKey) const;

//	template<typename F> void saveData(const Matrix2<F>& m, const char* subDir, const char* fileName) const;

	//template<typename T> void loadHDF5(const std::string& fileName, const std::string& datasetName, T& container) const;
	//template<typename T> void saveHDF5(const T& container, const std::string& fileName, const std::string& datasetName) const;
	template<typename FloatType> void loadHDF5(const std::string& fileName, const std::string& datasetName, std::vector<FloatType>& container) const;
	//template<typename FloatType> void loadHDF5(const std::string& fileName, const std::string& datasetName, std::vector<std::pair<FloatType, FloatType> >& container) const;
	template<typename FloatType> void loadHDF5(const std::string& fileName, const std::string& datasetName, Matrix2<FloatType>& container) const;
	template<typename FloatType> void saveHDF5(const std::vector<FloatType>& container, const std::string& fileName, const std::string& datasetName) const;
	//template<typename FloatType> void saveHDF5(const std::vector<std::pair<FloatType, FloatType> >& container, const std::string& fileName, const std::string& datasetName) const;
	template<typename FloatType> void saveHDF5(const Matrix2<FloatType>& container, const std::string& fileName, const std::string& datasetName) const;

	// These functions may be called by only one thread.
	template<typename T, typename U> void saveHDF5(const std::vector<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp);
	template<typename T, typename U> void saveHDF5(const Matrix2<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp);

	// Called by the producer.
	// This function blocks if waitPending = true and there is a pending request.
	// If waitPending = false and there is a pending request, this function does nothing.
	// This function is thread-safe.
	template<typename FloatType> void showFigure2D(int id,
						const char* figureName,
						std::vector<FloatType>& xList,
						std::vector<FloatType>& yList,
						bool waitPending = true,
						bool markPoints = false);

	// Called by the producer.
	// This function blocks if waitPending = true and there is a pending request.
	// If waitPending = false and there is a pending request, this function does nothing.
	// This function is thread-safe.
	template<typename T, typename U> void showFigure3D(
						int id,
						const char* figureName,
						const T* gridData,
						const U* pointList,
						bool waitPending = true,
						Figure::Visualization visualization = Figure::VISUALIZATION_DEFAULT,
						Figure::Colormap colormap = Figure::COLORMAP_DEFAULT);

	// Called by the consumer.
	// These functions are thread-safe.
	void handleShowFigure2DRequest();
	void handleShowFigure3DRequest();

	void executeProgram(std::string& programPath, std::vector<std::string>& programArgs);

	// These functions are thread-safe.
	void requestProcessingCancellation();
	bool processingCancellationRequested();

	static Matrix2<XZValue<float> >* emptyGridData;
	static std::vector<XZ<float> >* emptyPointList;

private:
	struct Figure2DData {
		Figure2DData()
			: showFigureRequested(false)
			, figureId(0)
			, figureName("Figure")
			, markPoints(false)
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
			: showFigureRequested(false)
			, newGridData(false)
			, newPointList(false)
			, figureId(0)
			, visualization(Figure::VISUALIZATION_DEFAULT)
			, colormap(Figure::COLORMAP_DEFAULT)
			, figureName("Figure")
		{ }
		QMutex mutex;
		QWaitCondition requestHandledCondition;
		bool showFigureRequested;
		bool newGridData;
		bool newPointList;
		int figureId;
		Figure::Visualization visualization;
		Figure::Colormap colormap;
		std::string figureName;
		GridDataType gridData;
		std::vector<PointType> pointList;
	};

	struct Flags {
		Flags()
			: processingCancellationRequested(false)
		{ }
		QMutex mutex;
		bool processingCancellationRequested;
	};

	Project(const Project&);
	Project& operator=(const Project&);

	Method::Type method_;
	USLab4a& mainWindow_;
	ConstParameterMapPtr taskParameterMap_;
	QString directory_;
	//QString taskFile_;
	Figure2DData figure2DData_;
	Figure3DData figure3DData_;
	Flags flags_;

	// These containers are here to avoid reallocations.
	std::vector<double> auxHDF5Vector_;
	Matrix2<double> auxHDF5Matrix_;
};



template<typename FloatType>
void
Project::loadHDF5(const std::string& fileName, const std::string& datasetName, std::vector<FloatType>& container) const
{
	QString filePath = directory_ + '/' + QString::fromStdString(fileName) + HDF5_FILE_SUFFIX;
	HDF5Util::load2(filePath.toStdString(), datasetName, container);
}

//template<typename FloatType>
//void
//Project::loadHDF5(const std::string& fileName, const std::string& datasetName, std::vector<std::pair<FloatType, FloatType> >& container) const
//{
//	QString filePath = directory_ + '/' + QString::fromStdString(fileName) + HDF5_FILE_SUFFIX;
//	HDF5Util::load2(filePath.toStdString(), datasetName, container);
//}

template<typename FloatType>
void
Project::loadHDF5(const std::string& fileName, const std::string& datasetName, Matrix2<FloatType>& container) const
{
	QString filePath = directory_ + '/' + QString::fromStdString(fileName) + HDF5_FILE_SUFFIX;
	HDF5Util::load2(filePath.toStdString(), datasetName, container);
}

template<typename FloatType>
void
Project::saveHDF5(const std::vector<FloatType>& container, const std::string& fileName, const std::string& datasetName) const
{
	QString filePath = directory_ + '/' + QString::fromStdString(fileName) + HDF5_FILE_SUFFIX;
	HDF5Util::save2(container, filePath.toStdString(), datasetName);
}

//template<typename FloatType>
//void
//Project::saveHDF5(const std::vector<std::pair<FloatType, FloatType> >& container, const std::string& fileName, const std::string& datasetName) const
//{
//	QString filePath = directory_ + '/' + QString::fromStdString(fileName) + HDF5_FILE_SUFFIX;
//	HDF5Util::save2(container, filePath.toStdString(), datasetName);
//}

template<typename FloatType>
void
Project::saveHDF5(const Matrix2<FloatType>& container, const std::string& fileName, const std::string& datasetName) const
{
	QString filePath = directory_ + '/' + QString::fromStdString(fileName) + HDF5_FILE_SUFFIX;
	HDF5Util::save2(container, filePath.toStdString(), datasetName);
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
Project::saveHDF5(const Matrix2<T>& container, const std::string& fileName, const std::string& datasetName, U copyOp)
{
	auxHDF5Matrix_.resize(container.n1(), container.n2());
	Util::copyUsingOperator(container.begin(), container.end(), auxHDF5Matrix_.begin(), copyOp);
	saveHDF5(auxHDF5Matrix_, fileName, datasetName);
}

template<typename FloatType>
void
Project::showFigure2D(
		int id,
		const char* figureName,
		std::vector<FloatType>& xList,
		std::vector<FloatType>& yList,
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
		figure2DData_.figureName = figureName;

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
		Figure::Visualization visualization,
		Figure::Colormap colormap)
{
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
		figure3DData_.figureName = figureName;
		if (gridData) {
			figure3DData_.gridData.resize(gridData->n1(), gridData->n2());
			Value::copyXZValueSequence(gridData->begin(), gridData->end(), figure3DData_.gridData.begin());
			figure3DData_.newGridData = true;
		} else {
			figure3DData_.newGridData = false;
		}
		if (pointList) {
			figure3DData_.pointList.resize(pointList->size());
			Value::copyXZSequence(pointList->begin(), pointList->end(), figure3DData_.pointList.begin());
			figure3DData_.newPointList = true;
		} else {
			figure3DData_.newPointList = false;
		}

		figure3DData_.visualization = visualization;
		figure3DData_.colormap = colormap;
	}
}

} // namespace Lab

#endif /* PROJECT_H_ */
