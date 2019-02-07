QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = us_lab4a
TEMPLATE = app

CONFIG += c++14 warn_on

SOURCES += src/main.cpp\
    src/Project.cpp \
    src/ProcessingThread.cpp \
    src/ProcessingNode.cpp \
    src/OGLFigureWidget.cpp \
    src/Controller.cpp \
    src/fft/FFTW.cpp \
    src/method/TestMethod.cpp \
    src/network_acquisition/ArrayAcqClient.cpp \
    src/util/MinstdPseudorandomNumberGenerator.cpp \
    src/util/LogSyntaxHighlighter.cpp \
    src/util/Log.cpp \
    src/util/KeyValueFileReader.cpp \
    src/util/HDF5Util.cpp \
    src/util/ParameterMap.cpp \
    src/method/Method.cpp \
    src/method/ShowImageMethod.cpp \
    src/method/SingleAcquisitionMethod.cpp \
    src/Figure3DWindow.cpp \
    src/qcustomplot/qcustomplot.cpp \
    src/Figure2DWindow.cpp \
    src/USLab4a.cpp \
    src/network_acquisition/PhasedArrayAcqClient.cpp \
    src/util/Util.cpp \
    src/network_sync/SyncServer.cpp \
    src/util/FileUtil.cpp \
    src/MultiLayer3DWindow.cpp \
    src/OGLMultiLayerWidget.cpp \
    src/method/MultiLayerImageMethod.cpp \
    src/method/VTKFileMultiImageMethod.cpp

HEADERS  += \
    src/Project.h \
    src/ProcessingThread.h \
    src/ProcessingNode.h \
    src/OGLFigureWidget.h \
    src/Controller.h \
    src/acquisition/STAAcquisition.h \
    src/acquisition/SavedSTAAcquisition.h \
    src/acquisition/MultiplexedAcquisition.h \
    src/configuration/STAConfiguration.h \
    src/fft/FFTW.h \
    src/method/TestMethod.h \
    src/method/STAMethod.h \
    src/method/Method.h \
    src/network_acquisition/RawBuffer.h \
    src/network_acquisition/NetworkSTAAcquisition.h \
    src/network_acquisition/NetworkMultiplexedAcquisition.h \
    src/network_acquisition/ArrayAcqProtocol.h \
    src/network_acquisition/ArrayAcqClient.h \
    src/processor/VectorialSTAProcessor.h \
    src/processor/SimpleSTAProcessor.h \
    src/processor/DefaultSTAProcessor.h \
    src/simulator/SimulatedSTAAcquisition.h \
    src/util/XZValue.h \
    src/util/XY.h \
    src/util/Waveform.h \
    src/util/Timer.h \
    src/util/Statistics.h \
    src/util/MinstdPseudorandomNumberGenerator.h \
    src/util/LogSyntaxHighlighter.h \
    src/util/Log.h \
    src/util/LinearInterpolator.h \
    src/util/KeyValueFileReader.h \
    src/util/Interpolator4X.h \
    src/util/IndexValue.h \
    src/util/HilbertEnvelope.h \
    src/util/HDF5Util.h \
    src/util/FFTWFilter2.h \
    src/util/FFTWFilter.h \
    src/util/Exception.h \
    src/util/ContainerDumper.h \
    src/util/Array.h \
    src/util/ParameterMap.h \
    src/global.h \
    src/acquisition/SavedGroupAcquisition.h \
    src/method/ShowImageMethod.h \
    src/method/SingleAcquisitionMethod.h \
    src/Figure3DWindow.h \
    src/qcustomplot/qcustomplot.h \
    src/Figure2DWindow.h \
    src/FigureWindowList.h \
    src/util/Util.h \
    src/util/KaiserWindow.h \
    src/util/Interpolator.h \
    src/util/Value.h \
    src/util/RealToComplexFFT.h \
    src/util/ComplexToRealIFFT.h \
    src/parallel/ParallelHilbertEnvelope.h \
    src/util/ImageGrid.h \
    src/util/CoherenceFactor.h \
    src/util/ArrayGeometry.h \
    src/util/DirectFFTWFilter.h \
    src/fft/FFTUtil.h \
    src/util/MeasurementList.h \
    src/util/ExecutionTimeMeasurement.h \
    src/util/FileUtil.h \
    src/USLab4a.h \
    src/network_acquisition/PhasedArrayAcqClient.h \
    src/network_acquisition/PhasedArrayAcqProtocol.h \
    src/method/SimRectangularFlatSourceMethod.h \
    src/simulator/NumericRectangularFlatSourceImpulseResponse.h \
    src/simulator/SimTransientRadiationPattern.h \
    src/util/XYZ.h \
    src/simulator/SimTransientAcousticField.h \
    src/util/ArrayUtil.h \
    src/util/Decimator.h \
    src/simulator/Simulated3DAcquisitionDevice.h \
    src/simulator/Simulated3DSTAAcquisition.h \
    src/util/XYZValue.h \
    src/util/Geometry.h \
    src/simulator/AnalyticRectangularFlatSourceImpulseResponse.h \
    src/simulator/ArrayOfRectangularFlatSourcesImpulseResponse.h \
    src/configuration/DeviceSectorialScanConfiguration.h \
    src/network_acquisition/NetworkDeviceSectorialScanAcquisition.h \
    src/method/DeviceSectorialScanMethod.h \
    src/acquisition/DeviceSectorialScanAcquisition.h \
    src/network_sync/SyncServer.h \
    src/method/NetworkSyncSTAMethod.h \
    src/OGLMultiLayerWidget.h \
    src/MultiLayer3DWindow.h \
    src/util/OGL.h \
    src/method/MultiLayerImageMethod.h \
    src/method/VTKFileMultiImageMethod.h \
    src/util/XYZValueFactor.h \
    src/method/STA3DMethod.h \
    src/processor/Vectorial3DSTAProcessor.h \
    src/method/T1R1SAFT3DMethod.h \
    src/simulator/Simulated3DT1R1SAFTAcquisition.h \
    src/processor/Vectorial3DT1R1SAFTProcessor.h \
    src/configuration/SA3DConfiguration.h \
    src/util/WindowFunction.h \
    src/util/Matrix.h \
    src/util/Tensor3.h \
    src/simulator/SimTransientPropagation.h \
    src/util/XYZValueArray.h \
    src/method/SingleVirtualSourceMethod.h \
    src/acquisition/TnRnAcquisition.h \
    src/simulator/Simulated3DTnRnAcquisition.h \
    src/configuration/TnRnConfiguration.h \
    src/processor/Vectorial3DTnRnProcessor.h \
    src/processor/ArrayProcessor.h \
    src/network_acquisition/NetworkTnRnAcquisition.h \
    src/acquisition/SavedTnRnAcquisition.h

FORMS    += \
    ui/Figure3DWindow.ui \
    ui/Figure2DWindow.ui \
    ui/USLab4a.ui \
    ui/MultiLayer3DWindow.ui

INCLUDEPATH += src \
    src/acquisition \
    src/configuration \
    src/fft \
    src/method \
    src/processor \
    src/simulator \
    src/util \
    src/network_acquisition \
    src/network_sync \
    src/qcustomplot \
    src/parallel
LIBS += -ltbb \
    -ltbbmalloc \
    -lhdf5_cpp \
    -lfftw3 \
    -lfftw3f \
    -lboost_system \
    -lrt \
    -lGLU
exists(/usr/include/hdf5/serial) {
    # Debian 9.
    INCLUDEPATH += /usr/include/hdf5/serial
    LIBS += -lhdf5_serial
} else {
    LIBS += -lhdf5
}
DEPENDPATH += src \
    src/acquisition \
    src/configuration \
    src/fft \
    src/method \
    src/processor \
    src/simulator \
    src/util \
    src/network_acquisition \
    src/qcustomplot \
    src/parallel
#QMAKE_CXXFLAGS += -std=c++14
QMAKE_CXXFLAGS_DEBUG = -march=native -O0 -g
QMAKE_CXXFLAGS_RELEASE = -march=native -O3
#QMAKE_CXXFLAGS_RELEASE = -march=native -O3 -ftree-vectorize -ftree-vectorizer-verbose=2
#QMAKE_LFLAGS += -std=c++14

MOC_DIR = tmp
OBJECTS_DIR = tmp
UI_DIR = tmp
