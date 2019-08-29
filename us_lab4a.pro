QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = us_lab4a
TEMPLATE = app

CONFIG += c++14
!win32 {
    CONFIG += warn_on
}

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
    src/configuration/STAConfiguration.h \
    src/fft/FFTW.h \
    src/method/TestMethod.h \
    src/method/STAMethod.h \
    src/method/Method.h \
    src/network_acquisition/RawBuffer.h \
    src/network_acquisition/NetworkSTAAcquisition.h \
    src/network_acquisition/ArrayAcqProtocol.h \
    src/network_acquisition/ArrayAcqClient.h \
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
    src/util/XYZ.h \
    src/util/ArrayUtil.h \
    src/util/Decimator.h \
    src/util/XYZValue.h \
    src/util/Geometry.h \
    src/configuration/DeviceSectorialScanConfiguration.h \
    src/network_acquisition/NetworkDeviceSectorialScanAcquisition.h \
    src/method/DeviceSectorialScanMethod.h \
    src/network_sync/SyncServer.h \
    src/method/NetworkSyncSTAMethod.h \
    src/OGLMultiLayerWidget.h \
    src/MultiLayer3DWindow.h \
    src/util/OGL.h \
    src/method/MultiLayerImageMethod.h \
    src/method/VTKFileMultiImageMethod.h \
    src/util/XYZValueFactor.h \
    src/method/STA3DMethod.h \
    src/method/T1R1SAFT3DMethod.h \
    src/configuration/SA3DConfiguration.h \
    src/util/WindowFunction.h \
    src/util/Matrix.h \
    src/util/Tensor3.h \
    src/util/XYZValueArray.h \
    src/method/SingleVirtualSourceMethod.h \
    src/configuration/TnRnConfiguration.h \
    src/network_acquisition/NetworkTnRnAcquisition.h \
    src/method/NetworkSyncSingleVirtualSourceMethod.h \
    src/method/SyntheticYSingleVirtualSourceMethod.h \
    src/method/method_table.h \
    src/user/acquisition/DeviceSectorialScanAcquisition.h \
    src/user/acquisition/SavedSTAAcquisition.h \
    src/user/acquisition/SavedTnRnAcquisition.h \
    src/user/acquisition/STAAcquisition.h \
    src/user/acquisition/TnRnAcquisition.h \
    src/user/processor/ArrayProcessor.h \
    src/user/processor/DefaultSTAProcessor.h \
    src/user/processor/SimpleSTAProcessor.h \
    src/user/processor/SynthYVectorial3DTnRnProcessor.h \
    src/user/processor/Vectorial3DSTAProcessor.h \
    src/user/processor/Vectorial3DT1R1SAFTProcessor.h \
    src/user/processor/Vectorial3DTnRnProcessor.h \
    src/user/processor/VectorialSTAProcessor.h \
    src/user/simulator/AnalyticRectangularFlatSourceImpulseResponse.h \
    src/user/simulator/ArrayOfRectangularFlatSourcesImpulseResponse.h \
    src/user/simulator/NumericRectangularFlatSourceImpulseResponse.h \
    src/user/simulator/SimTransientAcousticField.h \
    src/user/simulator/SimTransientRadiationPattern.h \
    src/user/simulator/Simulated3DAcquisitionDevice.h \
    src/user/simulator/Simulated3DSTAAcquisition.h \
    src/user/simulator/Simulated3DT1R1SAFTAcquisition.h \
    src/user/simulator/Simulated3DTnRnAcquisition.h \
    src/user/simulator/SimulatedSTAAcquisition.h \
    src/user/simulator/SimTransientPropagation.h

FORMS    += \
    ui/Figure3DWindow.ui \
    ui/Figure2DWindow.ui \
    ui/USLab4a.ui \
    ui/MultiLayer3DWindow.ui

DEPENDPATH += src \
    src/configuration \
    src/fft \
    src/method \
    src/util \
    src/network_acquisition \
    src/network_sync \
    src/qcustomplot \
    src/parallel \
    src/user/acquisition \
    src/user/processor \
    src/user/simulator

INCLUDEPATH += src \
    src/configuration \
    src/fft \
    src/method \
    src/util \
    src/network_acquisition \
    src/network_sync \
    src/qcustomplot \
    src/parallel \
    src/user/acquisition \
    src/user/processor \
    src/user/simulator

win32 {
    # Windows 10 - VS 2017 - Qt 5.12.1
    INCLUDEPATH += \
        c:/lib/tbb2019_20181203oss/include \
        c:/lib/fftw-3.3.5-dll64 \
        c:/lib/boost_1_69_0 \
        "c:/Program Files/HDF_Group/HDF5/1.10.4/include"
    DEFINES += NOMINMAX H5_BUILT_AS_DYNAMIC_LIB
    QMAKE_CXXFLAGS += /W3
    # winmm: timeBeginPeriod, timeEndPeriod
    LIBS += \
        c:/lib/tbb2019_20181203oss/lib/intel64/vc14/tbb.lib \
        c:/lib/tbb2019_20181203oss/lib/intel64/vc14/tbbmalloc.lib \
        "c:/Program Files/HDF_Group/HDF5/1.10.4/lib/szip.lib" \
        "c:/Program Files/HDF_Group/HDF5/1.10.4/lib/zlib.lib" \
        "c:/Program Files/HDF_Group/HDF5/1.10.4/lib/hdf5.lib" \
        "c:/Program Files/HDF_Group/HDF5/1.10.4/lib/hdf5_cpp.lib" \
        c:/lib/fftw-3.3.5-dll64/libfftw3-3.lib \
        c:/lib/fftw-3.3.5-dll64/libfftw3f-3.lib \
        c:/lib/boost_1_69_0/lib64-msvc-14.1/libboost_system-vc141-mt-x64-1_69.lib \
        c:/lib/boost_1_69_0/lib64-msvc-14.1/libboost_date_time-vc141-mt-x64-1_69.lib \
        c:/lib/boost_1_69_0/lib64-msvc-14.1/libboost_regex-vc141-mt-x64-1_69.lib \
        opengl32.lib \
        glu32.lib \
        winmm.lib
    # Add to PATH for execution:
    # c:/lib/tbb2019_20181203oss/bin/intel64/vc14;c:/lib/fftw-3.3.5-dll64;c:/lib/boost_1_69_0/lib64-msvc-14.1
} else {
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

    #QMAKE_CXXFLAGS += -std=c++14
    QMAKE_CXXFLAGS_DEBUG = -march=native -O0 -g
    QMAKE_CXXFLAGS_RELEASE = -march=native -O3
    #QMAKE_CXXFLAGS_RELEASE = -march=native -O3 -ftree-vectorize -ftree-vectorizer-verbose=2
    #QMAKE_LFLAGS += -std=c++14
}

MOC_DIR = tmp
OBJECTS_DIR = tmp
UI_DIR = tmp
