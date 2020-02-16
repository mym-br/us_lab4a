cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)

project(us_lab4a)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOUIC_SEARCH_PATHS ui)
set(CMAKE_AUTORCC ON)

find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
pkg_check_modules(GLU REQUIRED IMPORTED_TARGET glu)
pkg_check_modules(HDF5 REQUIRED IMPORTED_TARGET hdf5)

if(UNIX)
    if(NOT APPLE)
        set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native")
    endif()
endif()

set(us_lab4a_SRC
    src/Controller.cpp
    src/Controller.h
    src/Figure2DWidget.cpp
    src/Figure2DWidget.h
    src/Figure2DWindow.cpp
    src/Figure2DWindow.h
    src/Figure3DWindow.cpp
    src/Figure3DWindow.h
    src/FigureWindowList.h
    src/main.cpp
    src/Method.cpp
    src/Method.h
    src/MultiLayer3DWindow.cpp
    src/MultiLayer3DWindow.h
    src/OGLFigureWidget.cpp
    src/OGLFigureWidget.h
    src/OGLMultiLayerWidget.cpp
    src/OGLMultiLayerWidget.h
    src/ProcessingNode.cpp
    src/ProcessingNode.h
    src/ProcessingThread.cpp
    src/ProcessingThread.h
    src/Project.cpp
    src/Project.h
    src/USLab4a.cpp
    src/USLab4a.h
    src/fft/FFTUtil.h
    src/fft/FFTW.cpp
    src/fft/FFTW.h
    src/network_acquisition/ArrayAcqClient.cpp
    src/network_acquisition/ArrayAcqClient.h
    src/network_acquisition/ArrayAcqProtocol.h
    src/network_acquisition/NetworkAcquisition.h
    src/network_acquisition/PhasedArrayAcqClient.cpp
    src/network_acquisition/PhasedArrayAcqClient.h
    src/network_acquisition/PhasedArrayAcqProtocol.h
    src/network_acquisition/RawBuffer.h
    src/network_sync/SyncServer.cpp
    src/network_sync/SyncServer.h
    src/parallel/ParallelHilbertEnvelope.h
    src/util/ArrayGeometry.h
    src/util/ArrayUtil.h
    src/util/CoherenceFactor.h
    src/util/Colormap.cpp
    src/util/Colormap.h
    src/util/ComplexToRealIFFT.h
    src/util/ContainerDumper.h
    src/util/Decimator.h
    src/util/DirectFFTWFilter.h
    src/util/Exception.h
    src/util/ExecutionTimeMeasurement.h
    src/util/FermatPrinciple.h
    src/util/FFTWFilter.h
    src/util/FFTWFilter2.h
    src/util/FileUtil.cpp
    src/util/FileUtil.h
    src/util/Fitting.h
    src/util/Geometry.h
    src/util/HDF5Util.cpp
    src/util/HDF5Util.h
    src/util/HilbertEnvelope.h
    src/util/ImageGrid.h
    src/util/IndexValue.h
    src/util/Interpolator.h
    src/util/Interpolator4X.h
    src/util/IterationCounter.cpp
    src/util/IterationCounter.h
    src/util/KaiserWindow.h
    src/util/KeyValueFileReader.cpp
    src/util/KeyValueFileReader.h
    src/util/LinearInterpolator.h
    src/util/Log.cpp
    src/util/Log.h
    src/util/LogSyntaxHighlighter.cpp
    src/util/LogSyntaxHighlighter.h
    src/util/Matrix.h
    src/util/MeasurementList.h
    src/util/OGL.h
    src/util/ParameterMap.cpp
    src/util/ParameterMap.h
    src/util/PseudorandomNumberGenerator.cpp
    src/util/PseudorandomNumberGenerator.h
    src/util/RealToComplexFFT.h
    src/util/SIMD.h
    src/util/Statistics.h
    src/util/Stream.h
    src/util/String.h
    src/util/TemplateUtil.h
    src/util/Tensor3.h
    src/util/Timer.h
    src/util/Util.cpp
    src/util/Util.h
    src/util/Value.h
    src/util/Visualization.cpp
    src/util/Visualization.h
    src/util/Waveform.h
    src/util/WavefrontObjFileWriter.h
    src/util/WindowFunction.h
    src/util/XY.h
    src/util/XYZ.h
    src/util/XYZValue.h
    src/util/XYZValueArray.h
    src/util/XYZValueFactor.h
    src/util/XZ.h
    src/util/XZComplexValue.h
    src/util/XZComplexValueFactor.h
    src/util/XZValue.h
    src/util/XZValueFactor.h
    src/external/lzf/lzf.h
    src/external/lzf/lzf_c.c
    src/external/lzf/lzf_d.c
    src/external/lzf/lzf_filter.c
    src/external/lzf/lzf_filter.h
    src/external/lzf/lzfP.h
    src/user/acquisition/DeviceSectorialScanAcquisition.h
    src/user/acquisition/SavedSTAAcquisition.h
    src/user/acquisition/SavedTnRnAcquisition.h
    src/user/acquisition/STAAcquisition.h
    src/user/acquisition/TnRnAcquisition.h
    src/user/configuration/DeviceSectorialScanConfiguration.h
    src/user/configuration/SA3DConfiguration.h
    src/user/configuration/STAConfiguration.h
    src/user/configuration/TnRnConfiguration.h
    src/user/configuration/TwoMediumSTAConfiguration.h
    src/user/method/CylinderDetectionAndFermatMethod.h
    src/user/method/DeviceSectorialScanMethod.h
    src/user/method/method_includes.h
    src/user/method/method_table.h
    src/user/method/MultiLayerImageMethod.cpp
    src/user/method/MultiLayerImageMethod.h
    src/user/method/NetworkSyncSTAMethod.h
    src/user/method/NetworkSyncSingleVirtualSourceMethod.h
    src/user/method/STA3DMethod.h
    src/user/method/STAMethod.h
    src/user/method/ShowImageMethod.cpp
    src/user/method/ShowImageMethod.h
    src/user/method/SimCircularSourceMethod.h
    src/user/method/SimRectangularSourceMethod.h
    src/user/method/SingleAcquisitionMethod.cpp
    src/user/method/SingleAcquisitionMethod.h
    src/user/method/SingleVirtualSourceMethod.h
    src/user/method/SyntheticYSingleVirtualSourceMethod.h
    src/user/method/T1R1SAFT3DMethod.h
    src/user/method/TestMethod.cpp
    src/user/method/TestMethod.h
    src/user/method/VTKFileMultiImageMethod.cpp
    src/user/method/VTKFileMultiImageMethod.h
    src/user/network_acquisition/NetworkDeviceSectorialScanAcquisition.h
    src/user/network_acquisition/NetworkSTAAcquisition.h
    src/user/network_acquisition/NetworkTnRnAcquisition.h
    src/user/processor/ArcCylindricalWaveProcessor.h
    src/user/processor/ArrayProcessor.h
    src/user/processor/CCBFPitchCatchProcessor.h
    src/user/processor/CrossCorrArcCylindricalWaveProcessor.h
    src/user/processor/CrossCorrelationSTAProcessor.h
    src/user/processor/CylindricalWaveProcessor.h
    src/user/processor/DefaultSTAProcessor.h
    src/user/processor/SimpleSTAProcessor.h
    src/user/processor/SynthYVectorial3DTnRnProcessor.h
    src/user/processor/TangentCurveGeometryProcessor.h
    src/user/processor/TwoMediumSTAProcessor.h
    src/user/processor/Vectorial3DSTAProcessor.h
    src/user/processor/Vectorial3DT1R1SAFTProcessor.h
    src/user/processor/Vectorial3DTnRnProcessor.h
    src/user/processor/VectorialCombinedTwoMediumImagingProcessor.h
    src/user/processor/VectorialCombinedTwoMediumImagingProcessor2.h
    src/user/processor/VectorialCombinedTwoMediumImagingProcessor3.h
    src/user/processor/VectorialCombinedTwoMediumImagingProcessor4.h
    src/user/processor/VectorialSTAProcessor.h
    src/user/processor/VectorialTwoMediumSTAProcessor.h
    src/user/simulator/AnalyticCircularSourceImpulseResponse.h
    src/user/simulator/AnalyticRectangularSourceImpulseResponse.h
    src/user/simulator/ArrayOfRectangularSourcesImpulseResponse.h
    src/user/simulator/NumericCircularSourceImpulseResponse.h
    src/user/simulator/NumericRectangularSourceImpulseResponse.h
    src/user/simulator/SimTransientAcousticField.h
    src/user/simulator/SimTransientPropagation.h
    src/user/simulator/SimTransientRadiationPattern.h
    src/user/simulator/Simulated3DAcquisitionDevice.h
    src/user/simulator/Simulated3DSTAAcquisition.h
    src/user/simulator/Simulated3DT1R1SAFTAcquisition.h
    src/user/simulator/Simulated3DTnRnAcquisition.h
    src/user/simulator/SimulatedSTAAcquisition.h

    ui/Figure2DWindow.ui
    ui/Figure3DWindow.ui
    ui/MultiLayer3DWindow.ui
    ui/USLab4a.ui

    resource/us_lab4a.qrc
)

add_executable(us_lab4a ${us_lab4a_SRC})

target_include_directories(us_lab4a PRIVATE
    src
    src/fft
    src/network_acquisition
    src/network_sync
    src/parallel
    src/util
    src/external
    src/external/lzf
    src/user/acquisition
    src/user/configuration
    src/user/method
    src/user/network_acquisition
    src/user/processor
    src/user/simulator

    ${FFTW3_INCLUDE_DIRS}
    ${FFTW3F_INCLUDE_DIRS}
    ${GLU_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(us_lab4a
    Qt5::Core
    Qt5::Gui
    Qt5::OpenGL
    Qt5::Widgets

    PkgConfig::FFTW3
    PkgConfig::FFTW3F
    PkgConfig::GLU
    PkgConfig::HDF5

    tbb
    tbbmalloc
    #tbbmalloc_proxy
    hdf5_cpp
    boost_system
    pthread
)
