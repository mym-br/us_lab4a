cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

cmake_policy(SET CMP0104 NEW)

set(CMAKE_CUDA_ARCHITECTURES "89")
set(CMAKE_CUDA_STANDARD 17)
if(LINUX)
    set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++-12)
endif()

project(us_lab4a LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

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
pkg_search_module(HDF5 REQUIRED hdf5-serial hdf5)

if(UNIX)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:C,CXX>:-Wall>
        $<$<COMPILE_LANGUAGE:C,CXX>:-Wextra>
        $<$<COMPILE_LANGUAGE:C,CXX>:-march=native>
    )
endif()

option(LAB_ENABLE_SIMD "Enable SIMD." ON)
option(LAB_ENABLE_OPENCL "Enable OpenCL." OFF)
set(LAB_OPENCL_PLATFORM 0 CACHE STRING "OpenCL platform index.")
set(LAB_OPENCL_DEVICE 0 CACHE STRING "OpenCL device index.")
option(LAB_ENABLE_CUDA "Enable CUDA." OFF)
option(LAB_ENABLE_EXECUTION_TIME_MEASUREMENT "Enable execution time measurement." OFF)

if(LAB_ENABLE_SIMD)
    add_compile_definitions(LAB_ENABLE_SIMD)
endif()
add_compile_definitions(QT_NO_KEYWORDS)

if(CMAKE_COMPILER_IS_GNUCC AND (CMAKE_CXX_COMPILER_VERSION MATCHES "^8\."))
    set(CXX_FS_LIB stdc++fs)
endif()

if(LAB_ENABLE_OPENCL)
    set(OPENCL_LIB OpenCL)
    set(OPENCL_SRC
        src/user/processor/VectorialCombinedTwoMediumImagingOCLProcessor.h
        src/user/processor/VectorialCombinedTwoMediumImagingOCLProcessor2.h
        src/user/processor/VectorialSTAOCLProcessor.h
        src/user/simulator/NumericArrayOfRectangularSourcesOCLImpulseResponse.h
        src/user/simulator/NumericCircularSourceOCLImpulseResponse.h
        src/user/simulator/NumericRectangularSourceOCLImpulseResponse.h
        src/util/OCLAtomic.h
        src/util/OCLCoherenceFactor.h
        src/util/OCLGeometry.h
        src/util/OCLReduce.h
        src/util/OCLStatistics.h
        src/util/OCLUtil.h
    )
    add_compile_definitions(
        LAB_ENABLE_OPENCL
    )
    set(OPENCL_TARGET_DEFS
        CL_HPP_ENABLE_EXCEPTIONS
        CL_HPP_MINIMUM_OPENCL_VERSION=120
        CL_HPP_TARGET_OPENCL_VERSION=120
        CL_HPP_CL_1_2_DEFAULT_BUILD
        LAB_OPENCL_PLATFORM=${LAB_OPENCL_PLATFORM}
        LAB_OPENCL_DEVICE=${LAB_OPENCL_DEVICE}
        LAB_OPENCL_PROGRAM_BUILD_OPTIONS="-cl-std=CL1.2"
    )
endif()

if(LAB_ENABLE_EXECUTION_TIME_MEASUREMENT)
    add_compile_definitions(LAB_ENABLE_EXECUTION_TIME_MEASUREMENT)
endif()

set(uslab4a_SRC
    src/Method.cpp
    src/Method.h
    src/Project.h
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
    src/util/Math.cpp
    src/util/Math.h
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
    src/util/ThreadUtil.cpp
    src/util/ThreadUtil.h
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
    src/user/method/CylinderDetectionAndFermatMethod.cpp
    src/user/method/CylinderDetectionAndFermatMethod.h
    src/user/method/DeviceSectorialScanMethod.cpp
    src/user/method/DeviceSectorialScanMethod.h
    src/user/method/method_includes.h
    src/user/method/method_table.h
    src/user/method/MultiLayerImageMethod.cpp
    src/user/method/MultiLayerImageMethod.h
    src/user/method/NetworkSyncSTAMethod.cpp
    src/user/method/NetworkSyncSTAMethod.h
    src/user/method/NetworkSyncSingleVirtualSourceMethod.cpp
    src/user/method/NetworkSyncSingleVirtualSourceMethod.h
    src/user/method/ShowImageMethod.cpp
    src/user/method/ShowImageMethod.h
    src/user/method/SimCircularSourceMethod.cpp
    src/user/method/SimCircularSourceMethod.h
    src/user/method/SimRectangularSourceMethod.cpp
    src/user/method/SimRectangularSourceMethod.h
    src/user/method/SingleAcquisitionMethod.cpp
    src/user/method/SingleAcquisitionMethod.h
    src/user/method/SingleVirtualSourceMethod.cpp
    src/user/method/SingleVirtualSourceMethod.h
    src/user/method/STA3DMethod.cpp
    src/user/method/STA3DMethod.h
    src/user/method/STAMethod.cpp
    src/user/method/STAMethod.h
    src/user/method/SyntheticYSingleVirtualSourceMethod.cpp
    src/user/method/SyntheticYSingleVirtualSourceMethod.h
    src/user/method/T1R1SAFT3DMethod.cpp
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
    src/user/processor/VectorialSTAProcessor.h
    src/user/processor/VectorialTwoMediumSTAProcessor.h
    src/user/simulator/AnalyticCircularSourceImpulseResponse.h
    src/user/simulator/AnalyticRectangularSourceImpulseResponse.h
    src/user/simulator/ArrayOfRectangularSourcesImpulseResponse.h
    src/user/simulator/NumericArrayOfRectangularSourcesImpulseResponse.h
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

    ${OPENCL_SRC}
)

set(uslab4a_gui_SRC
    src/Controller.cpp
    src/Controller.h
    src/Project.cpp
    src/Figure2DWidget.cpp
    src/Figure2DWidget.h
    src/Figure2DWindow.cpp
    src/Figure2DWindow.h
    src/Figure3DWindow.cpp
    src/Figure3DWindow.h
    src/FigureWindowList.h
    src/LogSyntaxHighlighter.cpp
    src/LogSyntaxHighlighter.h
    src/main.cpp
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
    src/USLab4a.cpp
    src/USLab4a.h

    ui/Figure2DWindow.ui
    ui/Figure3DWindow.ui
    ui/MultiLayer3DWindow.ui
    ui/USLab4a.ui

    resource/us_lab4a.qrc
)

if(LAB_ENABLE_CUDA)
    enable_language(CUDA)
    #if(CMAKE_CUDA_FLAGS MATCHES "^ *$")
    #    set(CMAKE_CUDA_FLAGS "-arch compute_61 -code sm_61,compute_61" CACHE STRING "Flags used by the CUDA compiler during all build types." FORCE)
    #endif()
    add_compile_definitions(
        LAB_ENABLE_CUDA
    )
    add_library(usercuda STATIC
        src/user/processor/VectorialCombinedTwoMediumImagingCUDAProcessor.cpp
        src/user/processor/VectorialCombinedTwoMediumImagingCUDAProcessor.cu
        src/user/processor/VectorialCombinedTwoMediumImagingCUDAProcessor.h
        src/user/processor/VectorialCombinedTwoMediumImagingCUDAProcessor2.cpp
        src/user/processor/VectorialCombinedTwoMediumImagingCUDAProcessor2.cu
        src/user/processor/VectorialCombinedTwoMediumImagingCUDAProcessor2.h
        src/user/processor/VectorialSTACUDAProcessor.cpp
        src/user/processor/VectorialSTACUDAProcessor.cu
        src/user/processor/VectorialSTACUDAProcessor.h
        src/user/simulator/NumericArrayOfRectangularSourcesCUDAImpulseResponse.cu
        src/user/simulator/NumericArrayOfRectangularSourcesCUDAImpulseResponse.h
        src/user/simulator/NumericCircularSourceCUDAImpulseResponse.cu
        src/user/simulator/NumericCircularSourceCUDAImpulseResponse.h
        src/user/simulator/NumericRectangularSourceCUDAImpulseResponse.cu
        src/user/simulator/NumericRectangularSourceCUDAImpulseResponse.cuh
        src/user/simulator/NumericRectangularSourceCUDAImpulseResponse.h
        src/util/CUDACoherenceFactor.cuh
        src/util/CUDAGeometry.cuh
        src/util/CUDAReduce.cu
        src/util/CUDAReduce.cuh
        src/util/CUDAStatistics.cuh
        src/util/CUDAUtil.h
    )
    set(USER_CUDA_LIB usercuda)
    target_include_directories(usercuda PRIVATE
        src/fft
        src/util
        src/user/acquisition
        src/user/configuration
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
    set_target_properties(usercuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON)
    target_compile_features(usercuda PRIVATE cxx_std_17)
    target_compile_options(usercuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--compiler-options=-march=native>
                                            $<$<COMPILE_LANGUAGE:CUDA>:--compiler-options=-std=c++17>)
endif()

add_library(uslab4a STATIC ${uslab4a_SRC})

target_compile_definitions(uslab4a PRIVATE ${OPENCL_TARGET_DEFS})

target_include_directories(uslab4a PUBLIC
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
    ${HDF5_INCLUDE_DIRS}
)

add_executable(us_lab4a_gui ${uslab4a_gui_SRC})

target_link_libraries(us_lab4a_gui
    uslab4a

    Qt5::Core
    Qt5::Gui
    Qt5::OpenGL
    Qt5::Widgets

    PkgConfig::FFTW3
    PkgConfig::FFTW3F
    PkgConfig::GLU

    ${OPENCL_LIB}
    ${USER_CUDA_LIB}

    ${HDF5_LDFLAGS}
    hdf5_cpp
    hdf5

    ${CXX_FS_LIB}

    tbb
    tbbmalloc
    #tbbmalloc_proxy
    boost_system
    pthread
    GL
)

target_include_directories(us_lab4a_gui PUBLIC
    ${GLU_INCLUDE_DIRS}
)

add_executable(us_lab4a
    src/cli_main.cpp
    src/cli_Project.cpp
    src/cli_USLab4a.h
)

target_link_libraries(us_lab4a
    uslab4a

    PkgConfig::FFTW3
    PkgConfig::FFTW3F

    ${OPENCL_LIB}
    ${USER_CUDA_LIB}

    ${HDF5_LDFLAGS}
    hdf5_cpp
    hdf5

    ${CXX_FS_LIB}

    tbb
    tbbmalloc
    #tbbmalloc_proxy
    boost_system
    pthread
)
