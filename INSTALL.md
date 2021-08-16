
Requirements
------------

- CMake
- Boost
- FFTW (single and double precision)
- HDF5 (with C++ support)
- Qt 5.15
- TBB - Threading Building Blocks
- pkg-config
- GLU

A C++ compiler with support for C++17 (or newer specification) is required to
build the software.

The software has been tested in Debian 11 (Bullseye) x86_64. For Debian 11, the
following command installs the dependencies:

    apt install cmake g++ pkg-config \
      qt5-qmake qtbase5-dev \
      libtbb-dev libhdf5-dev libboost-dev libboost-system-dev libfftw3-dev

Dependencies for CUDA:

- NVIDIA closed source drivers.
- CUDA Toolkit.

Dependencies for OpenCL:

- The OpenCL library, headers and at least one installable client driver (ICD).

Building (Linux+GNU)
--------------------

    mkdir ../us_lab4a-build
    cd ../us_lab4a-build
    cmake -D CMAKE_BUILD_TYPE=Release ../us_lab4a
    cmake --build .

To enable CUDA, add
` -D LAB_ENABLE_CUDA=ON -D CMAKE_CUDA_FLAGS='-arch compute_75 -code sm_75' `
after "Release" in the cmake command above. Replace "75" with the compute
capability of your GPU (without the dot).
See https://developer.nvidia.com/cuda-gpus.

To enable OpenCL, add
` -D LAB_ENABLE_OPENCL=ON `
after "Release" in the cmake command above.

Running (Linux+GNU)
-------------------

    ./us_lab4a_gui &

- Click on "Select", choose a directory under `../us_lab4a/project`
  (e.g. sim_rectangular_source-single).
- Select a task (e.g. acoustic_field-transient-analytic).
- Select an experiment (e.g. config_1).
- Click on "Enable".
