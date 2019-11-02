
us_lab4a
========

us_lab4a is a tool for acoustic field simulation and imaging.

WARNING
-------

This is an experimental software.
Use at your own risk!

Features
--------

- Acoustic field simulation of circular sources, rectangular sources and
  arrays of rectangular sources (analytic and numeric methods).
  - Transient acoustic field.
  - Impulse response.
  - Transient propagation.
  - Radiation pattern.
- Simulated and experimental imaging using 1D or 2D arrays.
  - Synthetic transmit aperture.
  - SAFT (one transmitter and one receiver in each acquisition).
  - Single virtual source.
- Simple 2D and 3D figures.

Supported systems
-----------------

The software has been tested in Debian 10 x86_64.

The following libraries are required:

- Boost
- FFTW
- HDF5
- Qt 5
- TBB - Threading Building Blocks

License
-------

C++ code licensed under the GNU GPL 3 or later.

The scripts are in the public domain.

The .txt and .h5 files are in the public domain.

External code
-------------

The LZF filter for HDF5 is provided by the [h5py][] project
(files in src/external/lzf).

[h5py]: https://www.h5py.org/
