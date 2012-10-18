Installation
============

Sibelia is distributed both as binaries (Windows, Linux and OS X)
and as source codes. 

Building from the Source Code
-----------------------------

Building Sibelia requires:
* CMake
* GCC C++ compiler (version 4.2.0+ works fine)

The code is compatible with MSVS 2010 and MinGW,
as well as with default Linux and OS X environment.

To build and install the program, type following:

	cd build
	cmake ../src
	make install

If you want to install "Sibelia" to a location other than /usr/local/bin,
run CMake with -DCMAKE_INSTALL_PREFIX option set:

	cmake ../src -DCMAKE_INSTALL_PREFIX="<install destination>"