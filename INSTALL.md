
Installation
============

"Sibelia" is distributed both as binaries (Linux and OS X) and as source codes.
If you build "Sibelia" from source codes, you have to install it before usage.
However, if you download the binary distribution, you can use "Sibelia" without
installing it.

Building from the Source Code
-----------------------------

Building Sibelia requires:
* CMake
* GCC C++ compiler (version 4.6.0+ works fine)

The code is compatible with Linux and OS X environment.

To build and install the program under Linux, type following:

	cd build
	cmake ../src
	make
	sudo make install

If you do not have the root access or you want to install "Sibelia" to a
location other than /usr/local/bin, run CMake with -DCMAKE_INSTALL_PREFIX
option set:

	cd build
	cmake ../src -DCMAKE_INSTALL_PREFIX="<install destination>"
	make
	make install