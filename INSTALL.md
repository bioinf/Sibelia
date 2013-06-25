Installation
============

"Sibelia" is distributed both as binaries (Linux and OS X) and as source codes.
If you build "Sibelia" from source codes, you have to install it before usage.

System requirements
-------------------
"C-Sibelia" requires "Python" and "Perl" to be istalled in your system. Please
note that "C-Sibelia" runs only under "Python2" of version at least 2.7. The
annotation script requires "Java" runtime. Unfortunately, "C-Sibelia" is 
available only under "Linux", 

Building from the Source Code
-----------------------------

Building Sibelia requires:
* CMake
* GCC C++ compiler (version 4.2.0+ works fine)

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

You can install "Sibelia" only, without "C-Sibelia". To do this, use the option
"NOCSIBELIA" while running "CMake":

	cd build
	cmake ../src -NOCSIBELIA
	make
	sudo make install
