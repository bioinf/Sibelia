CC = g++
OPTFLAGS = 
CFLAGS = $(OPTFLAGS) -O3 -w
CLINKER = g++ 
# LIBDIR = -L/usr/local/lib
MLIB =  -lm
INCDIR =   -I.
TRGT_DIR = ../..
TRGT = glocal
OBJECTS = glocal.o io.o rightinfluence.o leftinfluence.o score.o

.cpp.o:
	$(CC) -Wno-deprecated $(CFLAGS) $(INCDIR) -c $*.cpp

$(TRGT): $(OBJECTS)
	$(CLINKER) $(OPTFLAGS) $(OBJECTS) -o $(TRGT_DIR)/$(TRGT) $(MLIB)

clean :
	rm -f *.o ./*~ *~ core
