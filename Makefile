FITSIOROOT = /home/babieca/MEGA/MEGAsync/CCD/Datos_Fe55/cfitsio/
CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS +=
ifeq "$(GCCNEWERTHAN47)" "1"
  CFLAGS += -std=c++11
else
  CFLAGS += -std=c++0x
endif
OBJECTS =  toyMC.o 
#HEADERS = globalConstants.h

ALL : toyMC.exe
	@echo "Listo!"

toyMC.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o toyMC.exe $(LIBS) $(GLIBS) $(CFLAGS)

toyMC.o : toyMC.C $(HEADERS)
	$(CPP) -c toyMC.C -o toyMC.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe


