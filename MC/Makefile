CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags) -O3
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
ifeq "$(GCCNEWERTHAN47)" "1"
  CFLAGS += -std=c++11
else
  CFLAGS += -std=c++0x
endif
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)  
GLIBS = 
GLIBS += 
OBJECTS = skExtract.o tinyxml2.o gConfig.o globalConstants.o 
HEADERS = globalConstants.h tinyxml2.h gConfig.h

ALL : skExtract.exe
	@echo "Listo!"

skExtract.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o skExtract.exe $(LIBS) $(GLIBS) $(CFLAGS)

skExtract.o : skExtract.cc $(HEADERS)
	$(CPP) -c skExtract.cc -o skExtract.o $(CFLAGS)

gConfig.o : gConfig.cc $(HEADERS)
	$(CPP) -c gConfig.cc -o gConfig.o $(CFLAGS)

tinyxml2.o : tinyxml2.cpp $(HEADERS)
	$(CPP) -c tinyxml2.cpp -o tinyxml2.o $(CFLAGS)

globalConstants.o : globalConstants.cc $(HEADERS)
	$(CPP) -c globalConstants.cc -o globalConstants.o $(CFLAGS)
	
clean:
	rm -f *~ *.o *.exe
