#!/bin/sh

EXEC = mummer

SOURCE = main.cpp
OBJECTS = main.o

CXXARGS = -march=native -O2 -Wall -pipe -lfasta -lsuffixarray -lpthread -D DEBUG

all: $(EXEC)

$(EXEC) : $(OBJECTS)
	g++ $< $(CXXARGS) -o $@

$(OBJECTS) : $(SOURCE)
	g++ $< $(CXXARGS) -c -o $@
	

install : $(EXEC)
	mv mummer /usr/bin/

clean:
	rm -f $(EXEC) $(OBJECTS)