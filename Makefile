#!/bin/sh

EXEC = mummer

SOURCE = main.cpp 

CXXARGS = -g -Wall -pipe -lfasta -lsuffixarray --std=gnu++0x

all: $(EXEC)

$(EXEC) : $(SOURCE)
	g++ $(CXXARGS) $< -o $@


install : $(EXEC)
	mv mummer /usr/bin/

clean:
	rm -f $(EXEC)