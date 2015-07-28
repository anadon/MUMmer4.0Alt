#!/bin/sh

EXEC = mummer

SOURCE = main.cpp 

CXXARGS = -g -Wall -pipe -lfasta -lpthread --std=gnu++0x 

all: $(EXEC)

$(EXEC) : $(SOURCE)
	g++ $(CXXARGS) $< -o $@

clean:
	rm -f $(EXEC)