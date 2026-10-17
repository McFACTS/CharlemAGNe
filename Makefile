# Declarations
# .PHONY all clean build tests
#
all: clean build tests

######## PACKAGE ########

#### Install ####
init:
	git submodule update --recursive --init
	mkdir -p build

build: clean
	cd build; cmake ..; make

run: build
	./build/src/charlemagne

######## Tests ########
tests: build

######## Cleanup #########
clean:
	rm -rf build/*
