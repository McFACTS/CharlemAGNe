# Declarations
# .PHONY all clean build tests
#
all: clean build tests

######## Settings ########
MAX_THREADS=20

######## PACKAGE ########

#### Install ####
init:
	git submodule update --recursive --init
	mkdir -p build

build: clean
	cd build; cmake ..; make

tests: clean
	cd build && cmake -DCHARLES_BUILD_TESTS=ON .. && make -j ${MAX_THREADS} && make test

run: build
	./build/src/charlemagne

######## Tests ########
tests: build

######## Cleanup #########
clean:
	rm -rf build/*
