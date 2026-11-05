# Declarations
.PHONY: all init configure build tests run clean
all: tests

######## Settings ########
VERSION = 0.1.0
DIST=dist/charlemagne-${VERSION}.tar.gz
MAX_THREADS=20
BUILD_DIR=build
CMAKE_FLAGS=-DCHARLES_BUILD_TESTS=ON -DBUILD_EXAMPLES=OFF
OPACITY_NAME=pagn_combined
OPACITY_PATH=${CURDIR}/python/charlemagne/data/pagn_combined.sp5
OPACITY_ID=1001

######## PACKAGE ########

#### Install ####
init:
	git submodule update --recursive --init
	mkdir -p ${BUILD_DIR}

configure:
	cmake ${CMAKE_FLAGS} -S . -B ${BUILD_DIR}

build: configure
	cmake --build ${BUILD_DIR} --target charlemagne -j ${MAX_THREADS}

tests: configure
	cmake --build ${BUILD_DIR} --target charlemagne_unit_tests -j ${MAX_THREADS}
	cd ${BUILD_DIR} && ctest --output-on-failure -R '^charlemagne\.unit\.'

run: build
	CHARLEMAGNE_OPACITY_NAME=${OPACITY_NAME} \
		CHARLEMAGNE_OPACITY_PATH=${OPACITY_PATH} \
		CHARLEMAGNE_OPACITY_ID=${OPACITY_ID} \
		./build/src/charlemagne

#### Python ####
version: clean
	sed --in-place "s/ version: .*/ version: '${VERSION}',/g" meson.build
	echo "__version__ = '${VERSION}'" > python/charlemagne/__version__.py

######## Tests ########
tests: build

######## Cleanup #########
clean:
	rm -rf ${BUILD_DIR}/*
