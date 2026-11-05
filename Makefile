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

### Python build things ###
WHL=wheelhouse
LOCAL_WHEELS=local-wheels
FREETHREADED := $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('Py_GIL_DISABLED') or 0)")
ifeq ($(FREETHREADED),1)
	NOLIMITED	= -Dpython.allow_limited_api=false
	BUILD_OPTS 	= -C setup-args=$(NOLIMITED)
	PIP_OPTS	= --config-settings=setup-args=$(NOLIMITED)
endif

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

$(DIST): clean build version
	- git add python meson.build __version__.py
	- git commit
	python3 -m build $(BUILD_OPTS)

install-dist: $(DIST)
	python3 -m pip install $(DIST) $(PIP_OPTS)

install: clean version
	pip install . $(PIP_OPTS)

######## Tests ########
#tests: build

######## Cleanup #########
clean:
	rm -rf ${BUILD_DIR}/*
