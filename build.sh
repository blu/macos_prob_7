#!/bin/bash

# alternatively, via xcode native cli interface (artifacts go to ./build/Release/GLEssentials.app/):
#
# $ xcodebuild -project GLEssentials.xcodeproj -configuration Release -target GLEssentials-OSX

EXTERN_PROB_PATH=../interview/prob_7
EXTERN_COMMON_PATH=../interview/common
MAKE_OUTPUT_PATH=problem_7.app/Contents/MacOS

if [[ $1 == "clean" ]]; then
	HOSTTYPE=${HOSTTYPE} LANG=en_US.US-ASCII EXTERN_COMMON_PATH=${EXTERN_COMMON_PATH} EXTERN_PROB_PATH=${EXTERN_PROB_PATH} make clean

	if [ -d ${MAKE_OUTPUT_PATH} ]; then
		rm ${MAKE_OUTPUT_PATH}/*.glsl?
		rmdir ${MAKE_OUTPUT_PATH}
	fi

	exit 0
fi

if [ ! -d ${MAKE_OUTPUT_PATH} ]; then
	mkdir ${MAKE_OUTPUT_PATH}
fi

HOSTTYPE=${HOSTTYPE} LANG=en_US.US-ASCII EXTERN_COMMON_PATH=${EXTERN_COMMON_PATH} EXTERN_PROB_PATH=${EXTERN_PROB_PATH} make

if (( $? == 0 )); then
	cp -v ${EXTERN_PROB_PATH}/*.glsl? ${MAKE_OUTPUT_PATH}
fi
