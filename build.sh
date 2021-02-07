#!/bin/bash

# alternatively, via xcode native cli interface (artifacts go to ./build/Release/GLEssentials.app/):
#
# $ xcodebuild -project GLEssentials.xcodeproj -configuration Release -target GLEssentials-OSX

shopt -s nullglob

EXTERN_PROB_PATH=../cg2_2014_demo/prob_7
EXTERN_COMMON_PATH=../cg2_2014_demo/common
MAKE_OUTPUT_PATH=problem_7.app/Contents/MacOS

if [[ $1 == "clean" ]]; then
	HOSTTYPE=${HOSTTYPE} LANG=en_US.US-ASCII EXTERN_COMMON_PATH=${EXTERN_COMMON_PATH} EXTERN_PROB_PATH=${EXTERN_PROB_PATH} make clean

	if [ -d ${MAKE_OUTPUT_PATH} ]; then
		FILES=(${MAKE_OUTPUT_PATH}/*.glsl? ${MAKE_OUTPUT_PATH}/kernel)

		if [[ ${#FILES[@]} -ne 0 ]]; then
			rm -rf ${FILES[@]}
		fi
		rmdir ${MAKE_OUTPUT_PATH}
	fi

	exit 0
fi

if [ ! -d ${MAKE_OUTPUT_PATH} ]; then
	mkdir ${MAKE_OUTPUT_PATH}
fi

HOSTTYPE=${HOSTTYPE} LANG=en_US.US-ASCII EXTERN_COMMON_PATH=${EXTERN_COMMON_PATH} EXTERN_PROB_PATH=${EXTERN_PROB_PATH} make

if (( $? == 0 )); then
	FILES=(${EXTERN_PROB_PATH}/*.glsl? ${EXTERN_PROB_PATH}/kernel)

	if [[ ${#FILES[@]} -ne 0 ]]; then
		cp -r ${FILES[@]} ${MAKE_OUTPUT_PATH}
	fi
fi
