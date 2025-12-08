#!/bin/bash

# This script targets M1 level of performance on a 60Hz display

realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

DIR=$(dirname $(realpath $0))
cd "${DIR}/../cg2_2014_demo/prob_7/"

/usr/bin/time "${DIR}/problem_7" -screen "2560 1440 20"  -frames 333  -group_size 128 # Hz = 60 / 3
/usr/bin/time "${DIR}/problem_7" -screen "2560 1440 30"  -frames 500  -group_size 128 # Hz = 60 / 2
/usr/bin/time "${DIR}/problem_7" -screen "2560 1440 60"  -frames 1000 -group_size 128 # Hz = 60 / 1
