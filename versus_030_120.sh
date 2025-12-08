#!/bin/bash

# This script targets M1 Max level of performance on a 120Hz display

realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

DIR=$(dirname $(realpath $0))
cd "${DIR}/../cg2_2014_demo/prob_7/"

# Juxtapose 30Hz vs 120Hz
# Run former for a 1/10th duration (16.6 sec) run, latter for a full duration (166.6 sec) run
/usr/bin/time "${DIR}/problem_7" -screen "2672 1600 30"  -frames 500   -group_size 128 # Hz = 120 / 4
/usr/bin/time "${DIR}/problem_7" -screen "2672 1600 120" -frames 20000 -group_size 128 # Hz = 120 / 1
