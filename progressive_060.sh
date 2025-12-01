#!/bin/bash

# This script targets M1 level of performance on a 60Hz display

cd problem_7.app/Contents/MacOS

/usr/bin/time ./problem_7 -screen "2560 1440 20"  -frames 333  -group_size 128 # Hz = 60 / 3
/usr/bin/time ./problem_7 -screen "2560 1440 30"  -frames 500  -group_size 128 # Hz = 60 / 2
/usr/bin/time ./problem_7 -screen "2560 1440 60"  -frames 1000 -group_size 128 # Hz = 60 / 1

cd -
