#!/bin/bash

# This script targets M1 Max level of performance on a 120Hz display

cd problem_7.app/Contents/MacOS

/usr/bin/time ./problem_7 -screen "2672 1600 20"  -frames 333  -group_size 128 # Hz = 120 / 6
/usr/bin/time ./problem_7 -screen "2672 1600 24"  -frames 400  -group_size 128 # Hz = 120 / 5
/usr/bin/time ./problem_7 -screen "2672 1600 30"  -frames 500  -group_size 128 # Hz = 120 / 4
/usr/bin/time ./problem_7 -screen "2672 1600 40"  -frames 666  -group_size 128 # Hz = 120 / 3
/usr/bin/time ./problem_7 -screen "2672 1600 60"  -frames 1000 -group_size 128 # Hz = 120 / 2
/usr/bin/time ./problem_7 -screen "2672 1600 120" -frames 2000 -group_size 128 # Hz = 120 / 1

cd -
