# COS 310 Topics in Computer Science: Intro to Computer Graphics ‒ extra material

Building
========

Clone [CG2 2014 demo](https://github.com/ChaosGroup/cg2_2014_demo) into the parent directory of this repo:

```
$ cd ..
$ git clone https://github.com/ChaosGroup/cg2_2014_demo
$ cd -
```

Invoke the build script from this repo:

```
$ ./build.sh
```

Invoke executable with `CWD` set to `../cg2_2014_demo/prob_7` as:

```
$ TO_EXE=$PWD
$ cd ../cg2_2014_demo/prob_7
$ $TO_EXE/problem_7 -help
usage: /path/to/problem_7 [<option> ...]
options (multiple args to an option must constitute a single string, eg. -foo "a b c"):
        -report_caps                    : report CL capabilities
        -discard_platform_version       : discard advertised platform version when producing platform report
        -discard_device_version         : discard advertised device version when producing device report
        -platform <index>               : use platform of specified index
        -device <index>                 : use device of specified index
        -use_images                     : use images instead of buffers as source arguments
        -report_kernel_time             : report CL kernel time
        -screen <width> <height> <Hz>   : set framebuffer of specified geometry and refresh
        -frames <unsigned_integer>      : set number of frames to run; default is max unsigned int
        -group_size <positive_integer>  : set workgroup size; must be even; default is CL_KERNEL_WORK_GROUP_SIZE
```

Reference Performance
=====================

Metal frontend and 128-strong OCL workgroup

| device                     | resolution @ Hz   |
| -------------------------- | ----------------- |
| Apple M1 (7-core GPU)      | 2560 x 1312 @ 60  |
| Apple M1 (8-core GPU)      | 2560 x 1440 @ 60  |
| Apple M2 Max (30-core GPU) | 2880 x 1800 @ 120 |
