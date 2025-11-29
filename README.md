# COS 310 Topics in Computer Science: Intro to Computer Graphics ‒ extra material

Building
========

Clone [CG2 2014 demo](https://github.com/ChaosGroup/cg2_2014_demo) into the parent directory of this repo:

```
$ cd ..
$ git clone https://github.com/ChaosGroup/cg2_2014_demo
```

Invoke the build script from this repo:

```
$ cd -
$ ./build.sh
```

Build artifacts are in `problem_7.app/Contents/MacOS`. Invoke executable as:

```
$ cd problem_7.app/Contents/MacOS
$ ./problem_7
```

Reference Performance
=====================

Metal frontend and 128-strong OCL workgroup

| device                     | resolution @ Hz   |
| -------------------------- | ----------------- |
| Apple M1 (7-core GPU)      | 2560 x 1312 @ 60  |
| Apple M1 (8-core GPU)      | 2560 x 1440 @ 60  |
| Apple M2 Max (30-core GPU) | 2880 x 1800 @ 120 |
