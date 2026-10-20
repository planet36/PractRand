# PractRand

This is a fork of [PractRand-pre0.95](https://pracrand.sourceforge.net/), a library and command-line suite for statistically testing random number generators.  The fork targets Linux and is updated to C++23.  It fixes many compiler and clang-tidy warnings and several upstream bugs, and its `RNG_test` runs faster.

## Requirements

- Linux on a little-endian CPU
- GCC with C++23 support, including `<print>` (tested with GCC 16)
- GNU Make

## Build

```sh
make -j $(nproc)
make install PREFIX=$HOME/.local
```

The Makefile builds with `-march=native`, so the binaries may not run on a different CPU.  Do not run `./configure`.  It is upstream's script, and it overwrites the Makefile.

## Usage

```sh
tools/RNG_test jsf32 -multithreaded -tlmax 1GB
some_generator | tools/RNG_test stdin64 -multithreaded
tools/RNG_test -help
```

[How to Test with PractRand](https://www.pcg-random.org/posts/how-to-test-with-practrand.html) explains how to read the results.
