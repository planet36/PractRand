# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

A fork of PractRand pre0.95 (https://pracrand.sourceforge.net/), a library and command-line suite for statistically testing random number generators.  The fork's main goal is fixing compiler and clang-tidy warnings, so most commits are small warning fixes or `.clang-tidy` adjustments.  Behavior changes to the RNGs or the statistical tests are not the point of the fork.

## Build and lint

`tools/RNG_test` is the only program the user runs.  `RNG_benchmark` and `RNG_output` are still built by default and linted, but they matter only insofar as they keep compiling cleanly.

```sh
make -j $(nproc) tools/RNG_test  # builds libPractRand.a and RNG_test only
make -j $(nproc)                 # also builds tools/RNG_benchmark and tools/RNG_output
make lint                        # clang-tidy over every library and tool source, using .clang-tidy
make clean
make install                     # copies the binaries to $(PREFIX)/bin, default /usr/local
```

The Makefile builds with `-std=c++20 -Wall -Wextra -Wpedantic -Wfatal-errors -O3 -flto=auto -march=native`.  `-Wfatal-errors` stops at the first error, so expect one diagnostic per build attempt.

To lint a single file, run clang-tidy with the same flags the `lint` target passes:

```sh
clang-tidy --quiet src/RNGs/jsf.cpp -- -MMD -MP -I include -std=c++20 -Wall -Wextra -Wpedantic
```

Any new tool must be named `tools/RNG_*.cpp` to be picked up by the Makefile's wildcard.  `tools/Test_calibration.cpp` does not match that pattern and is not built.

Do not run `./configure`.  It is upstream's script and it overwrites the hand-maintained `Makefile` with a generated POSIX one that targets `-std=c++14`.

## Testing

There is no unit-test suite.  To check a change, build and run `RNG_test`.

```sh
tools/RNG_test jsf32 -tlmax 1MB                  # test a built-in RNG by name
some_generator | tools/RNG_test stdin64          # test raw bytes piped in, read as 64-bit words
tools/RNG_test -help                             # full option list
```

A quick regression check for a warning fix is to run `tools/RNG_test <name> -seed 1234 -tlmax 16MB` before and after the change and compare the results.  The `-seed` option matters, because `RNG_test` otherwise picks a random seed.  The elapsed-time figures differ between runs and are not a regression.

## Architecture

- `include/PractRand/` and `src/` form the library, `libPractRand.a`.  `tools/` holds the command-line programs.  Most tool logic lives in headers under `tools/` (`RNG_from_name.h`, `TestManager.h`, `MultithreadedTestManager.h`, `Candidate_RNGs.h`), which each `RNG_*.cpp` includes.
- Every RNG comes in two forms.  `PractRand::RNGs::Raw::<name>` is a plain class with the state and `rawN()`, and `PractRand::RNGs::Polymorphic::<name>` wraps it as a subclass of `vRNG8`/`vRNG16`/`vRNG32`/`vRNG64`.  The macros `PRACTRAND_POLYMORPHIC_RNG_BASICS_H` and `PRACTRAND_LIGHT_WEIGHT_RNG` generate the boilerplate, so a header like `include/PractRand/RNGs/jsf32.h` is the template to copy.  The abstract base `vRNG` is in `rng_basics.h`.
- RNG state is visited through `walk_state(StateWalkingObject *)`, which drives seeding, serialization, and `print_state()`.
- RNG names given on the command line resolve through `RNG_Factories::RNG_factory_index` in `tools/RNG_from_name.h`.  That file also defines the `stdin`/`stdin8`/`stdin16`/`stdin32`/`stdin64` pseudo-RNGs that read from standard input, and it parses transform and composite names.
- Statistical tests subclass `PractRand::Tests::TestBaseclass` (`include/PractRand/tests.h`), with one header per test in `include/PractRand/Tests/`.  All of their implementations are in the single large file `src/tests.cpp`.
- `src/test_batteries.cpp` assembles tests into batteries (core and expanded) and applies "foldings", which run extra copies of the tests on sub-streams of the output.  `RNG_test`'s `-tf` and `-te` options select among these.
- `src/platform_specifics.cpp` holds the OS-dependent code for entropy and unique IDs.  `initialize_PractRand()` must be called before any threads start.

## Fixing warnings

- Fix one warning message per commit, across every place it occurs.
- Comment out unused code instead of deleting it.  Use `//` for a line or two and `#if 0` … `#endif` for a larger block, as `src/tests.cpp` already does.
- Leave the statistical arithmetic alone.  Warnings like `bugprone-integer-division` and `bugprone-incorrect-roundings` in `src/tests.cpp` and `src/math.cpp` sit in code that computes p-values, so a fix would change what `RNG_test` reports.
- Tests read the end of the previous block through negative indexes such as `data->as16[-1]`, and `TestBaseclass::get_blocks_to_repeat()` reserves room for them.  The `clang-diagnostic-array-bounds` warnings on those lines are intentional.
- The `memset(this, 0, sizeof(*this))` in the cipher RNG destructors wipes the key state.  It draws `bugprone-undefined-memory-manipulation` but is intentional.

## Code style

- Keep the upstream style: tabs for indentation and `#pragma once` header guards.  Many clang-tidy checks are disabled in `.clang-tidy` precisely because they conflict with this code's idioms (C arrays, magic numbers, pointer arithmetic, unions, and similar).  Disabling a check there is an accepted fix when it would force a rewrite of upstream code.
- Commit messages follow the existing pattern.  A `.clang-tidy` change is `Add -<check-name>`.  A warning fix is `Fix warning: <diagnostic> [<check or -W flag>]`, quoting the diagnostic with `XXX` in place of any name that varies between occurrences.
- Comments, doc blocks, and commit messages follow `COMMENT-STYLE.md`.  The rules that most often apply are one thought per sentence, two spaces after a sentence-ending period, a 96-column line limit, two hyphens for a dash in code comments, American spellings, and the serial comma.
