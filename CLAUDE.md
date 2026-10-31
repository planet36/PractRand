# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

A fork of PractRand pre0.95 (https://pracrand.sourceforge.net/), a library and command-line suite for statistically testing random number generators.  The fork targets Linux and is updated to C++23.  Most commits fix compiler and clang-tidy warnings, and others fix upstream bugs or speed up `RNG_test`.  A change must not alter what an RNG outputs or what a statistical test reports, unless it fixes a bug in one.

Only Linux with GCC on little-endian targets is supported, and `include/PractRand/endian.h` enforces little-endian with a `static_assert`.  The Windows, MSVC, Apple, and big-endian code was removed, so do not add branches for them back.

## Build and lint

`tools/RNG_test` is the only program the user runs.  `RNG_benchmark` and `RNG_output` are still built by default and linted, but they matter only insofar as they keep compiling cleanly.

```sh
make -j $(nproc) tools/RNG_test  # builds libPractRand.a and RNG_test only
make -j $(nproc)                 # also builds tools/RNG_benchmark and tools/RNG_output
make lint                        # clang-tidy over every source the build compiles, using .clang-tidy
make clean
make install                     # copies the binaries to $(PREFIX)/bin, default /usr/local
```

The Makefile's `CXXFLAGS` include `-Wfatal-errors`, which stops at the first error, so expect one diagnostic per build attempt.

To lint a single file, override the `lint` target's source lists.  That keeps clang-tidy's flags the same as the Makefile's.

```sh
make lint LIB_SRCS=src/RNGs/jsf.cpp BIN_SRCS=
make lint LIB_SRCS= BIN_SRCS=tools/RNG_test.cpp
```

After an LLVM upgrade, run `clang-tidy --verify-config` from the repo root.  clang-tidy silently ignores a check name in `.clang-tidy` that no longer exists, so a renamed check comes back on without notice.

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

- `include/PractRand/` and `src/` form the library, `libPractRand.a`.  `tools/` holds the command-line programs.  Most tool logic lives in headers under `tools/`, such as `RNG_from_name.h`, `TestManager.h`, `MultithreadedTestManager.h`, and `Candidate_RNGs.h`, which the `RNG_*.cpp` files include as needed.
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
- The cipher RNG destructors wipe the key state with `explicit_bzero(this, sizeof(*this))`.  Keep it rather than `memset`, which draws `bugprone-undefined-memory-manipulation` and `clang-diagnostic-nontrivial-memcall`, and which GCC removes as a dead store.

## Code style

@COMMENT-STYLE.md

- Keep the upstream style: tabs for indentation and `#pragma once` header guards.  Many clang-tidy checks are disabled in `.clang-tidy` precisely because they conflict with this code's idioms (C arrays, magic numbers, pointer arithmetic, unions, and similar).  Disabling a check there is an accepted fix when it would force a rewrite of upstream code.
- Put local `#include "..."` lines before system `#include <...>` lines, with a blank line between the two groups.
- Attach `&` and `*` to the type, as in `const std::string& name` and `vRNG* rng`.  A line that declares several pointers, such as `uint32_t *m = nullptr, *m2 = nullptr;`, keeps the sigil on each name.
- Never remove existing braces.  Adding braces is fine.
- Use signed integers for counts, sizes, indexes, and loop counters.  Keep unsigned integers for bit manipulation, raw RNG output, and byte data.  A class member keeps its fixed width when it changes signedness, so `uint32_t` becomes `int32_t`.
- Commit messages follow the existing pattern.  Disabling a check in `.clang-tidy` is `Add -<check-name>`.  A warning fix is `Fix warning: <diagnostic> [<check or -W flag>]`, quoting the diagnostic with `XXX` in place of any name that varies between occurrences.  When more than one part varies, each gets its own placeholder (`XXX`, `YYY`, and so on), as in `Fix warning: prefer 'std::numbers::XXX' to this literal, differs by 'YYY' [modernize-use-std-numbers]`.  A warning silenced with `// NOLINTNEXTLINE(<check>)` is `Suppress warning: <diagnostic> [<check>]`, with the same placeholders.
