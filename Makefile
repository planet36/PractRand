CXXFLAGS +=  -std=c++20 -Iinclude -pthread -DMULTITHREADING_SUPPORTED=1

# if optimized
CXXFLAGS += -O3 -march=native -mtune=native

# if multithreaded
# CXXFLAGS += -pthread -DMULTITHREADING_SUPPORTED=1
LDFLAGS += -pthread

LIB_PREFIX ?= lib
LIB_SUFFIX ?= .a
SHLIB_LDFLAGS ?= -shared -fPIC

LIB_NAME = PractRand

OS ?= $(shell uname -s)

ifeq ($(SHLIB_SUFFIX),)
ifeq ($(OS),Darwin)
	SHLIB_SUFFIX=.dylib
else
	SHLIB_SUFFIX=.so
endif
endif

CXXFLAGS += -O3 -flto=auto -march=native -mtune=native -fno-math-errno

STATIC_LIB=$(LIB_PREFIX)$(LIB_NAME)$(LIB_SUFFIX)
SHARED_LIB=$(LIB_PREFIX)$(LIB_NAME)$(SHLIB_SUFFIX)

ifeq ($(OS),Darwin)
	SHLIB_LDFLAGS += -undefined dynamic_lookup
else
	SHLIB_LDFLAGS += -Wl,-soname,$(SHARED_LIB)
endif


all: RNG_test RNG_benchmark $(STATIC_LIB)

.SUFFIXES: .cpp .o $LIB_SUFFIX $SHLIB_SUFFIX

src/math.o: src/math.cpp include/PractRand.h include/PractRand/config.h \
  include/PractRand/rng_basics.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h include/PractRand/RNGs/efiix64x384.h \
  include/PractRand/rng_internals.h include/PractRand/test_helpers.h
src/non_uniform.o: src/non_uniform.cpp include/PractRand/config.h \
  include/PractRand/rng_basics.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h include/PractRand/rng_internals.h \
  include/PractRand/RNGs/all.h include/PractRand/RNGs/arbee.h \
  include/PractRand/RNGs/chacha.h include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/RNGs/hc256.h include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/RNGs/jsf32.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/RNGs/mt19937.h \
  include/PractRand/RNGs/rarns16.h include/PractRand/RNGs/rarns32.h \
  include/PractRand/RNGs/rarns64.h include/PractRand/RNGs/salsa.h \
  include/PractRand/RNGs/sfc16.h include/PractRand/RNGs/sfc32.h \
  include/PractRand/RNGs/sfc64.h \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/RNGs/trivium.h include/PractRand/RNGs/xsm32.h \
  include/PractRand/RNGs/xsm64.h
src/platform_specifics.o: src/platform_specifics.cpp \
  include/PractRand/config.h include/PractRand/rng_basics.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h
src/rand.o: src/rand.cpp include/PractRand/config.h \
  include/PractRand/rng_basics.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h include/PractRand/rng_internals.h \
  include/PractRand/RNGs/all.h include/PractRand/RNGs/arbee.h \
  include/PractRand/RNGs/chacha.h include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/RNGs/hc256.h include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/RNGs/jsf32.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/RNGs/mt19937.h \
  include/PractRand/RNGs/rarns16.h include/PractRand/RNGs/rarns32.h \
  include/PractRand/RNGs/rarns64.h include/PractRand/RNGs/salsa.h \
  include/PractRand/RNGs/sfc16.h include/PractRand/RNGs/sfc32.h \
  include/PractRand/RNGs/sfc64.h \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/RNGs/trivium.h include/PractRand/RNGs/xsm32.h \
  include/PractRand/RNGs/xsm64.h
src/sha2.o: src/sha2.cpp include/PractRand/config.h \
  include/PractRand/endian.h include/PractRand/sha2.h
src/test_batteries.o: src/test_batteries.cpp \
  include/PractRand/Tests/BCFN.h include/PractRand/test_helpers.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/tests.h include/PractRand/Tests/BCFN_MT.h \
  include/PractRand/Tests/BRank.h include/PractRand/Tests/Birthday.h \
  include/PractRand/Tests/CoupGap.h include/PractRand/Tests/DistC6.h \
  include/PractRand/Tests/DistFreq4.h include/PractRand/Tests/FPF.h \
  include/PractRand/Tests/FPMulti.h include/PractRand/Tests/Gap16.h \
  include/PractRand/Tests/NearSeq.h include/PractRand/Tests/Pat5.h \
  include/PractRand/Tests/coup16.h include/PractRand/Tests/mod3.h \
  include/PractRand/Tests/transforms.h \
  include/PractRand/test_batteries.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h
src/tests.o: src/tests.cpp include/PractRand/Tests/BCFN.h \
  include/PractRand/test_helpers.h include/PractRand/rng_basics.h \
  include/PractRand/config.h include/PractRand/tests.h \
  include/PractRand/Tests/BCFN_MT.h include/PractRand/Tests/BRank.h \
  include/PractRand/Tests/Birthday.h include/PractRand/Tests/CoupGap.h \
  include/PractRand/Tests/DistC6.h include/PractRand/Tests/DistFreq4.h \
  include/PractRand/Tests/FPF.h include/PractRand/Tests/FPMulti.h \
  include/PractRand/Tests/Gap16.h include/PractRand/Tests/NearSeq.h \
  include/PractRand/Tests/Pat5.h include/PractRand/Tests/coup16.h \
  include/PractRand/Tests/mod3.h include/PractRand/Tests/transforms.h \
  include/PractRand/test_batteries.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h include/PractRand/rng_internals.h
src/RNGs/arbee.o: src/RNGs/arbee.cpp include/PractRand/RNGs/arbee.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/endian.h include/PractRand/rng_internals.h
src/RNGs/chacha.o: src/RNGs/chacha.cpp include/PractRand/RNGs/chacha.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h
src/RNGs/efiix.o: src/RNGs/efiix.cpp include/PractRand/RNGs/arbee.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/rng_internals.h
src/RNGs/hc256.o: src/RNGs/hc256.cpp include/PractRand/RNGs/hc256.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h
src/RNGs/isaac.o: src/RNGs/isaac.cpp include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/rng_internals.h
src/RNGs/jsf.o: src/RNGs/jsf.cpp include/PractRand/RNGs/jsf32.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/rng_internals.h
src/RNGs/mt19937.o: src/RNGs/mt19937.cpp include/PractRand/RNGs/jsf32.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/mt19937.h include/PractRand/rng_internals.h
src/RNGs/rarns.o: src/RNGs/rarns.cpp include/PractRand/RNGs/rarns16.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/rarns32.h include/PractRand/RNGs/rarns64.h \
  include/PractRand/rng_internals.h
src/RNGs/salsa.o: src/RNGs/salsa.cpp include/PractRand/RNGs/salsa.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/endian.h include/PractRand/rng_internals.h
src/RNGs/sfc.o: src/RNGs/sfc.cpp include/PractRand/RNGs/sfc16.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/sfc32.h include/PractRand/RNGs/sfc64.h \
  include/PractRand/rng_internals.h
src/RNGs/sha2_based_pool.o: src/RNGs/sha2_based_pool.cpp \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h include/PractRand/sha2.h
src/RNGs/trivium.o: src/RNGs/trivium.cpp include/PractRand/RNGs/trivium.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/endian.h include/PractRand/rng_internals.h
src/RNGs/xsm.o: src/RNGs/xsm.cpp include/PractRand/config.h \
  include/PractRand/rng_basics.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h include/PractRand/rng_internals.h \
  include/PractRand/RNGs/xsm32.h include/PractRand/RNGs/xsm64.h
src/RNGs/other/fibonacci.o: src/RNGs/other/fibonacci.cpp \
  include/PractRand/RNGs/mt19937.h include/PractRand/rng_basics.h \
  include/PractRand/config.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/other/fibonacci.h \
  include/PractRand/rng_internals.h
src/RNGs/other/indirection.o: src/RNGs/other/indirection.cpp \
  include/PractRand/RNGs/arbee.h include/PractRand/rng_basics.h \
  include/PractRand/config.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/other/indirection.h \
  include/PractRand/rng_internals.h
src/RNGs/other/mult.o: src/RNGs/other/mult.cpp \
  include/PractRand/RNGs/other/mult.h include/PractRand/rng_helpers.h \
  include/PractRand/config.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_basics.h include/PractRand/rng_internals.h
src/RNGs/other/rng_sets.o: src/RNGs/other/rng_sets.cpp \
  include/PractRand/RNGs/mt19937.h include/PractRand/rng_basics.h \
  include/PractRand/config.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/other/fibonacci.h \
  include/PractRand/rng_internals.h
src/RNGs/other/simple.o: src/RNGs/other/simple.cpp \
  include/PractRand/RNGs/other/simple.h include/PractRand/rng_helpers.h \
  include/PractRand/config.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_basics.h include/PractRand/rng_internals.h
src/RNGs/other/transform.o: src/RNGs/other/transform.cpp \
  include/PractRand/RNGs/other/transform.h \
  include/PractRand/rng_helpers.h include/PractRand/config.h \
  include/PractRand/rng_adaptors.h include/PractRand/rng_basics.h \
  include/PractRand/rng_internals.h include/PractRand/tests.h
src/RNGs/arbee.o: src/RNGs/arbee.cpp include/PractRand/RNGs/arbee.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/endian.h include/PractRand/rng_internals.h
src/RNGs/chacha.o: src/RNGs/chacha.cpp include/PractRand/RNGs/chacha.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h
src/RNGs/efiix.o: src/RNGs/efiix.cpp include/PractRand/RNGs/arbee.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/rng_internals.h
src/RNGs/hc256.o: src/RNGs/hc256.cpp include/PractRand/RNGs/hc256.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h
src/RNGs/isaac.o: src/RNGs/isaac.cpp include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/rng_internals.h
src/RNGs/jsf.o: src/RNGs/jsf.cpp include/PractRand/RNGs/jsf32.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/rng_internals.h
src/RNGs/mt19937.o: src/RNGs/mt19937.cpp include/PractRand/RNGs/jsf32.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/mt19937.h include/PractRand/rng_internals.h
src/RNGs/rarns.o: src/RNGs/rarns.cpp include/PractRand/RNGs/rarns16.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/rarns32.h include/PractRand/RNGs/rarns64.h \
  include/PractRand/rng_internals.h
src/RNGs/salsa.o: src/RNGs/salsa.cpp include/PractRand/RNGs/salsa.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/endian.h include/PractRand/rng_internals.h
src/RNGs/sfc.o: src/RNGs/sfc.cpp include/PractRand/RNGs/sfc16.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/sfc32.h include/PractRand/RNGs/sfc64.h \
  include/PractRand/rng_internals.h
src/RNGs/sha2_based_pool.o: src/RNGs/sha2_based_pool.cpp \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/rng_internals.h include/PractRand/sha2.h
src/RNGs/trivium.o: src/RNGs/trivium.cpp include/PractRand/RNGs/trivium.h \
  include/PractRand/rng_basics.h include/PractRand/config.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/endian.h include/PractRand/rng_internals.h
src/RNGs/xsm.o: src/RNGs/xsm.cpp include/PractRand/config.h \
  include/PractRand/rng_basics.h include/PractRand/rng_helpers.h \
  include/PractRand/rng_adaptors.h include/PractRand/rng_internals.h \
  include/PractRand/RNGs/xsm32.h include/PractRand/RNGs/xsm64.h
tools/RNG_test.o: tools/RNG_test.cpp include/PractRand_full.h \
  include/PractRand/config.h include/PractRand/rng_basics.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/test_batteries.h include/PractRand/test_helpers.h \
  include/PractRand/tests.h include/PractRand/RNGs/all.h \
  include/PractRand/RNGs/arbee.h include/PractRand/RNGs/chacha.h \
  include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/RNGs/hc256.h include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/RNGs/jsf32.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/RNGs/mt19937.h \
  include/PractRand/RNGs/rarns16.h include/PractRand/RNGs/rarns32.h \
  include/PractRand/RNGs/rarns64.h include/PractRand/RNGs/salsa.h \
  include/PractRand/RNGs/sfc16.h include/PractRand/RNGs/sfc32.h \
  include/PractRand/RNGs/sfc64.h \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/RNGs/trivium.h include/PractRand/RNGs/xsm32.h \
  include/PractRand/RNGs/xsm64.h \
  include/PractRand/RNGs/other/fibonacci.h \
  include/PractRand/RNGs/other/indirection.h \
  include/PractRand/RNGs/other/mult.h \
  include/PractRand/RNGs/other/simple.h \
  include/PractRand/RNGs/other/special.h \
  include/PractRand/RNGs/other/transform.h \
  include/PractRand/rng_internals.h tools/RNG_from_name.h \
  tools/TestManager.h tools/MultithreadedTestManager.h \
  tools/multithreading.h tools/Candidate_RNGs.h tools/SeedingTester.h \
  include/PractRand/Tests/Birthday.h include/PractRand/Tests/DistFreq4.h \
  include/PractRand/Tests/FPF.h include/PractRand/Tests/FPMulti.h \
  include/PractRand/Tests/Gap16.h
tools/RNG_benchmark.o: tools/RNG_benchmark.cpp include/PractRand.h \
  include/PractRand/config.h include/PractRand/rng_basics.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/RNGs/all.h include/PractRand/RNGs/arbee.h \
  include/PractRand/RNGs/chacha.h include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/RNGs/hc256.h include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/RNGs/jsf32.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/RNGs/mt19937.h \
  include/PractRand/RNGs/rarns16.h include/PractRand/RNGs/rarns32.h \
  include/PractRand/RNGs/rarns64.h include/PractRand/RNGs/salsa.h \
  include/PractRand/RNGs/sfc16.h include/PractRand/RNGs/sfc32.h \
  include/PractRand/RNGs/sfc64.h \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/RNGs/trivium.h include/PractRand/RNGs/xsm32.h \
  include/PractRand/RNGs/xsm64.h tools/Candidate_RNGs.h \
  tools/RNG_from_name.h include/PractRand/RNGs/other/fibonacci.h \
  include/PractRand/RNGs/other/indirection.h \
  include/PractRand/RNGs/other/mult.h \
  include/PractRand/RNGs/other/simple.h \
  include/PractRand/RNGs/other/transform.h \
  tools/measure_RNG_performance.h
tools/RNG_output.o: tools/RNG_output.cpp include/PractRand_full.h \
  include/PractRand/config.h include/PractRand/rng_basics.h \
  include/PractRand/rng_helpers.h include/PractRand/rng_adaptors.h \
  include/PractRand/test_batteries.h include/PractRand/test_helpers.h \
  include/PractRand/tests.h include/PractRand/RNGs/all.h \
  include/PractRand/RNGs/arbee.h include/PractRand/RNGs/chacha.h \
  include/PractRand/RNGs/efiix16x48.h \
  include/PractRand/RNGs/efiix32x48.h \
  include/PractRand/RNGs/efiix64x48.h include/PractRand/RNGs/efiix8x48.h \
  include/PractRand/RNGs/hc256.h include/PractRand/RNGs/isaac32x256.h \
  include/PractRand/RNGs/isaac64x256.h include/PractRand/RNGs/jsf32.h \
  include/PractRand/RNGs/jsf64.h include/PractRand/RNGs/mt19937.h \
  include/PractRand/RNGs/rarns16.h include/PractRand/RNGs/rarns32.h \
  include/PractRand/RNGs/rarns64.h include/PractRand/RNGs/salsa.h \
  include/PractRand/RNGs/sfc16.h include/PractRand/RNGs/sfc32.h \
  include/PractRand/RNGs/sfc64.h \
  include/PractRand/RNGs/sha2_based_pool.h \
  include/PractRand/RNGs/trivium.h include/PractRand/RNGs/xsm32.h \
  include/PractRand/RNGs/xsm64.h \
  include/PractRand/RNGs/other/fibonacci.h \
  include/PractRand/RNGs/other/indirection.h \
  include/PractRand/RNGs/other/mult.h \
  include/PractRand/RNGs/other/simple.h \
  include/PractRand/RNGs/other/special.h \
  include/PractRand/RNGs/other/transform.h tools/RNG_from_name.h \
  tools/Candidate_RNGs.h tools/SeedingTester.h

obj =  \
 src/math.o \
 src/non_uniform.o \
 src/platform_specifics.o \
 src/rand.o \
 src/sha2.o \
 src/test_batteries.o \
 src/tests.o \
 src/RNGs/arbee.o \
 src/RNGs/chacha.o \
 src/RNGs/efiix.o \
 src/RNGs/hc256.o \
 src/RNGs/isaac.o \
 src/RNGs/jsf.o \
 src/RNGs/mt19937.o \
 src/RNGs/rarns.o \
 src/RNGs/salsa.o \
 src/RNGs/sfc.o \
 src/RNGs/sha2_based_pool.o \
 src/RNGs/trivium.o \
 src/RNGs/xsm.o \
 src/RNGs/other/fibonacci.o \
 src/RNGs/other/indirection.o \
 src/RNGs/other/mult.o \
 src/RNGs/other/rng_sets.o \
 src/RNGs/other/simple.o \
 src/RNGs/other/transform.o

RNG_test: $(obj)  \
 src/RNGs/arbee.o \
 src/RNGs/chacha.o \
 src/RNGs/efiix.o \
 src/RNGs/hc256.o \
 src/RNGs/isaac.o \
 src/RNGs/jsf.o \
 src/RNGs/mt19937.o \
 src/RNGs/rarns.o \
 src/RNGs/salsa.o \
 src/RNGs/sfc.o \
 src/RNGs/sha2_based_pool.o \
 src/RNGs/trivium.o \
 src/RNGs/xsm.o tools/RNG_test.o
	$(CXX) $(LDFLAGS) $^ -o $@

RNG_benchmark: $(obj)  \
 src/RNGs/arbee.o \
 src/RNGs/chacha.o \
 src/RNGs/efiix.o \
 src/RNGs/hc256.o \
 src/RNGs/isaac.o \
 src/RNGs/jsf.o \
 src/RNGs/mt19937.o \
 src/RNGs/rarns.o \
 src/RNGs/salsa.o \
 src/RNGs/sfc.o \
 src/RNGs/sha2_based_pool.o \
 src/RNGs/trivium.o \
 src/RNGs/xsm.o tools/RNG_benchmark.o
	$(CXX) $(LDFLAGS) $^ -o $@

RNG_output: $(obj)  \
 src/RNGs/arbee.o \
 src/RNGs/chacha.o \
 src/RNGs/efiix.o \
 src/RNGs/hc256.o \
 src/RNGs/isaac.o \
 src/RNGs/jsf.o \
 src/RNGs/mt19937.o \
 src/RNGs/rarns.o \
 src/RNGs/salsa.o \
 src/RNGs/sfc.o \
 src/RNGs/sha2_based_pool.o \
 src/RNGs/trivium.o \
 src/RNGs/xsm.o tools/RNG_output.o
	$(CXX) $(LDFLAGS) $^ -o $@

$(STATIC_LIB): $(obj)
	$(AR) rcs $@ $^

$(SHARED_LIB): $(obj)
	$(CXX) $(SHLIB_LDFLAGS) -o $@ $^

clean:
	rm -f RNG_test $(obj)  \
 src/RNGs/arbee.o \
 src/RNGs/chacha.o \
 src/RNGs/efiix.o \
 src/RNGs/hc256.o \
 src/RNGs/isaac.o \
 src/RNGs/jsf.o \
 src/RNGs/mt19937.o \
 src/RNGs/rarns.o \
 src/RNGs/salsa.o \
 src/RNGs/sfc.o \
 src/RNGs/sha2_based_pool.o \
 src/RNGs/trivium.o \
 src/RNGs/xsm.o 

.PHONY: all clean

