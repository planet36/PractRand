# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# https://pracrand.sourceforge.net/installation.txt
# https://www.pcg-random.org/posts/how-to-test-with-practrand.html
# https://espadrine.github.io/blog/posts/a-primer-on-randomness.html

export LC_ALL := C

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

LIB_SRCS = $(wildcard src/*.cpp src/RNGs/*.cpp src/RNGs/other/*.cpp)
LIB_DEPS = $(LIB_SRCS:.cpp=.d)
LIB_OBJS = $(LIB_SRCS:.cpp=.o)
LIB = libPractRand.a

BIN_SRCS = $(wildcard tools/RNG_*.cpp)
BIN_DEPS = $(BIN_SRCS:.cpp=.d)
BINS = $(basename $(BIN_SRCS))

CPPFLAGS = -MMD -MP
CPPFLAGS += -Iinclude

CXXFLAGS = -std=c++26
CXXFLAGS += -pipe -Wall -Wextra -Wpedantic -Wfatal-errors
CXXFLAGS += -O3 -flto=auto -march=native -fno-math-errno

LDLIBS += -pthread

all: $(BINS)

$(LIB): $(LIB_OBJS)
	ar rcs $@ $^

# The built-in recipe for the implicit rule uses $^ instead of $<
# https://www.gnu.org/software/make/manual/html_node/Static-Usage.html
# Static Pattern Rule
$(BINS): tools/RNG_% : tools/RNG_%.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIB)

install: $(BINS) | $(DESTDIR)$(BINDIR)
	@cp -v -f -- $^ $(DESTDIR)$(BINDIR)

$(DESTDIR)$(BINDIR):
	@mkdir -v -p -- $(DESTDIR)$(BINDIR)

clean:
	@$(RM) --verbose -- $(LIB_DEPS) $(LIB_OBJS) $(LIB) $(BIN_DEPS) $(BINS)

lint:
	-clang-tidy --quiet $(LIB_SRCS) $(BIN_SRCS) -- $(CPPFLAGS) $(CXXFLAGS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: all clean install lint

-include $(LIB_DEPS) $(BIN_DEPS)
