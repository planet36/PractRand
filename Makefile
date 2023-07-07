# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

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

CPPFLAGS += -MMD -MP
CPPFLAGS += -Iinclude

CXXFLAGS += -O3 -flto=auto -march=native

LDLIBS += -pthread

all: $(BINS)

$(LIB): $(LIB_OBJS)
	ar rcs $@ $^

# https://www.gnu.org/software/make/manual/html_node/Static-Usage.html
# Static Pattern Rule
$(BINS): tools/RNG_% : tools/RNG_%.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIB)

install: $(BINS)
	@mkdir -v -p -- $(DESTDIR)$(BINDIR)
	@cp -v -f -- $(BINS) $(DESTDIR)$(BINDIR)

clean:
	@$(RM) --verbose -- $(LIB_DEPS) $(LIB_OBJS) $(LIB) $(BIN_DEPS) $(BINS)

lint:
	-clang-tidy --quiet $(LIB_SRCS) $(BIN_SRCS) -- $(CPPFLAGS) $(CXXFLAGS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: all clean install lint

-include $(LIB_DEPS) $(BIN_DEPS)
