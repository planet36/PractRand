# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

export LC_ALL := C

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

LIB_SRCS = $(wildcard src/*.cpp src/RNGs/*.cpp src/RNGs/other/*.cpp)
LIB_DEPS = $(LIB_SRCS:.cpp=.d)
LIB_OBJS = $(LIB_SRCS:.cpp=.o)
LIB = libPractRand.a
BINS = RNG_benchmark RNG_output RNG_test

CPPFLAGS += -MMD -MP
CPPFLAGS += -Iinclude

CXXFLAGS += -O3 -flto=auto -march=native

LDLIBS += -pthread

all: $(BINS)

$(LIB): $(LIB_OBJS)
	ar rcs $@ $^

# https://www.gnu.org/software/make/manual/html_node/Static-Usage.html
# Static Pattern Rule
$(BINS): RNG_% : tools/RNG_%.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIB)

install: $(BINS)
	@mkdir -v -p -- $(DESTDIR)$(BINDIR)
	@cp -v -f -- $(BINS) $(DESTDIR)$(BINDIR)

clean:
	@$(RM) --verbose -- $(LIB_DEPS) $(LIB_OBJS) $(LIB) $(BINS) $(addsuffix .d,$(BINS))

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: all clean install

-include $(LIB_DEPS)
