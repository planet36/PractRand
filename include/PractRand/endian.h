#pragma once

#include "config.h"

#include <bit>

static_assert(std::endian::native == std::endian::little, "PractRand requires a little-endian target");

namespace PractRand {
#if 0
	static inline uint16_t invert_endianness16(uint16_t v) {return (v >> 8) | (v << 8);}
	static inline uint32_t invert_endianness32(uint32_t v) {
		v = ((v & 0xFF00FF00) >> 8) | ((v & 0x00ff00ff) << 8);
		return (v >> 16) | (v << 16);
	}
	static inline uint64_t invert_endianness64(uint64_t v) {
		v = ((v & 0xFF00FF00FF00FF00ULL) >> 8) | ((v & 0x00ff00ff00ff00ffULL) << 8);
		v = ((v & 0xFFFF0000FFFF0000ULL) >> 16) | ((v & 0x0000ffff0000ffffULL) << 16);
		return (v >> 32) | (v << 32);
	}
#endif
#if 0
	union split_int_16 {
		uint16_t whole;
		struct blah {
			uint8_t low8;
			uint8_t high8;
		} split;
	};
	union split_int_32 {
		uint32_t whole;
		struct blah {
			split_int_16 low16;
			split_int_16 high16;
		} split;
	};
	union split_int_64 {
		uint64_t whole;
		struct blah {
			split_int_32 low32;
			split_int_32 high32;
		} split;
	};
#endif
}
