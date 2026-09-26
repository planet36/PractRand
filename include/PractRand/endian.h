#pragma once

#include "config.h"
#include <bit>

static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big,
	"PractRand requires a little-endian or big-endian target");

namespace PractRand {
	static inline Uint16 invert_endianness16(Uint16 v) {return (v >> 8) | (v << 8);}
	static inline Uint32 invert_endianness32(Uint32 v) {
		v = ((v & 0xFF00FF00) >> 8) | ((v & 0x00ff00ff) << 8);
		return (v >> 16) | (v << 16);
	}
	static inline Uint64 invert_endianness64(Uint64 v) {
		v = ((v & 0xFF00FF00FF00FF00ULL) >> 8) | ((v & 0x00ff00ff00ff00ffULL) << 8);
		v = ((v & 0xFFFF0000FFFF0000ULL) >> 16) | ((v & 0x0000ffff0000ffffULL) << 16);
		return (v >> 32) | (v << 32);
	}
#if 0
	static inline Uint16 little_endian_conversion16 ( Uint16 v ) {
		if constexpr (std::endian::native == std::endian::little) return v;
		else return invert_endianness16(v);
	}
	static inline Uint32 little_endian_conversion32 ( Uint32 v ) {
		if constexpr (std::endian::native == std::endian::little) return v;
		else return invert_endianness32(v);
	}
#endif
	static inline Uint64 little_endian_conversion64 ( Uint64 v ) {
		if constexpr (std::endian::native == std::endian::little) return v;
		else return invert_endianness64(v);
	}
#if 0
	static inline Uint16 big_endian_conversion16 ( Uint16 v ) {
		if constexpr (std::endian::native == std::endian::big) return v;
		else return invert_endianness16(v);
	}
	static inline Uint32 big_endian_conversion32 ( Uint32 v ) {
		if constexpr (std::endian::native == std::endian::big) return v;
		else return invert_endianness32(v);
	}
	static inline Uint64 big_endian_conversion64 ( Uint64 v ) {
		if constexpr (std::endian::native == std::endian::big) return v;
		else return invert_endianness64(v);
	}
#endif
#if 0
#if defined PRACTRAND_TARGET_IS_LITTLE_ENDIAN
	union split_int_16 {
		Uint16 whole;
		struct blah {
			Uint8 low8;
			Uint8 high8;
		} split;
	};
	union split_int_32 {
		Uint32 whole;
		struct blah {
			split_int_16 low16;
			split_int_16 high16;
		} split;
	};
	union split_int_64 {
		Uint64 whole;
		struct blah {
			split_int_32 low32;
			split_int_32 high32;
		} split;
	};
#elif defined PRACTRAND_TARGET_IS_BIG_ENDIAN
	union split_int_16 {
		Uint16 whole;
		struct blah {
			Uint8 high8;
			Uint8 low8;
		} split;
	};
	union split_int_32 {
		Uint32 whole;
		struct blah {
			split_int_16 high16;
			split_int_16 low16;
		} split;
	};
	union split_int_64 {
		Uint64 whole;
		struct blah {
			split_int_32 high32;
			split_int_32 low32;
		} split;
	};
#endif
#endif
}
