#ifndef _practrand_rng_helpers_h
#define _practrand_rng_helpers_h

#ifndef __PRACTRAND_CONFIG_H__
#include "PractRand/config.h"
#endif //__PRACTRAND_CONFIG_H__

namespace PractRand {
	namespace Internals {
		//inline
#if defined _MSC_VER && !defined PRACTRAND_NO_INTRINSICS
		static inline Uint8  rotate8(Uint8  value, int bits) { return _rotl8(value, bits); }
		static inline Uint16 rotate16(Uint16 value, int bits) { return _rotl16(value, bits); }
		static inline Uint32 rotate32(Uint32 value, int bits) { return _rotl(value, bits); }
		static inline Uint64 rotate64(Uint64 value, int bits) { return _rotl64(value, bits); }
#else
		// compilers are smart enough to recognize this for fixed shifts, but reports suggest they tend to add unnecessary bitmasking for variable shifts
		static inline Uint8  rotate8(Uint8  value, int bits) { return (value << bits) | (value >> (8 - bits)); }
		static inline Uint16 rotate16(Uint16 value, int bits) { return (value << bits) | (value >> (16 - bits)); }
		static inline Uint32 rotate32(Uint32 value, int bits) { return (value << bits) | (value >> (32 - bits)); }
		static inline Uint64 rotate64(Uint64 value, int bits) { return (value << bits) | (value >> (64 - bits)); }
#endif
		static inline Uint8  rotate(Uint8  value, int bits) { return rotate8(value, bits); }
		static inline Uint16 rotate(Uint16 value, int bits) { return rotate16(value, bits); }
		static inline Uint32 rotate(Uint32 value, int bits) { return rotate32(value, bits); }
		static inline Uint64 rotate(Uint64 value, int bits) { return rotate64(value, bits); }
		int count_ones8(Uint8 a);// now known as popcount, change name, add alias?
		int count_ones16(Uint16 a);
		int count_ones32(Uint32 a);
		int count_ones64(Uint64 a);
		inline int count_ones_(Uint8  a) { return count_ones8(a); }
		inline int count_ones_(Uint16 a) { return count_ones16(a); }
		inline int count_ones_(Uint32 a) { return count_ones32(a); }
		inline int count_ones_(Uint64 a) { return count_ones64(a); }
		Uint8  reverse_bits8(Uint8);
		Uint16 reverse_bits16(Uint16);
		Uint32 reverse_bits32(Uint32);
		Uint64 reverse_bits64(Uint64);
		unsigned long count_low_zeroes32(Uint32);
		unsigned long count_low_zeroes64(Uint64);
		unsigned long count_high_zeroes32(Uint32);
		unsigned long count_high_zeroes64(Uint64);
		long ilog2_32(Uint32);//returns -1 on zero
		long ilog2_64(Uint64);//returns -1 on zero
	}
	class StateWalkingObject {
	public:
		virtual ~StateWalkingObject() {}
		virtual void handle(bool  &) = 0;
		/*virtual void handle(unsigned char &) = 0;
		virtual void handle(unsigned short &) = 0;
		virtual void handle(unsigned int &) = 0;
		virtual void handle(unsigned long &) = 0;
		virtual void handle(unsigned long long &) = 0;*/
		virtual void handle(Uint8 &) = 0;
		virtual void handle(Uint16&) = 0;
		virtual void handle(Uint32&) = 0;
		virtual void handle(Uint64&) = 0;

		virtual void handle(float &) = 0;
		virtual void handle(double &) = 0;

		//purposes: 
		// 1. serialization (serialize, deserialize, measure state size)
		// 2. seeding of last resort (when an algorithm provides no seeding function, a generic one using this lets it work anyway)
		// 3. avalanche style testing (measure hamming distance between two states, generate a state at a fixed hamming distances from another state, etc)
		// 4. ???
		enum {
			FLAG_READ_ONLY = 1, // does not make changes
			FLAG_WRITE_ONLY = 2,// result does not depend upon prior state
			FLAG_CLUMSY = 4,    // may violate invariants (if also FLAG_READ_ONLY then only wants to see state visible to clumsy writers)
			FLAG_SEEDER = 8,    // some RNGs may have extra invariants enforced only on seeded states (minimum distance away from other seeded states on cycle)
		};
		virtual Uint32 get_properties() const = 0;
		bool is_read_only() const { return (get_properties() & FLAG_READ_ONLY) ? true : false; }
		bool is_write_only() const { return (get_properties() & FLAG_WRITE_ONLY) ? true : false; }
		bool is_clumsy() const { return (get_properties() & FLAG_CLUMSY) ? true : false; }
		bool is_seeder() const {return (get_properties() & FLAG_SEEDER) ? true : false;}

		/*void handle(signed char      &v) {handle((unsigned char)v);}
		void handle(signed short     &v) {handle((unsigned short)v);}
		void handle(signed int       &v) {handle((unsigned int)v);}
		void handle(signed long      &v) {handle((unsigned long)v);}
		void handle(signed long long &v) {handle((unsigned long long)v);}*/
		void handle(Sint8 &v) {handle((Uint8 &)v);}
		void handle(Sint16&v) {handle((Uint16&)v);}
		void handle(Sint32&v) {handle((Uint32&)v);}
		void handle(Sint64&v) {handle((Uint64&)v);}

		StateWalkingObject &operator<<(Uint8 &v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Uint16&v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Uint32&v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Uint64&v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Sint8 &v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Sint16&v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Sint32&v) {handle(v);return *this;}
		StateWalkingObject &operator<<(Sint64&v) { handle(v); return *this; }
		StateWalkingObject &operator<<(float  &v){ handle(v); return *this; }
		StateWalkingObject &operator<<(double &v){handle(v);return *this;}
#if defined __SIZEOF_INT128__
		void handle(unsigned __int128 &v) {
			Uint64 *ptr = reinterpret_cast<Uint64*>(&v);
			// first low, then high
#if defined PRACTRAND_TARGET_IS_LITTLE_ENDIAN
			handle(ptr[0]);
			handle(ptr[1]);
#elif defined PRACTRAND_TARGET_IS_BIG_ENDIAN
			handle(ptr[1]);
			handle(ptr[0]);
#else
#error unknown endianness
#endif
		}
		void handle(signed __int128 &v) { handle((unsigned __int128&)v); }
		StateWalkingObject &operator<<(unsigned __int128 &v) { handle(v); return *this; }
		StateWalkingObject &operator<<(signed __int128 &v) { handle(v); return *this; }
#endif
	};
	namespace RNGs {
		class vRNG;
	}
	namespace Internals {
		Uint32 rand_i32_fast_implementation(Uint32 random_value, Uint32 max);//attempts to map 32 random bits to (approximately) uniform integers in the range [0, max) -- note that this is impossible to do perfectly, especially for large values of max
		//Uint64 randi64_fast_implementation(Uint64 random_value, Uint64 max);//attempts to map 64 random bits to (approximately) uniform integers in the range [0, max) -- note that this is impossible to do perfectly, especially for large values of max
		float rand_float_implementation(Uint32 random_bits);//maps 32 random bits to a random ieee754 single-precision floating point number ; output is in the range [0,1)
		double rand_double_implementation(Uint64 random_bits);//maps 64 random bits to a random ieee754 double-precision floating point number ; output is in the range [0,1)

		StateWalkingObject *int_to_rng_seeder(Uint64 seed_low, Uint64 seed_high=0);//must be deleted after use
		StateWalkingObject *vrng_to_rng_seeder(RNGs::vRNG *);//must be deleted after use
		StateWalkingObject *get_autoseeder(const void *);//must be deleted after use
		//non_uniform.cpp
		// I'm going to get rid of some of these after I have more time to test them out
		// each converts raw random bits to gaussian distribution
		double generate_gaussian_lookup3_linear(Uint64 raw64); // 277 M/sec on test computer
		double generate_gaussian_popcount_lookup_linear(Uint64, Uint64);// 248 M/sec on test computer
		double generate_gaussian_popcount_lookup_quadratic(Uint64, Uint64);// 233 M/sec on test computer
		double generate_gaussian_popcount_lookup_cubic(Uint64, Uint64);// 222 M/sec on test computer
		double generate_gaussian_popcount_sum(Uint64, Uint64);// 228 M/sec on test computer
		double generate_gaussian_popcount_sum2(Uint64, Uint64, Uint64, Uint64, Uint64);// 114 M/sec on test computer
		double generate_gaussian_high_quality(Uint64 a, Uint64 b, Uint64 c);//
		bool _generate_gaussian_ziggurat(double &result, Uint64, Uint64);// ziggurat-based implementation needs slightly different interface
		template<typename RNG> double generate_with_ziggurat(RNG& rng) { // 173 M/sec on test computer
			while (true) {
				double rv;  Uint64 a = rng.raw64(); Uint64 b = rng.raw64();
				if (PractRand::Internals::_generate_gaussian_ziggurat(rv, a, b)) return rv;
			}
		}

		//this is used to figure out the shift constants for specific integer sizes, it doesn't need to support non-powers-of-2, nor very large values
		template <unsigned int N> struct StaticLog2 {};
		template <> struct StaticLog2<1> { enum { value = 0 }; };
		template <> struct StaticLog2<2> { enum { value = 1 }; };
		template <> struct StaticLog2<4> { enum { value = 2 }; };
		template <> struct StaticLog2<8> { enum { value = 3 }; };
		template <> struct StaticLog2<16> { enum { value = 4 }; };
		template <> struct StaticLog2<32> { enum { value = 5 }; };
		template <> struct StaticLog2<64> { enum { value = 6 }; };
		template <> struct StaticLog2<128> { enum { value = 7 }; };
		template <> struct StaticLog2<256> { enum { value = 8 }; };
		template <> struct StaticLog2<512> { enum { value = 9 }; };
		template <> struct StaticLog2<1024> { enum { value = 10 }; };
		template <> struct StaticLog2<2048> { enum { value = 11 }; };
		template <> struct StaticLog2<4096> { enum { value = 12 }; };
		template <> struct StaticLog2<8192> { enum { value = 13 }; };
		template <> struct StaticLog2<16384> { enum { value = 14 }; };
		template <> struct StaticLog2<32768> { enum { value = 15 }; };
		template <> struct StaticLog2<65536> { enum { value = 16 }; };

		//rand.cpp
		void test_random_access(PractRand::RNGs::vRNG* rng, PractRand::RNGs::vRNG* known_good, Uint64 period_low64, Uint64 period_high64);
		//platform_specific.cpp
		bool add_entropy_automatically(PractRand::RNGs::vRNG* entropy_pool, int milliseconds = 0);
		Uint64 issue_unique_identifier();
		Sint64 high_resolution_time();
		//math.cpp
		void fast_forward_lcg128(Uint64 how_far_low, Uint64 how_far_high, Uint64& value_low, Uint64& value_high, Uint64 mul_low, Uint64 mul_high, Uint64 add_low, Uint64 add_high);
		Uint64 fast_forward_lcg64(Uint64 how_far, Uint64 val, Uint64 mul, Uint64 add);
		Uint32 fast_forward_lcg32(Uint32 how_far, Uint32 val, Uint32 mul, Uint32 add);
		Uint32 fast_forward_lcg32c(Uint32 how_far, Uint32 val, Uint32 mul, Uint32 add, Uint32 mod);
	}
}

#define PRACTRAND__RANDI32_IMPLEMENTATION(max)     \
	Uint32 mask, tmp;\
	max -= 1;\
	mask = max;\
	mask |= mask >> 1; mask |= mask >>  2; mask |= mask >> 4;\
	mask |= mask >> 8; mask |= mask >> 16;\
	do {\
		tmp = raw32() & mask;\
	} while (tmp > max);\
	return tmp;
#define PRACTRAND__RANDI64_IMPLEMENTATION(max)     \
	Uint64 mask, tmp;\
	max -= 1;\
	mask = max;\
	mask |= mask >> 1; mask |= mask >>  2; mask |= mask >>  4;\
	mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;\
	do {\
		tmp = raw64() & mask;\
	} while (tmp > max);\
	return tmp;

#define PRACTRAND__POLYMORPHIC_RNG_BASICS_C8(RNG) \
	Uint8  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return implementation.raw8();}\
	Uint16 PractRand::RNGs::Polymorphic:: RNG ::raw16() {Uint16 r = implementation.raw8(); return (Uint16(implementation.raw8()) << 8) | r;}\
	Uint32 PractRand::RNGs::Polymorphic:: RNG ::raw32() {\
		Uint32 r = implementation.raw8();\
		r = r | (Uint32(implementation.raw8()) << 8);\
		r = r | (Uint32(implementation.raw8()) << 16);\
		return r | (Uint32(implementation.raw8()) << 24);\
	}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::raw64() {\
		Uint64 r = raw32();\
		return r | (Uint64(raw32()) << 32);\
	}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject *walker) {\
		implementation.walk_state(walker);\
	}
#define PRACTRAND__POLYMORPHIC_RNG_BASICS_C16(RNG) \
	Uint8  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return Uint8(implementation.raw16());}\
	Uint16 PractRand::RNGs::Polymorphic:: RNG ::raw16() {return implementation.raw16();}\
	Uint32 PractRand::RNGs::Polymorphic:: RNG ::raw32() {\
		Uint16 r = implementation.raw16();\
		return r | (Uint32(implementation.raw16()) << 16);\
	}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::raw64() {\
		Uint64 r = implementation.raw16();\
		r = r | (Uint32(implementation.raw16()) << 16);\
		r = r | (Uint64(implementation.raw16()) << 32);\
		return r | (Uint64(implementation.raw16()) << 48);\
	}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject *walker) {\
		implementation.walk_state(walker);\
	}
#define PRACTRAND__POLYMORPHIC_RNG_BASICS_C32(RNG) \
	Uint8  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return Uint8 (implementation.raw32());}\
	Uint16 PractRand::RNGs::Polymorphic:: RNG ::raw16() {return Uint16(implementation.raw32());}\
	Uint32 PractRand::RNGs::Polymorphic:: RNG ::raw32() {return Uint32(implementation.raw32());}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::raw64() {\
		Uint64 r = implementation.raw32();\
		return (r << 32) | implementation.raw32();\
	}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject *walker) {\
		implementation.walk_state(walker);\
	}
#define PRACTRAND__POLYMORPHIC_RNG_BASICS_C64(RNG) \
	Uint8  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return Uint8 (implementation.raw64());}\
	Uint16 PractRand::RNGs::Polymorphic:: RNG ::raw16() {return Uint16(implementation.raw64());}\
	Uint32 PractRand::RNGs::Polymorphic:: RNG ::raw32() {return Uint32(implementation.raw64());}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::raw64() {return implementation.raw64();}\
	Uint64 PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject *walker) {\
		implementation.walk_state(walker);\
	}
#define PRACTRAND__POLYMORPHIC_RNG_BASICS_H(RNG) public:\
		enum {OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL,OUTPUT_BITS = Raw:: RNG ::OUTPUT_BITS,FLAGS = Raw:: RNG ::FLAGS};\
		Raw:: RNG implementation;\
		RNG (Uint64 seed_low, Uint64 seed_high=0) {seed(seed_low, seed_high);}\
		RNG (vRNG *seeder) {seed(seeder);}\
		RNG (SEED_AUTO_TYPE ) {autoseed();}\
		RNG (SEED_NONE_TYPE ) {}\
		Uint8  raw8 ();\
		Uint16 raw16();\
		Uint32 raw32();\
		Uint64 raw64();\
		using vRNG::seed;\
		Uint64 get_flags() const;\
		std::string get_name() const;\
		void walk_state(StateWalkingObject *walker);

#if defined PRACTRAND_NO_LIGHT_WEIGHT_RNGS
#define PRACTRAND__LIGHT_WEIGHT_RNG(RNG)
#define PRACTRAND__LIGHT_WEIGHT_ENTROPY_POOL(RNG)
#else // ! PRACTRAND_NO_LIGHT_WEIGHT_RNGS
#ifndef __PRACTRAND_RNG_ADAPTORS_H__
#include "rng_adaptors.h"
#endif//__PRACTRAND_RNG_ADAPTORS_H__
#define PRACTRAND__LIGHT_WEIGHT_RNG(RNG) 	\
	namespace LightWeight {\
		typedef PractRand::RNGs::Adaptors::RAW_TO_LIGHT_WEIGHT_RNG<PractRand::RNGs::Raw:: RNG > RNG;\
	}
#endif//PRACTRAND_PROVIDE_LIGHT_WEIGHT_RNGS


#endif //_practrand_rng_helpers_h
