#pragma once

#include <cstdint>

//#include <vector>

#define PRACTRAND_RANDF_IMPLEMENTATION(RNG)  {return  ((RNG).raw32() & ((uint32_t(1) << 24)-1)) *  float(1.0/16777216.0);}
#define PRACTRAND_RANDLF_IMPLEMENTATION(RNG)  {return ((RNG).raw64() & ((uint64_t(1) << 53)-1)) * (1.0/9007199254740992.0);}

#define PRACTRAND_RANDI_IMPLEMENTATION(max)     \
	uint32_t mask, tmp;\
	(max) -= 1;\
	mask = max;\
	mask |= mask >> 1; mask |= mask >>  2; mask |= mask >> 4;\
	mask |= mask >> 8; mask |= mask >> 16;\
	do {\
		tmp = raw32() & mask;\
	} while (tmp > (max));\
	return tmp;
#define PRACTRAND_RANDLI_IMPLEMENTATION(max)     \
	uint64_t mask, tmp;\
	(max) -= 1;\
	mask = max;\
	mask |= mask >> 1; mask |= mask >>  2; mask |= mask >>  4;\
	mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;\
	do {\
		tmp = raw64() & mask;\
	} while (tmp > (max));\
	return tmp;

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PRACTRAND_POLYMORPHIC_RNG_BASICS_C8(RNG) \
	uint8_t  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return implementation.raw8();}\
	uint16_t PractRand::RNGs::Polymorphic:: RNG ::raw16() {uint16_t r = implementation.raw8(); return (uint16_t(implementation.raw8()) << 8) | r;}\
	uint32_t PractRand::RNGs::Polymorphic:: RNG ::raw32() {\
		uint32_t r = implementation.raw8();\
		r = r | (uint32_t(implementation.raw8()) << 8);\
		r = r | (uint32_t(implementation.raw8()) << 16);\
		return r | (uint32_t(implementation.raw8()) << 24);\
	}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::raw64() {\
		uint64_t r = raw32();\
		return r | (uint64_t(raw32()) << 32);\
	}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject* walker) {\
		implementation.walk_state(walker);\
	}
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PRACTRAND_POLYMORPHIC_RNG_BASICS_C16(RNG) \
	uint8_t  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return uint8_t(implementation.raw16());}\
	uint16_t PractRand::RNGs::Polymorphic:: RNG ::raw16() {return implementation.raw16();}\
	uint32_t PractRand::RNGs::Polymorphic:: RNG ::raw32() {\
		uint16_t r = implementation.raw16();\
		return r | (uint32_t(implementation.raw16()) << 16);\
	}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::raw64() {\
		uint64_t r = implementation.raw16();\
		r = r | (uint32_t(implementation.raw16()) << 16);\
		r = r | (uint64_t(implementation.raw16()) << 32);\
		return r | (uint64_t(implementation.raw16()) << 48);\
	}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject* walker) {\
		implementation.walk_state(walker);\
	}
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PRACTRAND_POLYMORPHIC_RNG_BASICS_C32(RNG) \
	uint8_t  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return uint8_t (implementation.raw32());}\
	uint16_t PractRand::RNGs::Polymorphic:: RNG ::raw16() {return uint16_t(implementation.raw32());}\
	uint32_t PractRand::RNGs::Polymorphic:: RNG ::raw32() {return implementation.raw32();}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::raw64() {\
		uint64_t r = implementation.raw32();\
		return (r << 32) | implementation.raw32();\
	}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject* walker) {\
		implementation.walk_state(walker);\
	}
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PRACTRAND_POLYMORPHIC_RNG_BASICS_C64(RNG) \
	uint8_t  PractRand::RNGs::Polymorphic:: RNG ::raw8 () {return uint8_t (implementation.raw64());}\
	uint16_t PractRand::RNGs::Polymorphic:: RNG ::raw16() {return uint16_t(implementation.raw64());}\
	uint32_t PractRand::RNGs::Polymorphic:: RNG ::raw32() {return uint32_t(implementation.raw64());}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::raw64() {return implementation.raw64();}\
	uint64_t PractRand::RNGs::Polymorphic:: RNG ::get_flags() const {return implementation.FLAGS;}\
	void PractRand::RNGs::Polymorphic:: RNG ::walk_state(StateWalkingObject* walker) {\
		implementation.walk_state(walker);\
	}



namespace PractRand {
	namespace RNGs {
		class vRNG;
	}
	namespace Internals {
		//inline
		//non_uniform.cpp
		double generate_gaussian_fast(uint64_t raw64);//fast CDF-based hybrid method
		//double generate_gaussian_high_quality(uint64_t raw192[3]);//slow PDF-based ziggurat method
		//rand.cpp
		void test_random_access(PractRand::RNGs::vRNG* rng, PractRand::RNGs::vRNG* known_good, uint64_t period_low64, uint64_t period_high64);
		//platform_specific.cpp
		bool add_entropy_automatically( PractRand::RNGs::vRNG* entropy_pool, int milliseconds=0 );
		uint64_t issue_unique_identifier();
		uint64_t high_resolution_time();
		//math.cpp
		void fast_forward_lcg128 ( uint64_t how_far_low, uint64_t how_far_high, uint64_t& value_low, uint64_t& value_high, uint64_t mul_low, uint64_t mul_high, uint64_t add_low, uint64_t add_high );
		uint64_t fast_forward_lcg64 ( uint64_t how_far, uint64_t val, uint64_t mul, uint64_t add );
		uint32_t fast_forward_lcg32 ( uint32_t how_far, uint32_t val, uint32_t mul, uint32_t add );
		uint32_t fast_forward_lcg32c ( uint32_t how_far, uint32_t val, uint32_t mul, uint32_t add, uint32_t mod );
		/*class XorshiftMatrix {
			//matrix representing a state transition function for an RNG that uses only shifts and xors
			std::vector<bool> bits;
			int size;
		public:
			XorshiftMatrix(int size_, bool identity);
			void apply(const std::vector<bool> &input, std::vector<bool> &output);
			XorshiftMatrix operator*(const XorshiftMatrix &other) const;
			bool operator==(const XorshiftMatrix &other) const;
			bool operator!=(const XorshiftMatrix &other) const {return !(*this == other);}
			XorshiftMatrix exponent(uint64_t exponent_value) const;
			XorshiftMatrix exponent2Xminus1(uint64_t X) const;
			bool verify_period_factorization(const std::vector<uint64_t> &factors) const;
			bool get(int in_index, int out_index) const {return bits[in_index + out_index * size];}
			void set(int in_index, int out_index, bool value) {bits[in_index + out_index * size] = value;}
			void toggle(int in_index, int out_index, bool value) {bits[in_index + out_index * size] = !bits[in_index + out_index * size];}
		};*/
	}
}
