#include "PractRand/RNGs/xsm32.h"
#include "PractRand/RNGs/xsm64.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <bit>
#include <cstdint>
#include <string>

using namespace PractRand;
using namespace PractRand::Internals;

//polymorphic:
PRACTRAND_POLYMORPHIC_RNG_BASICS_C32(xsm32)
void PractRand::RNGs::Polymorphic::xsm32::seed(uint64_t s) { implementation.seed(s); }
void PractRand::RNGs::Polymorphic::xsm32::seed(vRNG* seeder_rng) { implementation.seed(seeder_rng); }
void PractRand::RNGs::Polymorphic::xsm32::seek_forward128(uint64_t how_far_low64, [[maybe_unused]] uint64_t how_far_high64) { implementation.seek_forward(how_far_low64); }
void PractRand::RNGs::Polymorphic::xsm32::seek_backward128(uint64_t how_far_low64, [[maybe_unused]] uint64_t how_far_high64) {implementation.seek_backward(how_far_low64);}
std::string PractRand::RNGs::Polymorphic::xsm32::get_name() const {return "xsm32";}

PRACTRAND_POLYMORPHIC_RNG_BASICS_C64(xsm64)
void PractRand::RNGs::Polymorphic::xsm64::seed(uint64_t s) { implementation.seed(s); }
void PractRand::RNGs::Polymorphic::xsm64::seed(uint64_t seed_low, uint64_t seed_high) { implementation.seed(seed_low, seed_high); }
void PractRand::RNGs::Polymorphic::xsm64::seed(vRNG* seeder_rng) { implementation.seed(seeder_rng); }
void PractRand::RNGs::Polymorphic::xsm64::seek_forward128(uint64_t how_far_low64, uint64_t how_far_high64) { implementation.seek_forward(how_far_low64, how_far_high64); }
void PractRand::RNGs::Polymorphic::xsm64::seek_backward128(uint64_t how_far_low64, uint64_t how_far_high64) { implementation.seek_backward(how_far_low64, how_far_high64); }
std::string PractRand::RNGs::Polymorphic::xsm64::get_name() const { return "xsm64"; }

//raw:

uint32_t PractRand::RNGs::Raw::xsm32::raw32() {
	//const uint32_t K = 0xE99569B5;  //	11101001100101010110100110110101 473  205 - best value found so far
	//const uint32_t K =	0x299529B5;	//	00101001100101010010100110110101 470  202
	//const uint32_t K =	0x6595a395;	//	01100101100101011010001110010101 465  204 - original high quality version
	//const uint32_t K = 0xD251CF2D;  //	11010010010100011100111100101101 463  196
	//const uint32_t K = 0x08304E2D;  //	00001000001100000100111000101101 463  194
	//const uint32_t K =	0xCA566EB5;	//	11001010010101100110111010101101 458  201
	//const uint32_t K = 0xD250CE2D;  //	11010010010100001100111000101101 456  194
	//const uint32_t K =	0x21042185;	//	00100001000001000010000110000101 421  188 - now the middle quality version
	//const uint32_t K = 0x21052045;  //	00100001000001010010000001000101 373  161 - original reduced quality version
	//const uint32_t K =	0x01002045;  //	00000001000000000010000001000101 363  157 - lowest quality version tested

	uint32_t tmp = lcg_high ^ std::rotl(lcg_high + lcg_low, 9); tmp ^= std::rotl(tmp + lcg_adder_high, 19);
	tmp *= 0xD251CF2D;
	step_forwards();
	tmp = tmp ^ std::rotl(tmp + lcg_high, 16);
	tmp *= uint32_t(0x299529B5);
	tmp ^= tmp >> 16;
	return tmp;
}
void PractRand::RNGs::Raw::xsm32::step_backwards() {
	bool carry = lcg_low < lcg_adder_low;
	lcg_low -= lcg_adder_low;
	lcg_high -= lcg_low + lcg_adder_high + (carry ? 1 : 0);
}
void PractRand::RNGs::Raw::xsm32::step_forwards() {
	uint32_t tmp = lcg_low + lcg_adder_high;
	lcg_low += lcg_adder_low;
	lcg_high += tmp + ((lcg_low < lcg_adder_low) ? 1 : 0);
}
void PractRand::RNGs::Raw::xsm32::seed(uint64_t s) {
	//guaranteed that no two seeds are closer than 2**63 apart on the same cycle

	lcg_adder_low = uint32_t(s) | 1;
	lcg_adder_high = uint32_t(s >> 32);
	lcg_high = lcg_adder_high + uint32_t(s << 31);//the one bit of seed that did not effect which cycle it's on
	lcg_low = lcg_adder_low;
	step_forwards();
}
void PractRand::RNGs::Raw::xsm32::seed(vRNG* seeder_rng) {
	//guaranteed that no two *distinct* seeded states are within 2**48 of each other on the same cycle ; 2**79 distinct seeded states are possible

	lcg_adder_low = seeder_rng->raw32() | 1;
	lcg_adder_high = seeder_rng->raw32();
	lcg_high = lcg_adder_high + (seeder_rng->raw32() << 16);
	lcg_low = lcg_adder_low;
}
void PractRand::RNGs::Raw::xsm32::walk_state(StateWalkingObject* walker) {
	walker->handle(lcg_low);
	walker->handle(lcg_high);
	walker->handle(lcg_adder_low);
	walker->handle(lcg_adder_high);
	if (walker->is_clumsy() && !walker->is_read_only()) {
		lcg_adder_low |= 1;
	}
}
void PractRand::RNGs::Raw::xsm32::seek_forward (uint64_t how_far) {
	if (!how_far) return;
	uint64_t lcg_state = lcg_low | (uint64_t(lcg_high) << 32);
	uint64_t lcg_adder = lcg_adder_low + (uint64_t(lcg_adder_high) << 32);
	lcg_state = fast_forward_lcg64(how_far, lcg_state, 0x0000000100000001ULL, lcg_adder);
	lcg_low = uint32_t(lcg_state);
	lcg_high = uint32_t(lcg_state >> 32);
}
void PractRand::RNGs::Raw::xsm32::seek_backward(uint64_t how_far) {
	seek_forward(0 - how_far);
}


uint64_t PractRand::RNGs::Raw::xsm64::raw64() {
	const uint64_t K = 0xA3EC647659359ACDULL;//1010001111101100011001000111011001011001001101011001101011001101

	uint64_t tmp = lcg_high ^ std::rotl(lcg_high + lcg_low, 16); tmp ^= std::rotl(tmp + lcg_adder_high, 40);
	tmp *= K;
	step_forwards();
	tmp = tmp ^ std::rotl(tmp + lcg_high, 32);
	tmp *= K;
	tmp ^= tmp >> 32;
	return tmp;
}
void PractRand::RNGs::Raw::xsm64::step_backwards() {
	bool carry = lcg_low < lcg_adder_low;
	lcg_low -= lcg_adder_low;
	lcg_high -= lcg_low + lcg_adder_high + (carry ? 1 : 0);
}
void PractRand::RNGs::Raw::xsm64::step_forwards() {
	uint64_t tmp = lcg_low + lcg_adder_high;
	lcg_low += lcg_adder_low;
	lcg_high += tmp + ((lcg_low < lcg_adder_low) ? 1 : 0);
}
void PractRand::RNGs::Raw::xsm64::seed(uint64_t seed_low, uint64_t seed_high) {
	//no two seeds are closer than 2**127 from each other on the same cycle

	lcg_adder_low = (seed_low << 1) | 1;
	lcg_adder_high = (seed_low >> 63) | (seed_high << 1);//every bit of seed except the highest bit gets used in the adder

	lcg_low = lcg_adder_low;
	lcg_high = lcg_adder_high ^ ((seed_high >> 63) << 63);//and the highest bit of seed is used to determine which end of the cycle we start at
	step_forwards();
}
void PractRand::RNGs::Raw::xsm64::seed(vRNG* seeder_rng) {
	//guarantees that no two distinct seeded states are within 2**95 of each other on the same cycle ; 2**160 distinct seeded states are possible
	uint64_t s1 = 0, s2 = 0, s3 = 0;
	s1 = seeder_rng->raw64();
	s2 = seeder_rng->raw64();
	s3 = seeder_rng->raw64();
	seed(s1, s2);
	lcg_high += s3 << 31;
}
void PractRand::RNGs::Raw::xsm64::walk_state(StateWalkingObject* walker) {
	walker->handle(lcg_low);
	walker->handle(lcg_high);
	walker->handle(lcg_adder_low);
	walker->handle(lcg_adder_high);
	if (walker->is_clumsy() && !walker->is_read_only()) {
		lcg_adder_low |= 1;
	}
}
void PractRand::RNGs::Raw::xsm64::seek_forward (uint64_t how_far_low, uint64_t how_far_high) {
	if (!how_far_low && !how_far_high) return;
	fast_forward_lcg128(how_far_low, how_far_high,
		lcg_low, lcg_high, 1, 1, lcg_adder_low, lcg_adder_high
	);
}
void PractRand::RNGs::Raw::xsm64::seek_backward(uint64_t how_far_low, uint64_t how_far_high) {
	if (!how_far_low && !how_far_high) return;
	how_far_low = 1 + ~how_far_low;
	how_far_high = ~how_far_high + (how_far_low ? 0 : 1);
	seek_forward(how_far_low, how_far_high);
}

