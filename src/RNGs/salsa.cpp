
#include "PractRand/RNGs/salsa.h"
#include "PractRand/endian.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <bit>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <string>

using namespace PractRand;
using namespace PractRand::Internals;

/*
Salsa matrix structure:
00-03	constant0, seed0, seed1, seed2
04-07	seed3, constant1, position0, position1
08-11	IV0, IV1, constant2, seed0/4
12-15	seed1/5, seed2/6, seed3/7, constant3

The constants go on the diagonal, everything else fits its way in around that,
with the position & IV going between the first half of the seed and the second half of the seed.
If the seed is short then the 2nd half is treated as equal to the first half.
*/
constexpr int POS_INDEX0 = 8;
constexpr int POS_INDEX1 = 9;
constexpr int IV_INDEX0 = 6;
constexpr int IV_INDEX1 = 7;
constexpr int SEED_INDEX_A = 1;
constexpr int SEED_INDEX_B = 11;
constexpr int CONST_INDEX_0 = 0;
constexpr int CONST_INDEX_1 = 5;
constexpr int CONST_INDEX_2 = 10;
constexpr int CONST_INDEX_3 = 15;
constexpr int POSITION_OVERFLOW_INDEX = CONST_INDEX_0;

//polymorphic:
PRACTRAND_POLYMORPHIC_RNG_BASICS_C32(salsa)
std::string PractRand::RNGs::Polymorphic::salsa::get_name() const {
	std::ostringstream tmp;
	tmp << "salsa(" << implementation.get_rounds() << ")";
	return tmp.str();
}
void PractRand::RNGs::Polymorphic::salsa::seed(uint64_t s) {implementation.seed(s);}
void PractRand::RNGs::Polymorphic::salsa::seed(uint32_t seed_and_iv[10], bool extend_cycle_) {implementation.seed(seed_and_iv, extend_cycle_);}
void PractRand::RNGs::Polymorphic::salsa::seed_short(uint32_t seed_and_iv[6], bool extend_cycle_) {implementation.seed(seed_and_iv, extend_cycle_);}
void PractRand::RNGs::Polymorphic::salsa::seek_forward128 (uint64_t how_far_low64, uint64_t how_far_high64) {implementation.seek_forward (how_far_low64, how_far_high64);}
void PractRand::RNGs::Polymorphic::salsa::seek_backward128(uint64_t how_far_low64, uint64_t how_far_high64) {implementation.seek_backward(how_far_low64, how_far_high64);}
void PractRand::RNGs::Polymorphic::salsa::set_rounds(int rounds_) {implementation.set_rounds(rounds_);}
int PractRand::RNGs::Polymorphic::salsa::get_rounds() const {return implementation.get_rounds();}


//raw:
PractRand::RNGs::Raw::salsa::~salsa() {std::memset(this, 0, sizeof(*this));}
static void salsa_mix_core(uint32_t& a, uint32_t& b, uint32_t& c, uint32_t& d) {
	b ^= std::rotl(a + d, 7);
	c ^= std::rotl(b + a, 9);
	d ^= std::rotl(c + b, 13);
	a ^= std::rotl(d + c, 18);
}
//	const uint32_t *constants = short_seed ? salsa_short_seed_constants : salsa_long_seed_constants;
void PractRand::RNGs::Raw::salsa::_core() {
	for (int i = 0; i < 16; i++) outbuf[i] = state[i];
	if (extend_cycle) outbuf[POSITION_OVERFLOW_INDEX] += position_overflow;
#define QUARTERROUND(i1,i2,i3,i4) salsa_mix_core( outbuf[i1], outbuf[i2], outbuf[i3], outbuf[i4] );
	for (int round = 1; round < rounds; round+=2) {
		QUARTERROUND( 0, 4, 8,12)//columns
		QUARTERROUND( 5, 9,13, 1)
		QUARTERROUND(10,14, 2, 6)
		QUARTERROUND(15, 3, 7,11)
		QUARTERROUND( 0, 1, 2, 3)//rows
		QUARTERROUND( 5, 6, 7, 4)
		QUARTERROUND(10,11, 8, 9)
		QUARTERROUND(15,12,13,14)
	}
	if (rounds & 1) {
		QUARTERROUND( 0, 4, 8,12)//columns
		QUARTERROUND( 5, 9,13, 1)
		QUARTERROUND(10,14, 2, 6)
		QUARTERROUND(15, 3, 7,11)
	}
#undef QUARTERROUND
	for (int i = 0; i < 16; i++) outbuf[i] += state[i];
	if (extend_cycle) outbuf[POSITION_OVERFLOW_INDEX] += position_overflow;
}
uint32_t PractRand::RNGs::Raw::salsa::_refill_and_raw32() {
	_advance_1();
	_core();
	used = 1;
	return outbuf[0];
}
void PractRand::RNGs::Raw::salsa::_advance_1() {
	if (!++state[POS_INDEX0]) {
		if (!++state[POS_INDEX1]) position_overflow++;
	}
}
//void PractRand::RNGs::Raw::salsa::_reverse_1();
void PractRand::RNGs::Raw::salsa::_set_position(uint64_t low, uint64_t high) {
	used = low & 15;
	low >>= 4;
	low |= high << 60;
	high >>= 4;
	state[POS_INDEX0] = uint32_t(low);
	state[POS_INDEX1] = uint32_t(low >> 32);
	position_overflow = uint32_t(high);
	_core();
}
void PractRand::RNGs::Raw::salsa::_get_position(uint64_t& low, uint64_t& high) const {
	low = used + (uint64_t(state[POS_INDEX0]) << 4) + (uint64_t(state[POS_INDEX1]) << 36);
	high = (state[POS_INDEX1] >> 28) + (uint64_t(position_overflow) << 4);
}
void PractRand::RNGs::Raw::salsa::seed(uint64_t s) {
	uint32_t seed_and_iv[10] = {0};
	seed_and_iv[0] = uint32_t(s);
	seed_and_iv[1] = uint32_t(s >> 32);
	seed(seed_and_iv, true);
}
const uint32_t salsa_short_seed_constants[4] = {
	(uint32_t(101) << 0) + (uint32_t(120) << 8) + (uint32_t(112) << 16) + (uint32_t( 97) << 24),
	(uint32_t(110) << 0) + (uint32_t(100) << 8) + (uint32_t( 32) << 16) + (uint32_t( 49) << 24),
	(uint32_t( 54) << 0) + (uint32_t( 45) << 8) + (uint32_t( 98) << 16) + (uint32_t(121) << 24),
	(uint32_t(116) << 0) + (uint32_t(101) << 8) + (uint32_t( 32) << 16) + (uint32_t(107) << 24),
};
const uint32_t salsa_long_seed_constants[4] = {
	(uint32_t(101) << 0) + (uint32_t(120) << 8) + (uint32_t(112) << 16) + (uint32_t( 97) << 24),
	(uint32_t(110) << 0) + (uint32_t(100) << 8) + (uint32_t( 32) << 16) + (uint32_t( 51) << 24),
	(uint32_t( 50) << 0) + (uint32_t( 45) << 8) + (uint32_t( 98) << 16) + (uint32_t(121) << 24),
	(uint32_t(116) << 0) + (uint32_t(101) << 8) + (uint32_t( 32) << 16) + (uint32_t(107) << 24),
};
void PractRand::RNGs::Raw::salsa::seed(const uint32_t seed_and_iv[10], bool extend_cycle_) {
	const uint32_t* constants = salsa_long_seed_constants;
	state[CONST_INDEX_0] = constants[0];
	state[CONST_INDEX_1] = constants[1];
	state[CONST_INDEX_2] = constants[2];
	state[CONST_INDEX_3] = constants[3];
	position_overflow = 0;
	extend_cycle = extend_cycle_;
	for (int i = 0; i < 4; i++) state[SEED_INDEX_A + i] = seed_and_iv[i];
	for (int i = 0; i < 4; i++) state[SEED_INDEX_B + i] = seed_and_iv[i+4];
	state[POS_INDEX0] = 0;
	state[POS_INDEX1] = 0;
	state[IV_INDEX0] = seed_and_iv[8];
	state[IV_INDEX1] = seed_and_iv[9];
	_core();
	used = 0;
}
void PractRand::RNGs::Raw::salsa::seed_short(const uint32_t seed_and_iv[6], bool extend_cycle_) {
	const uint32_t* constants = salsa_short_seed_constants;
	state[CONST_INDEX_0] = constants[0];
	state[CONST_INDEX_1] = constants[1];
	state[CONST_INDEX_2] = constants[2];
	state[CONST_INDEX_3] = constants[3];
	position_overflow = 0;
	extend_cycle = extend_cycle_;
	for (int i = 0; i < 4; i++) state[SEED_INDEX_A + i] = seed_and_iv[i];
	for (int i = 0; i < 4; i++) state[SEED_INDEX_B + i] = seed_and_iv[i];
	state[POS_INDEX0] = 0;
	state[POS_INDEX1] = 0;
	state[IV_INDEX0] = seed_and_iv[4];
	state[IV_INDEX1] = seed_and_iv[5];
	_core();
	used = 0;
}
void PractRand::RNGs::Raw::salsa::walk_state(StateWalkingObject* walker) {
	for (auto& i : state) walker->handle(i);
	walker->handle(used);
	walker->handle(extend_cycle);
	if (extend_cycle) walker->handle(position_overflow);
	if (walker->is_seeder()) {
		const uint32_t* constants = salsa_long_seed_constants;
		state[CONST_INDEX_0] = constants[0];
		state[CONST_INDEX_1] = constants[1];
		state[CONST_INDEX_2] = constants[2];
		state[CONST_INDEX_3] = constants[3];
		state[POS_INDEX0] = 0;
		state[POS_INDEX1] = 0;
		state[IV_INDEX0] = 0;
		state[IV_INDEX1] = 0;
		extend_cycle = true;
		position_overflow = 0;
	}
	else { walker->handle(rounds); }

	if (used >= 16) used = 16;
	if (!walker->is_read_only()) {
		_core();
		used &= 15;
	}
}
void PractRand::RNGs::Raw::salsa::seek_forward (uint64_t how_far_low, uint64_t how_far_high) {
	uint64_t pos_low = 0, pos_high = 0;
	_get_position(pos_low, pos_high);
	uint64_t new_pos_low = pos_low + how_far_low;
	if (new_pos_low < pos_low) how_far_high++;
	uint64_t new_pos_high = pos_high + how_far_high;
	_set_position(new_pos_low, new_pos_high);
}
void PractRand::RNGs::Raw::salsa::seek_backward(uint64_t how_far_low, uint64_t how_far_high) {
	seek_forward(~how_far_low, ~how_far_high);
	raw32();
}
void PractRand::RNGs::Raw::salsa::set_rounds(int rounds_) {
	if (rounds_ < 1 || rounds_ > 255) issue_error("salsa rounds out of range");
	if (rounds == rounds_) return;
	rounds = rounds_;
	//_core();
}
/*
static void test_salsa ( uint32_t rounds, const uint32_t *seed_and_iv, bool short_seed, uint32_t expected0, uint32_t index, uint32_t expected1) {
	PractRand::RNGs::Raw::salsa rng;
	rng.set_rounds(rounds);
	if (!short_seed) rng.seed(seed_and_iv, false);
	else rng.seed_short(seed_and_iv, false);
	uint32_t observed0 = rng.raw32();
	uint32_t observed1;
	if (!index) observed1 = observed0;
	else {
		for (uint32_t i = 1; i < index; i++) rng.raw32();
		observed1 = rng.raw32();
	}

	if (expected0 != observed0 || expected1 != observed1) {
		PractRand::issue_error("Salsa CSPRNG failed self-test");
	}
}
*/
void PractRand::RNGs::Raw::salsa::self_test() {
	PractRand::RNGs::Raw::salsa engine;
	uint32_t seed_and_iv[10] = {0};
	for (int i = 0; i < 4; i++) seed_and_iv[i+0] = i * 0x04040404 + 0x04030201;
	for (int i = 0; i < 4; i++) seed_and_iv[i+4] = i * 0x04040404 + 0xCCCBCAC9;
	for (int i = 0; i < 2; i++) seed_and_iv[i+8] = i * 0x04040404 + 0x68676665;
	engine.set_rounds(20);
	engine.seed(seed_and_iv, false);
	uint64_t N = (uint64_t(109)<<0)+(uint64_t(110)<<8)+(uint64_t(111)<<16)+(uint64_t(112)<<24)+(uint64_t(113)<<32)+(uint64_t(114)<<40)+(uint64_t(115)<<48)+(uint64_t(116)<<56);
	engine.seek_forward( N << 4, N >> 60);
	uint64_t E = (uint64_t(69)<<0)+(uint64_t(37)<<8)+(uint64_t(68)<<16)+(uint64_t(39)<<24)+(uint64_t(41)<<32)+(uint64_t(15)<<40)+(uint64_t(107)<<48)+(uint64_t(193)<<56);
	if (uint32_t(E) != engine.raw32()) issue_error("salsa::self_test() failed\n");
	if (uint32_t(E>>32) != engine.raw32()) issue_error("salsa::self_test() failed\n");

	engine.set_rounds(20);
	for (int i = 0; i < 6; i++) seed_and_iv[i] = 0;
	engine.seed_short(seed_and_iv, false);
	uint64_t E2 = 0x4c12ebcfaead1365ULL;
	if (uint32_t(E2) != engine.raw32()) issue_error("salsa::self_test() failed a\n");
	if (uint32_t(E2>>32) != engine.raw32()) issue_error("salsa::self_test() failed b\n");
}


