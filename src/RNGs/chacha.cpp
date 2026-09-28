
#include "PractRand/RNGs/chacha.h"
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
ChaCha matrix structure:
00-03	constant0, constant1, constant2, constant3
04-07	seed0, seed1, seed2, seed3
08-11	seed0/4, seed1/5, seed2/6, seed3/7
12-15	position0, position1, IV0, IV1

top row is constants
2nd row is the first 16 bytes of seed
3rd row is the rest of the seed (for long seeds) or a duplicate of the seed (for short seeds)
4th row is the position & initialization vector
*/
constexpr int POS_INDEX0 = 12 - 4;
constexpr int POS_INDEX1 = 13 - 4;
constexpr int IV_INDEX0 = 14 - 4;
constexpr int IV_INDEX1 = 15 - 4;
constexpr int POSITION_OVERFLOW_INDEX = 0;

//polymorphic:
PRACTRAND_POLYMORPHIC_RNG_BASICS_C32(chacha)
std::string PractRand::RNGs::Polymorphic::chacha::get_name() const {
	std::ostringstream tmp;
	tmp << "chacha(" << implementation.get_rounds() << ")";
	return tmp.str();
}
void PractRand::RNGs::Polymorphic::chacha::seed(uint64_t s) {implementation.seed(s);}
void PractRand::RNGs::Polymorphic::chacha::seed(uint32_t seed_and_iv[10], bool extend_cycle_) {implementation.seed(seed_and_iv, extend_cycle_);}
void PractRand::RNGs::Polymorphic::chacha::seed_short(uint32_t seed_and_iv[6], bool extend_cycle_) {implementation.seed_short(seed_and_iv, extend_cycle_);}
void PractRand::RNGs::Polymorphic::chacha::seek_forward128 (uint64_t how_far_low64, uint64_t how_far_high64) {implementation.seek_forward (how_far_low64, how_far_high64);}
void PractRand::RNGs::Polymorphic::chacha::seek_backward128(uint64_t how_far_low64, uint64_t how_far_high64) {implementation.seek_backward(how_far_low64, how_far_high64);}
void PractRand::RNGs::Polymorphic::chacha::set_rounds(int rounds_) {implementation.set_rounds(rounds_);}
int PractRand::RNGs::Polymorphic::chacha::get_rounds() const {return implementation.get_rounds();}

//raw:
PractRand::RNGs::Raw::chacha::~chacha() {explicit_bzero(this, sizeof(*this));}
static constexpr uint32_t chacha_long_seed_constants[4] = {
	//"expand 32-byte k"
	(uint32_t('e') << 0) + (uint32_t('x') << 8) + (uint32_t('p') << 16) + (uint32_t('a') << 24),
	(uint32_t('n') << 0) + (uint32_t('d') << 8) + (uint32_t(' ') << 16) + (uint32_t('3') << 24),
	(uint32_t('2') << 0) + (uint32_t('-') << 8) + (uint32_t('b') << 16) + (uint32_t('y') << 24),
	(uint32_t('t') << 0) + (uint32_t('e') << 8) + (uint32_t(' ') << 16) + (uint32_t('k') << 24)
};
static constexpr uint32_t chacha_short_seed_constants[4] = {
	//"expand 16-byte k"
	(uint32_t('e') << 0) + (uint32_t('x') << 8) + (uint32_t('p') << 16) + (uint32_t('a') << 24),
	(uint32_t('n') << 0) + (uint32_t('d') << 8) + (uint32_t(' ') << 16) + (uint32_t('1') << 24),
	(uint32_t('6') << 0) + (uint32_t('-') << 8) + (uint32_t('b') << 16) + (uint32_t('y') << 24),
	(uint32_t('t') << 0) + (uint32_t('e') << 8) + (uint32_t(' ') << 16) + (uint32_t('k') << 24)
};

void PractRand::RNGs::Raw::chacha::_advance_1() {
	if (!++state[POS_INDEX0]) {
		if (!++state[POS_INDEX1]) position_overflow++;
	}
}
void PractRand::RNGs::Raw::chacha::_set_position(uint64_t low, uint64_t high) {
	used = low & 15;
	low >>= 4;
	low |= high << 60;
	high >>= 4;
	state[POS_INDEX0] = uint32_t(low);
	state[POS_INDEX1] = uint32_t(low >> 32);
	position_overflow = uint32_t(high);
	_core();
}
void PractRand::RNGs::Raw::chacha::_get_position(uint64_t& low, uint64_t& high) const {
	low = used + (uint64_t(state[POS_INDEX0]) << 4) + (uint64_t(state[POS_INDEX1]) << 36);
	high = (state[POS_INDEX1] >> 28) + (uint64_t(position_overflow) << 4);
}
static void chacha_mix_core(uint32_t& a, uint32_t& b, uint32_t& c, uint32_t& d) {
	a += b; d = std::rotl(d ^ a, 16);
	c += d; b = std::rotl(b ^ c, 12);
	a += b; d = std::rotl(d ^ a, 8);
	c += d; b = std::rotl(b ^ c, 7);
}
void PractRand::RNGs::Raw::chacha::_core() {
	const uint32_t* constants = short_seed ? chacha_short_seed_constants : chacha_long_seed_constants;

	for (int i = 0; i < 4; i++) outbuf[i] = constants[i];
	for (int i = 4; i < 16; i++) outbuf[i] = state[i-4];
	if (extend_cycle) outbuf[POSITION_OVERFLOW_INDEX] += position_overflow;
	for (int round = 1; round < rounds; round+=2) {
		chacha_mix_core(outbuf[0], outbuf[4], outbuf[8], outbuf[12]);
		chacha_mix_core(outbuf[1], outbuf[5], outbuf[9], outbuf[13]);
		chacha_mix_core(outbuf[2], outbuf[6], outbuf[10], outbuf[14]);
		chacha_mix_core(outbuf[3], outbuf[7], outbuf[11], outbuf[15]);
		chacha_mix_core(outbuf[0], outbuf[5], outbuf[10], outbuf[15]);
		chacha_mix_core(outbuf[1], outbuf[6], outbuf[11], outbuf[12]);
		chacha_mix_core(outbuf[2], outbuf[7], outbuf[8], outbuf[13]);
		chacha_mix_core(outbuf[3], outbuf[4], outbuf[9], outbuf[14]);
	}
	if (rounds & 1) {
		chacha_mix_core(outbuf[0], outbuf[5], outbuf[10], outbuf[15]);
		chacha_mix_core(outbuf[1], outbuf[6], outbuf[11], outbuf[12]);
		chacha_mix_core(outbuf[2], outbuf[7], outbuf[8], outbuf[13]);
		chacha_mix_core(outbuf[3], outbuf[4], outbuf[9], outbuf[14]);
	}
	for (int i = 0; i < 4; i++) outbuf[i] += constants[i];
	for (int i = 4; i < 16; i++) outbuf[i] += state[i-4];
	if (extend_cycle) outbuf[POSITION_OVERFLOW_INDEX] += position_overflow;
}
uint32_t PractRand::RNGs::Raw::chacha::_refill_and_raw32() {
	_advance_1();
	_core();
	used = 1;
	return outbuf[0];
}
void PractRand::RNGs::Raw::chacha::seed(uint64_t s) {
	uint32_t seed_and_iv[10] = {0};
	seed_and_iv[0] = uint32_t(s);
	seed_and_iv[1] = uint32_t(s >> 32);
	seed(seed_and_iv, true);
}
void PractRand::RNGs::Raw::chacha::seed(const uint32_t seed_and_iv[10], bool extend_cycle_) {
	short_seed = false;
	position_overflow = 0;
	extend_cycle = extend_cycle_;
	for (int i = 0; i < 8; i++) state[i] = seed_and_iv[i];
	state[POS_INDEX0] = 0;
	state[POS_INDEX1] = 0;
	state[IV_INDEX0] = seed_and_iv[8];
	state[IV_INDEX1] = seed_and_iv[9];
	_core();
	used = 0;
}
void PractRand::RNGs::Raw::chacha::seed_short(const uint32_t seed_and_iv[6], bool extend_cycle_) {
	short_seed = true;
	position_overflow = 0;
	extend_cycle = extend_cycle_;
	for (int i = 0; i < 8; i++) state[i] = seed_and_iv[i & 3];
	state[POS_INDEX0] = 0;
	state[POS_INDEX1] = 0;
	state[IV_INDEX0] = seed_and_iv[4];
	state[IV_INDEX1] = seed_and_iv[5];
	_core();
	used = 0;
}
void PractRand::RNGs::Raw::chacha::walk_state(StateWalkingObject* walker) {
	for (auto& i : state) walker->handle(i);
	walker->handle(used);
	walker->handle(extend_cycle);
	if (extend_cycle) walker->handle(position_overflow);
	walker->handle(short_seed);
	if (walker->is_seeder()) {
		short_seed = false;
		extend_cycle = true;
		position_overflow = 0;
		state[POS_INDEX0] = 0;
		state[POS_INDEX1] = 0;
		state[IV_INDEX0] = 0;
		state[IV_INDEX1] = 0;
	}
	else { walker->handle(rounds); }


	if (used >= 16) used = 16;
	if (!walker->is_read_only()) {
		_core();
		used &= 15;
	}
}
void PractRand::RNGs::Raw::chacha::seek_forward (uint64_t how_far_low, uint64_t how_far_high) {
	uint64_t pos_low = 0, pos_high = 0;
	_get_position(pos_low, pos_high);
	uint64_t new_pos_low = pos_low + how_far_low;
	if (new_pos_low < pos_low) how_far_high++;
	uint64_t new_pos_high = pos_high + how_far_high;
	_set_position(new_pos_low, new_pos_high);
}
void PractRand::RNGs::Raw::chacha::seek_backward(uint64_t how_far_low, uint64_t how_far_high) {
	seek_forward(~how_far_low, ~how_far_high);
	raw32();
}
void PractRand::RNGs::Raw::chacha::set_rounds(int rounds_) {
	if (rounds_ < 1 || rounds_ > 255) issue_error("chacha rounds out of range");
	if (rounds == rounds_) return;
	rounds = rounds_;
	//_core();
}
static void test_chacha ( uint32_t rounds, const uint32_t* seed_and_iv, bool short_seed, uint32_t expected0, uint32_t index, uint32_t expected1) {
	PractRand::RNGs::Raw::chacha rng;
	rng.set_rounds(rounds);
	if (!short_seed) rng.seed(seed_and_iv, false);
	else rng.seed_short(seed_and_iv, false);
	uint32_t observed0 = rng.raw32();
	uint32_t observed1 = 0;
	if (!index) { observed1 = observed0; }
	else {
		for (uint32_t i = 1; i < index; i++) rng.raw32();
		observed1 = rng.raw32();
	}

	if (expected0 != observed0 || expected1 != observed1) {
		PractRand::issue_error("ChaCha CSPRNG failed self-test");
	}
}
void PractRand::RNGs::Raw::chacha::self_test() {
	uint32_t seed_and_iv[10] = {0};
	test_chacha(  8, seed_and_iv, false, 0x2fef003e, 16, 0x0dfaaed2);
	test_chacha( 12, seed_and_iv, false, 0x6a9af49b, 16, 0x4188d50b);
	test_chacha( 20, seed_and_iv, false, 0xade0b876, 16, 0xbee7079f);
	seed_and_iv[0] = 0x80;
	test_chacha(  8, seed_and_iv,  true, 0x1ee8b1be, 16, 0x7d720502);
	seed_and_iv[4] = 0x01;
	test_chacha(  8, seed_and_iv,  true, 0x3dc23b69, 16, 0xc9814bee);
	test_chacha(  8, seed_and_iv, false, 0x6bcbf0bb, 16, 0x281b6aa7);
	seed_and_iv[9] = 0x01;
	test_chacha(  8, seed_and_iv, false, 0xcc135949, 16, 0x82f734d3);
	seed_and_iv[8] = 0x12;
	test_chacha(  8, seed_and_iv, false, 0xd5a08851, 3, 0xab4a48a9);
}




