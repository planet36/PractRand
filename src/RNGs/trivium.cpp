#include "PractRand/RNGs/trivium.h"
#include "PractRand/endian.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"
#include <bit>
#include <cstdint>
#include <cstring>
#include <string>

using namespace PractRand;


//polymorphic:
PRACTRAND_POLYMORPHIC_RNG_BASICS_C64(trivium)
void PractRand::RNGs::Polymorphic::trivium::seed(uint64_t s) { implementation.seed(s); }
void PractRand::RNGs::Polymorphic::trivium::seed_fast(uint64_t s) { implementation.seed_fast(s, s); }
void PractRand::RNGs::Polymorphic::trivium::seed(vRNG* seeder_rng) { implementation.seed(seeder_rng); }
void PractRand::RNGs::Polymorphic::trivium::seed(const uint8_t* seed_and_iv, int length) { implementation.seed(seed_and_iv, length); }
std::string PractRand::RNGs::Polymorphic::trivium::get_name() const {return "trivium";}

static uint64_t shift_array64( uint64_t vec[2], unsigned long bits ) {
	bits -= 64;
	if (!(bits % 64)) return vec[bits/64];
	return (vec[bits / 64] << (bits & 63)) | (vec[1 + bits / 64] >> (64-(bits & 63)));
}
//raw:
PractRand::RNGs::Raw::trivium::~trivium() {std::memset(this, 0, sizeof(*this));}
uint64_t PractRand::RNGs::Raw::trivium::raw64() {//LOCKED, do not change
	uint64_t tmp_a = shift_array64(c, 66) ^ shift_array64(c,111);
	uint64_t tmp_b = shift_array64(a, 66) ^ shift_array64(a, 93);
	uint64_t tmp_c = shift_array64(b, 69) ^ shift_array64(b, 84);
	uint64_t new_a = tmp_a ^ shift_array64(a, 69) ^ (shift_array64(c,110) & shift_array64(c,109));
	uint64_t new_b = tmp_b ^ shift_array64(b, 78) ^ (shift_array64(a, 92) & shift_array64(a, 91));
	uint64_t new_c = tmp_c ^ shift_array64(c, 87) ^ (shift_array64(b, 83) & shift_array64(b, 82));
	a[1] = a[0]; a[0] = new_a;
	b[1] = b[0]; b[0] = new_b;
	c[1] = c[0]; c[0] = new_c;

	return tmp_a ^ tmp_b ^ tmp_c;
}
void PractRand::RNGs::Raw::trivium::seed(uint64_t s) {//LOCKED, do not change
	//Triviums standard seeding algorithm adapted to PractRand interface
	uint8_t vec[8];
	for (auto& i : vec) {
		i = uint8_t(s);
		s >>= 8;
	}
	seed(vec, 8);
}
void PractRand::RNGs::Raw::trivium::seed_fast(uint64_t s1, uint64_t s2, int quality) {//LOCKED, do not change
	//a non-standard simplified 128-bit seeding algorithm
	a[0] = s1; a[1] = 0;
	b[0] = s2; b[1] = 0;
	c[0] = 0;
	c[1] = uint64_t(7) << (128-111);
	for (int i = 0; i < quality; i++) raw64();
}
void PractRand::RNGs::Raw::trivium::seed(vRNG* seeder_rng) {//LOCKED, do not change
	uint64_t s1 = seeder_rng->raw64();
	uint64_t s2 = seeder_rng->raw64();
	seed_fast(s1, s2, 3);
}
void PractRand::RNGs::Raw::trivium::seed(const uint8_t* seed_and_iv, int length) {//LOCKED, do not change
	//standard algorithm for Trivium, not a good match for PractRand
	if (length > 20) issue_error("trivium seeded with invalid length");

	union SeedVec{
		uint8_t as8[16];
		uint64_t as64[2];
	};
	SeedVec s{};
	int elen = length > 10 ? 10 : length;
	for (int i = 0; i < 6; i++) s.as8[i] = 0;
	for (int i = 0; i < elen; i++) s.as8[i+6] = seed_and_iv[i];
	for (int i = elen; i < 10; i++) s.as8[i+6] = 0;
	for (int i = 0; i < 2; i++) this->a[i] = s.as64[1-i];

	length -= 10;
	seed_and_iv += 10;
	length = length > 0 ? length : 0;
	SeedVec iv{};
	for (int i = 0; i < 16-length; i++) iv.as8[i] = 0;
	for (int i = 16-length; i < 16; i++) iv.as8[i] = seed_and_iv[i-(16-length)];
	for (int i = 0; i < 2; i++) this->b[i] = iv.as64[1-i];

	c[0] = 0;
	c[1] = uint64_t(7) << (128-111);

	//for (int i = 0; i < 1152; i+=1) raw1();
	for (int i = 0; i < 18; i++) raw64();
	/*
		(# of outputs discarded) vs (log2 of # of seeds needed to detect interseed correlation)
			4 - 13
			5 - 15
			6 - 25
			7 - 30
			8 -
	*/
}
void PractRand::RNGs::Raw::trivium::walk_state(StateWalkingObject* walker) {
	//LOCKED, do not change
	walker->handle(a[0]);
	walker->handle(a[1]);
	walker->handle(b[0]);
	walker->handle(b[1]);
	walker->handle(c[0]);
	walker->handle(c[1]);
}

static void validate_trivium_result(uint64_t output, uint64_t reference) {
	//convert format of reference
	//	it will always be in the wrong endianness, regardless of platform
	//	due to how it was cut & pasted
	reference = std::byteswap(reference);
	//raise an error if they don't match
	if (output != reference) {
		//uint64_t r = output;
		//printf("%02x %02x %02x %02x %02x %02x %02x %02x\n", uint8_t(r>>0), uint8_t(r>>8), uint8_t(r>>16), uint8_t(r>>24), uint8_t(r>>32), uint8_t(r>>40), uint8_t(r>>48), uint8_t(r>>56));
		//for (int i = 0; i < 16; i++) printf("%X", int(uint8_t(r >> (i*4)) & 15));
		//printf("\n");
		PractRand::issue_error("PractRand::RNGs::Raw::trivium failed validation");
	}
}
void PractRand::RNGs::Raw::trivium::self_test() {
	Raw::trivium rng{};
	uint8_t seed_and_iv[10+10] = {0};
	rng.seed(seed_and_iv, 14);
	validate_trivium_result(rng.raw64(), 0xFBE0BF265859051BULL);
	seed_and_iv[0] = 0x80; rng.seed(seed_and_iv, 14); seed_and_iv[0] = 0x00;
	validate_trivium_result(rng.raw64(), 0x38EB86FF730D7A9CULL);
	seed_and_iv[9] = 0x80; rng.seed(seed_and_iv, 14); seed_and_iv[9] = 0x00;
	validate_trivium_result(rng.raw64(), 0x5D492E77F8FE62D7ULL);
	seed_and_iv[13] = 0x10; rng.seed(seed_and_iv, 14); seed_and_iv[13] = 0x00;
	validate_trivium_result(rng.raw64(), 0xB0820A503ABB0329ULL);
	seed_and_iv[17] = 0x01; rng.seed(seed_and_iv, 18); seed_and_iv[17] = 0x00;
	validate_trivium_result(rng.raw64(), 0x9A5C56169E7FA406ULL);
}
