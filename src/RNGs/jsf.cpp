#include <string>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>

//#include "PractRand/RNGs/jsf8.h"
//#include "PractRand/RNGs/jsf16.h"
#include "PractRand/RNGs/jsf32.h"
#include "PractRand/RNGs/jsf64.h"

using namespace PractRand;

//polymorphic:
PRACTRAND__POLYMORPHIC_RNG_BASICS_C64(jsf64)
void PractRand::RNGs::Polymorphic::jsf64::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
void PractRand::RNGs::Polymorphic::jsf64::seed_fast(Uint64 s) {implementation.seed_fast(s);}
std::string PractRand::RNGs::Polymorphic::jsf64::get_name() const {return "jsf64";}

PRACTRAND__POLYMORPHIC_RNG_BASICS_C32(jsf32)
void PractRand::RNGs::Polymorphic::jsf32::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
void PractRand::RNGs::Polymorphic::jsf32::seed_fast(Uint64 s) {implementation.seed_fast(s);}
std::string PractRand::RNGs::Polymorphic::jsf32::get_name() const {return "jsf32";}

//raw:
Uint32 PractRand::RNGs::Raw::jsf32::raw32() {//LOCKED, do not change
	Uint32 e = a - ((b << 27) | (b >> 5));
	a = b ^ ((c << 17) | (c >> 15));
	b = c + d;
	c = d + e;
	d = e + a;
	return d;
}
void PractRand::RNGs::Raw::jsf32::seed(Uint64 seed_low, Uint64 seed_high) {//LOCKED, do not change
	//LOCKED, do not change
	//exception: changed in 0.81 to match Robert Jenkins code
	//note: only actually matches Robert Jenkins code for seeds that fit in a 32 bit value, as he only took a 32 bit seed
	//exception: changed in 0.96 to support 128 bit seeds - seeds that fit in a 64 bit value should produce the same result as before
	a = 0xf1ea5eed ^ Uint32(seed_low >> 32);
	b = Uint32(seed_low >> 0);
	c = b ^ Uint32(seed_low >> 32) ^ Uint32(seed_high >> 0);
	d = b ^ Uint32(seed_high >> 32);
	for (int i = 0; i < 20; i++) raw32();//20
	/*
		number of outputs skipped vs number of seeds needed to detect bias via expanded battery w/ extra folding:
		2: 2**10
		3: 2**12
		4: 2**18
		5: 2**23
		6: 2**35
		7: >= 2**42
		conclusion:
			The 20 outputs skipped by the standard algorithm are more than sufficient.  
			8 should be enough for most purposes, 
			12 should be enough for any purpose that could possibly be satisfied by seeding from a 64 bit integer.  
			(this PRNG won't be secure no matter how good seeding is, so the best we can do is 2**64 seeds producing uncorrelated results)
	*/
}
void PractRand::RNGs::Raw::jsf32::seed_fast(Uint64 s) {
	a = 0xf1ea5eed ^ Uint32(s >> 32);
	b = Uint32(s);
	c = b ^ Uint32(s >> 32);
	d = b;
	for (int i = 0; i < 8; i++) raw32();
}
void PractRand::RNGs::Raw::jsf32::seed(vRNG *seeder_rng) {//custom seeding
	a = b = seeder_rng->raw32();
	c = seeder_rng->raw32();
	d = ~c;
	for (int i = 0; i < 4; i++) raw32();//4
}
bool PractRand::RNGs::Raw::jsf32::is_state_bad() const {//added in 0.95 to consolidate checks for bad cycles in one place
	//this code checks for cycles of length 1
	//I should search for other short cycle lengths as well, but that's impractical at this time
	static const Uint32 bad_jsf32_state_values[6][4] = {
		{ 0, 0, 0, 0 },
		{ 0x77777777, 0x55555555, 0x11111111, 0x44444444 },
		{ 0x5591F2E3, 0x69EBA6CD, 0x2A171E3D, 0x3FD48890 },
		{ 0x47CB8D56, 0xAE9B35A7, 0x5C78F4A8, 0x522240FF },
		{ 0x71AAC8F9, 0x66B4F5D3, 0x1E950B8F, 0x481FEA44 },
		{ 0xAB23E5C6, 0xD3D74D9A, 0x542E3C7A, 0x7FA91120 },
	};
	//if ((d ^ (d >> 1)) & 0x80800001) return false;//not used - this doesn't filter out enough values, and requires extra operations
	if (!(d & (~bad_jsf32_state_values[0][3] & ~bad_jsf32_state_values[1][3] & ~bad_jsf32_state_values[2][3] & ~bad_jsf32_state_values[3][3]))) {
		if (a == bad_jsf32_state_values[0][0] && b == bad_jsf32_state_values[0][1] && c == bad_jsf32_state_values[0][2] && d == bad_jsf32_state_values[0][3]) return true;
		if (a == bad_jsf32_state_values[1][0] && b == bad_jsf32_state_values[1][1] && c == bad_jsf32_state_values[1][2] && d == bad_jsf32_state_values[1][3]) return true;
		if (a == bad_jsf32_state_values[2][0] && b == bad_jsf32_state_values[2][1] && c == bad_jsf32_state_values[2][2] && d == bad_jsf32_state_values[2][3]) return true;
		if (a == bad_jsf32_state_values[3][0] && b == bad_jsf32_state_values[3][1] && c == bad_jsf32_state_values[3][2] && d == bad_jsf32_state_values[3][3]) return true;
	}
	if (!(d & (~bad_jsf32_state_values[4][3] & ~bad_jsf32_state_values[5][3]))) {
		if (a == bad_jsf32_state_values[4][0] && b == bad_jsf32_state_values[4][1] && c == bad_jsf32_state_values[4][2] && d == bad_jsf32_state_values[4][3]) return true;
		if (a == bad_jsf32_state_values[5][0] && b == bad_jsf32_state_values[5][1] && c == bad_jsf32_state_values[5][2] && d == bad_jsf32_state_values[5][3]) return true;
	}
	return false;
}
void PractRand::RNGs::Raw::jsf32::seed(Uint32 seed1, Uint32 seed2, Uint32 seed3, Uint32 seed4) {//custom seeding
		//LOCKED, do not change
	//exception to the locked status - 
	//   when more bad cycles are found, more code might be added to prohibit them
	//exception: changed in 0.87 to to reduce correlation between similar seeds
	//exception: changed in 0.95 to add more zero-length cycles, and reorganize them
	a = seed1;
	b = seed2;
	c = seed3;
	d = seed4;
	if (is_state_bad()) d++;
	for (int i = 0; i < 12; i++) raw32();//12
}
void PractRand::RNGs::Raw::jsf32::walk_state(StateWalkingObject *walker) {
	//LOCKED, do not change
	//exception to the locked status - 
	//   when more bad cycles are found, more code might be added to prohibit them
	//exception: in 0.87 changed how the all-zeroes case is handled for consistency
	//exception: changed in 0.95 to add more zero-length cycles, and reorganize them
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
	walker->handle(d);
	if (is_state_bad()) d++;
}
void PractRand::RNGs::Raw::jsf32::self_test() {
	Raw::jsf32 rng;
	rng.seed(0ull);
	if (rng.raw32() != 0x1a9b6c07) issue_error("jsf32::self_test() failed(a0)");
	if (rng.raw32() != 0x9a550895) issue_error("jsf32::self_test() failed(a1)");
	if (rng.raw32() != 0xf12be876) issue_error("jsf32::self_test() failed(a2)");
	if (rng.raw32() != 0x0902ba19) issue_error("jsf32::self_test() failed(a3)");
	if (rng.raw32() != 0x20f1a244) issue_error("jsf32::self_test() failed(a4)");

	rng.seed(0x9a550895ull);
	if (rng.raw32() != 0x11eac194) issue_error("jsf32::self_test() failed(b0)");
	if (rng.raw32() != 0x42165952) issue_error("jsf32::self_test() failed(b1)");

	rng.seed(0x42165952, 0x12345678, 1234567890, 10092);
	if (rng.raw32() != 0xdde46ccb) issue_error("jsf32::self_test() failed(c0)");
	if (rng.raw32() != 0xa793b5e4) issue_error("jsf32::self_test() failed(c1)");
	if (rng.raw32() != 0xe882f402) issue_error("jsf32::self_test() failed(c2)");
	if (rng.raw32() != 0x2edff2e8) issue_error("jsf32::self_test() failed(c3)");
	if (rng.raw32() != 0x119346a9) issue_error("jsf32::self_test() failed(c4)");

	//rng.seed(0, 0, 0, 0);
	//rng.seed(0x77777777, 0x55555555, 0x11111111, 0x44444444);
	//rng.seed(0x5591F2E3, 0x69EBA6CD, 0x2A171E3D, 0x3FD48890);
	//rng.seed(0x47CB8D56, 0xAE9B35A7, 0x5C78F4A8, 0x522240FF);
	//rng.seed(0x71AAC8F9, 0x66B4F5D3, 0x1E950B8F, 0x481FEA44);
	//rng.seed(0xAB23E5C6, 0xD3D74D9A, 0x542E3C7A, 0x7FA91120);
}


Uint64 PractRand::RNGs::Raw::jsf64::raw64() {
	//LOCKED, do not change
	Uint64 e = a - ((b << 39) | (b >> 25));
	a = b ^ ((c << 11) | (c >> 53));
	b = c + d;
	c = d + e;
	d = e + a;
	return d;
}
void PractRand::RNGs::Raw::jsf64::seed(Uint64 seed_low, Uint64 seed_high) {
	//LOCKED, do not change
	//changed in v0.96 to support 128 bit seeds ; seeds that fit in 64 bits should produce the same results as before
	a = 0xf1ea5eed;
	b = c = d = seed_low;
	c ^= seed_high;
	d ^= (seed_high >> 32) | (seed_high << 32);
	for (int i = 0; i < 20; i++) raw64();
}
void PractRand::RNGs::Raw::jsf64::seed_fast(Uint64 s) {
	a = 0xf1ea5eed;
	b = c = d = s;
	for (int i = 0; i < 8; i++) raw64();
}
void PractRand::RNGs::Raw::jsf64::walk_state(StateWalkingObject *walker) {
	//LOCKED, do not change
	//exception: changed in 0.87 to to reduce correlation between similar seeds
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
	walker->handle(d);
	if (!(a|b) && !(c|d)) d = 1;
	for (int i = 0; i < 12; i++) raw64();//12
}
void PractRand::RNGs::Raw::jsf64::self_test() {
	Raw::jsf64 rng;
	rng.seed(0ull);
	if (rng.raw64() != 0x76914495e6291d20ull) issue_error("jsf64::self_test() failed(a0)");
	if (rng.raw64() != 0x11596dd4917e4a2full) issue_error("jsf64::self_test() failed(a1)");
	if (rng.raw64() != 0x0d2ce75bc2869b29ull) issue_error("jsf64::self_test() failed(a2)");
	if (rng.raw64() != 0x0d066cc51c74176bull) issue_error("jsf64::self_test() failed(a3)");
	if (rng.raw64() != 0xd3672b8f73777390ull) issue_error("jsf64::self_test() failed(a4)");

	rng.seed(0x123456789a550895ull);
	if (rng.raw64() != 0xccd9a28ed7b4e456ull) issue_error("jsf64::self_test() failed(b0)");
	if (rng.raw64() != 0xa18b8f58aa157087ull) issue_error("jsf64::self_test() failed(b1)");
}



