#include <string>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>
#include "PractRand/endian.h"

#include "PractRand/RNGs/mrc64.h"
#include "PractRand/RNGs/mrc32.h"
#include "PractRand/RNGs/mrc16.h"

using namespace PractRand;
using namespace PractRand::Internals;

//polymorphic:
PRACTRAND__POLYMORPHIC_RNG_BASICS_C16(mrc16)
void PractRand::RNGs::Polymorphic::mrc16::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::mrc16::get_name() const { return "mrc16"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C32(mrc32)
void PractRand::RNGs::Polymorphic::mrc32::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::mrc32::get_name() const { return "mrc32"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C64(mrc64)
void PractRand::RNGs::Polymorphic::mrc64::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::mrc64::get_name() const { return "mrc64"; }

//raw:
Uint64 PractRand::RNGs::Raw::mrc64::raw64() {
	Uint64 old;
	enum { SHIFT = 21 };// 21
	const Uint64 K = 0x9e3779b97f4a7c15ull;

	//*
	old = a * K; // this one is currently the official mrc64 algorithm, usable at 16, 32, and 64 bits
	a = b + counter++;
	b = rotate64(b, SHIFT) ^ old;
	return old + a;//*/

	/*
	old = a + counter++;
	a = b * K;
	b = rotate64(b, SHIFT) ^ old;
	return old + a;//*/

	/*
	old = a * K + counter++;
	a = rotate64(b, 11) + b;
	b = rotate64(b, SHIFT) ^ old;
	return old + a;//*/

	/*
	old = a + b;
	a = (b * K) ^ counter++;
	b = rotate(old, SHIFT);
	return a;//*/


	/*
	old = a;
	//old = a + counter++;
	//old = a + b;
	//old = a + b + counter++;
	a = (b * K) ^ rotate(a, SHIFT);
	b = old;
	return old;//*/

	/*
	old = a + b;
	a = b * 0x9e3779b97f4a7c15ull;
	b = rotate64(old, 23) + ++counter;
	return a;//*/
	/*
	old = a * 0x9e3779b97f4a7c15ull;
	a = rotate(a, SHIFT) * 9;
	b += old;
	a ^= b;
	return b;//*/
}
void PractRand::RNGs::Raw::mrc64::seed(Uint64 seed_low, Uint64 seed_high) {
	a = seed_low;
	b = seed_high;
	counter = 1;
	for (int i = 0; i < 6; i++) raw64();
}
void PractRand::RNGs::Raw::mrc64::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	walker->handle(counter);
}


Uint32 PractRand::RNGs::Raw::mrc32::raw32() {
	Uint32 old;
	enum { SHIFT = 19 };
	const Uint32 K = 0x7f4a7c15;

	old = a * K; // this one is currently the official mrc32 algorithm
	a = b + counter++;
	b = rotate32(b, SHIFT) ^ old;
	return old + a;//*/

	//old = a;
	//old = a + counter++;
	//old = a + b;
	/*old = a + b + counter++;
	a = (b * K) ^ rotate(a, SHIFT);
	b = old;
	return old;*/

	//three basic variations:

	// 1. fastest, passes output tests but has horrible avalanche properties
	/*
	Uint32 old = a + b;
	a = (b * 0x7f4a7c15) ^ counter++;
	b = rotate32(old, 19);
	return a;//*/

	// 2. much better avalanche properties, but is irreversible so cycle length is just ~2^48
	/*
	Uint32 old = a + b + counter++;
	a = (b * 0x7f4a7c15) ^ rotate32(a, 12);
	b = rotate32(old, 19);
	return a;//*/

	// 3. now reversible, but... notice that returning 'a' instead of 'old' makes it fail binary rank tests?  I'm not sure I trust it.  
	/*
	Uint32 old = a + counter++;
	a = (b * 0x7f4a7c15) ^ rotate32(a, 24);// output tests: 8:19, 9:19, 10:37, 11:37-, 12:37, 13:37, 14:37, 15:37, 16:37, 19:37, 21:36, 23:36, 24:37, 25:28, 26:29
	b = old;// seeding tests (6 discarded): 8:24, 9:27, 10:27, 11:26, 12:26, 13:27, 14:24, 15:21, 16:9, 17:19, 18:18, 19:22, 20:22, 21:20, 22:23, 23:20, 24:17
	return old;//*/
}
void PractRand::RNGs::Raw::mrc32::seed(Uint64 seed_low, Uint64 seed_high) {
	a = seed_low;
	b = seed_low >> 32;
	counter = 1;
	for (int i = 0; i < 7; i++) raw32();
}
void PractRand::RNGs::Raw::mrc32::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	walker->handle(counter);
}


Uint16 PractRand::RNGs::Raw::mrc16::raw16() {
	Uint16 old;
	enum { SHIFT = 10 };
	const Uint16 K = 0xA965ull; // 0xA965ull; // aka 1010100101100101
	//const Uint16 K = 0xB4CDull; // 1011010011001101

	//testing different constants:
	//0xA965ull multiplier
	//	_	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15		the shift constants
	//	s6	17	21	22	22	23	20	20	20	20	18	18	14	17	13	12		
	//	s7	20	25	29	27	26	24	27	23	26	28	23	17	17	14	12		log2 of the number of seeds necessary to detect interseed correlation, when seeding skips 7 outputs
	//	s8	21	26	31	32	32	35	34	26	31	33	25	28	18	17	13		same but skipping 8 ; still not enough, but getting close
	//	s9	20	26	30	32	35	41	37	36	>41	>40	33	>34	21	17	13		some of those individually took multiple days to calculate, so I got fed up and started aborting them earlier
	//	s14	21	26	30	33	35									34	18		skipped ahead a bit here, to investigate assymptotic behavior of bad shifts
	//	out	14	18	20	33	36	42	37	30	34	40	39	39	23	18	16
	// seeding performance scales WAY too poorly for smaller shifts, as if this wasn't actually a chaotic system, but I don't see why any non-zero shift would be non-chaotic
	//0xB4CDull multiplier
	//	_	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15		the shift constants
	//	s6	17	20	22	22	21	20	20	18	
	//	s7	19	24	31	26	31	25	28	24	
	//	s8	20	27	32
	//	out	14	18	20	33	35	38	33	22	28	41	37	38	23	18	16

	old = a * K; // this one is currently the official mrc16 algorithm, usable at 16, 32, and 64 bits
	a = b + counter++;
	b = rotate16(b, SHIFT) ^ old;
	return old + a;//*/
}
void PractRand::RNGs::Raw::mrc16::seed(Uint64 seed_low, Uint64 seed_high) {
	a = seed_low;
	b = seed_low >> 16;
	counter = seed_low >> 32;
	for (int i = 0; i < 9; i++) raw16();//12 - I'd say 14, but it already performs better on seeding testing than it does on output testing, so further increases are likely wasted
}
void PractRand::RNGs::Raw::mrc16::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	walker->handle(counter);
}



