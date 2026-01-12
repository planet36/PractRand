#include <string>
#if defined _MSC_VER
#include <intrin.h> // todo: remove
#endif
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>
#include "PractRand/endian.h"

#include "PractRand/RNGs/mrsf64.h"
#include "PractRand/RNGs/mrsf32.h"

using namespace PractRand;
using namespace PractRand::Internals;

//polymorphic:
PRACTRAND__POLYMORPHIC_RNG_BASICS_C64(mrsf64)
void PractRand::RNGs::Polymorphic::mrsf64::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
void PractRand::RNGs::Polymorphic::mrsf64::seed_fast(Uint64 s) { implementation.seed_fast(s); }
std::string PractRand::RNGs::Polymorphic::mrsf64::get_name() const { return "mrsf64"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C32(mrsf32)
void PractRand::RNGs::Polymorphic::mrsf32::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
void PractRand::RNGs::Polymorphic::mrsf32::seed_fast(Uint64 s) { implementation.seed_fast(s); }
std::string PractRand::RNGs::Polymorphic::mrsf32::get_name() const { return "mrsf32"; }

//#include <intrin.h>
/*
	basic algorithm:
		Word old = a;
		a = b * K;
		b = rotate(b, SHIFT) ^ old;
		return old + a;
	current shift (leftward rotation) values used:
		23 for 64 bit integers
		13 for 32 bit integers
		10 for 16 bit integers (16 bit version is not recommended, does not meet basic quality standards)
		undefined for 8 bit integers
	current multiplication constants used:
		0xAE3769b9D3519D65ull for 64 bit integers
		0xD3519D65ul for 32 bit integers
		0xB635 for 16 bit inegers

*/


//raw:
Uint64 PractRand::RNGs::Raw::mrsf64::raw64() {
	// this one is the currently mtsf64 algorithm, 32 bit variants are also decent
	Uint64 old = a;
	a = b * 0xAE3769b9D3519D65ull; // 1010111000110111011010011011100111010011010100011001110101100101
	b = rotate64(b, 23) ^ old;
	return old + a;
	/*
		output testing is simply not viable with 64 bit integers, so -ttseed64 only
		I'm going to go for a shift between 16 and 31, but will let seeding testing determine within that range, and will at least look at the results from outside of that range as sanity checks

		multiplier 0xAE3769b9D3519D65ull
		__		8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56
		s5		23	22	21	23	20	24	22	20	22	21	20	21	19	19	23	23	21	20	20	17	18	16	16	15	15	18	22	21	23	18	18	18	17	18	18	17	19	17	15	14	14
		s6		>38															>39				39	39	32	34	33	26

		Testing s6 is too slow for much.  Still, it's enough to show that the s6 results correlate with the s5 results.  
		So I'm going to go with 23 for my shift.  It's in a 3-way tie for 2nd place on the s5 results - 13 beats it, though not by much, and 22 and 11 match it.  But 13 feels like it's maybe on too small.  Same for 11.  And 22 has slightly worse neighbors on average.  So 23 it is.  
		And I've confirmed it doesn't have any problems on s6 either.  So at the very least, it's not bad.  

	*/


	/*
	// for testing a number of minor variants of the algorithm
	enum { ORDER = 3, SHIFT = 23 };
	const Uint64 K = 0xAE3769b9D3519D65ull;
	Uint64 old;
	if (ORDER == 0) {
		old = a;
		a = b * K;
		b = rotate64(b, SHIFT) ^ old;
		return old + a;
	}
	else if (ORDER == 1) {// SLOW
		old = a * K;
		a = b;
		b = rotate64(b, SHIFT) ^ old;
		return old + a;
	}
	else if (ORDER == 2) {
		old = a * K;
		a = rotate64(a, SHIFT) ^ b;
		b = old;
		return a+b;
	}
	else if (ORDER == 3) {
		old = a;
		a = rotate64(a, SHIFT) ^ b;
		b = old * K;
		return old + a;
	}//*/

}
void PractRand::RNGs::Raw::mrsf64::seed(Uint64 seed_low, Uint64 seed_high) {
	// zero is a forbidden state
	if (!seed_low && !seed_high) {
		seed_low = ~seed_low;
		seed_high = ~seed_high;
	}
	a = seed_low;
	b = seed_high;
	for (int i = 0; i < 8; i++) raw64();//8
}
void PractRand::RNGs::Raw::mrsf64::seed_fast(Uint64 s) {
	a = s;
	b = ~a;
	for (int i = 0; i < 6; i++) raw64();//6
}
void PractRand::RNGs::Raw::mrsf64::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	if (!a && !b) b = 1;// zero is a forbidden state
}



Uint32 PractRand::RNGs::Raw::mrsf32::raw32() {
	//*
	// this one is the current mtsf32 algorithm
	// moving the multiplication one line earlier seems to slow it down but improve output quality and seeding quality... despite the fact it's a completely equivalent state transition function (logically equivalent to an output function tweak)
	// this multiplier and shift outperform others I've checked, but it was not an exhaustive search - I took a brief glance at every possible shift, but only tried a few multipliers
	Uint32 old = a;
	a = b * 0xD3519D65ul;
	b = rotate32(b, 13) ^ old;
	return old + a;//*/

	/*


		K			4		5		6		7		8		9		10		11		12		13		14		15		16		17		18		19		20		21		22		23		24		25		26		27		28
		0x9B52A6C5	20/20	25/25	23/23	28/28	36/35	37/38	37/42	39/34																			?36?
			seed4	15/19	17/22	21/22	17/22	18/22	16/21	16/20	17/20	15/21	16/18	14/17	16/20	14/16	14/19	14/19	15/17	14/15	14/14	14/17	13/16	13/13	13/14	12/13	12/13	9/9
					+4		+5		+1		+5		+4		+5		+4		+3		+6		+2		+3		+4		+2		+5		+5		+2		+1		0		+3		+3		0		+1		+1		+1		0

		0xD3519D65	20/20	22/22	26/26	27/27	40/36	40/38	42/42	42/44			43/>43
			seed4	18/24	22/23	19/23	18/21	15/21	16/21	18/20	16/21	16/21	18/20	15/20	18/18	14/16	14/19	14/18	15/17	15/16	14/14	14/17	14/16	13/14	13/15	13/14	12/12	9/9


		delta		0/0		-3/-3	+3/+3	-1/-1	+4/+1	+3/0	+5/0	+3/+10
			d.sums			0/0						+6/0					+
			seed4	+3/+5	+5/+1	-2/+1	+1/-1	-3/-1	0/0		+2/0	-1/+1	+1/0	+2/+2	+1/+3	+2/-2	0/0		0/0		0/-1	0/0		+1/+1	0/0		0/0		+1/0	0/+1	0/+1	+1/+1	0/-1	0/0
			s4sums			+6/+7					-2/-2					+2/+1					+5/+3					0/-1					+1/+1					+1/+1					+1/+1

		0x187FF007	20/21	24/24	26/26	27/28	36/36	36/37	37/37

	*/
	/*enum {
		ORDER = 1
		, SHIFT = 13
		//,K = 0x9B52A6C5ul
		, K = 0xD3519D65ul
		//, K = 0x187FF007ul //a deliberately bad multiplier - 00011000011111111111000000000111
	};
	if (ORDER == 0) {
		Uint32 old = a;
		a = b * K;
		b = rotate32(b, SHIFT) ^ old;
		return old + a;
	}
	else if (ORDER == 1) {
		Uint32 old = a * K;
		a = b;
		b = rotate32(b, SHIFT) ^ old;
		return old + a;
	}
	else if (ORDER == 2) {
		Uint32 old = a * K;
		a = rotate32(a, SHIFT) ^ b;
		b = old;
		return a+b;
	}
	else if (ORDER == 3) {
		Uint32 old = a;
		a = rotate32(a, SHIFT) ^ b;
		b = old * K;
		return a + old;
	}//*/
}
void PractRand::RNGs::Raw::mrsf32::seed(Uint64 seed_low, Uint64 seed_high) {
	if (!seed_low) seed_low = ~seed_low;// zero is a forbidden state for the current version of the algorithm
	a = Uint32(seed_low);
	b = Uint32(seed_low >> 32);
	for (int i = 0; i < 8; i++) raw32();//8
}
void PractRand::RNGs::Raw::mrsf32::seed_fast(Uint64 s) {
	if (!s) s = ~s;
	a = Uint32(s);
	b = Uint32(s >> 32);
	for (int i = 0; i < 6; i++) raw32();//6
}
void PractRand::RNGs::Raw::mrsf32::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	if (!a && !b) a = b = ~0ull;// zero is a forbidden state for the current version of the algorithm
}


