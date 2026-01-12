#include <cstring>
#include <string>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>

#include "PractRand/RNGs/efiix8min.h"
#include "PractRand/RNGs/efiix8x48.h"
#include "PractRand/RNGs/efiix16x48.h"
#include "PractRand/RNGs/efiix32x48.h"
#include "PractRand/RNGs/efiix64x48.h"

#include "PractRand/RNGs/arbee.h"
//#include "PractRand/RNGs/trivium.h"

using namespace PractRand;
using namespace PractRand::Internals;

PRACTRAND__POLYMORPHIC_RNG_BASICS_C8(efiix8min)
void PractRand::RNGs::Polymorphic::efiix8min::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::efiix8min::get_name() const { return "efiix8min"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C8(efiix8x48)
void PractRand::RNGs::Polymorphic::efiix8x48::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::efiix8x48::get_name() const { return "efiix8x48"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C16(efiix16x48)
void PractRand::RNGs::Polymorphic::efiix16x48::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::efiix16x48::get_name() const { return "efiix16x48"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C32(efiix32x48)
void PractRand::RNGs::Polymorphic::efiix32x48::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::efiix32x48::get_name() const { return "efiix32x48"; }
PRACTRAND__POLYMORPHIC_RNG_BASICS_C64(efiix64x48)
void PractRand::RNGs::Polymorphic::efiix64x48::seed(Uint64 seed_low, Uint64 seed_high) { implementation.seed(seed_low, seed_high); }
std::string PractRand::RNGs::Polymorphic::efiix64x48::get_name() const { return "efiix64x48"; }

//#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word old = a + i; \
	Word iterated = iteration_table[i % ITERATION_SIZE]; \
	Word indirect = indirection_table[c % INDIRECTION_SIZE]; \
	indirection_table[c % INDIRECTION_SIZE] = iterated; \
	iteration_table[i++ % ITERATION_SIZE] = old; \
	a = b + c; \
	b = c ^ rotate ## BITS(iterated, BITS / 2); \
	c = indirect + rotate ## BITS(old, SHIFT_AMOUNT); \
	return iterated + indirect;


//#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word old = rotate ## BITS(a, 2) ^ b; \
	Word iterated = iteration_table[i % ITERATION_SIZE]; \
	Word indirect = indirection_table[c % INDIRECTION_SIZE]; \
	indirection_table[c % INDIRECTION_SIZE] = iterated; \
	iteration_table[i % ITERATION_SIZE] = old; \
	a = b + i++; \
	b = c + old; \
	c = indirect + rotate ## BITS(c, SHIFT_AMOUNT); \
	return iterated + indirect;

//#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word iterated = iteration_table[i % ITERATION_SIZE]; \
	Word indirect = indirection_table[c % INDIRECTION_SIZE]; \
	indirection_table[c % INDIRECTION_SIZE] = iterated + a; \
	iteration_table[i % ITERATION_SIZE] = rotate ## BITS(indirect, b % BITS); \
	Word old = a ^ b; \
	a = b + i++; \
	b = c + indirect; \
	c = old + rotate ## BITS(c, SHIFT_AMOUNT); \
	return b ^ iterated;

#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word iterated = iteration_table[i % ITERATION_SIZE]; \
	Word indirect = indirection_table[c % INDIRECTION_SIZE]; \
	indirection_table[c % INDIRECTION_SIZE] = iterated + a; \
	iteration_table[i % ITERATION_SIZE] = indirect; \
	Word old = a ^ b; \
	a = b + i++; \
	b = c + indirect; \
	c = old + rotate ## BITS(c, SHIFT_AMOUNT); \
	return b ^ iterated;
//current algorithm
//perfect statistically (unless ITERATION_SIZE and BITS are both very small), reasonably fast
//safe from a variety of attacks
//but the "b ^ iterated" and "iterated + a" have too much in common, and the later can come right back out of the indirection pool right away
//not much weakness, but some
//alternatively, can return "b" alone - it's horrible at 8/16 bit, but acceptable at 32 and good at 64

//later... It all looks solid to me right now
//in particular, I said there was too much in common between "b ^ iterated" and "iterated + a", but I now disagree with my former self
//"iterated + a" uses the old "a", while in terms of the old values the other is more "iterated ^ (c + indirect)".  
//whereas if you try to compare to the next pass or previous pass, "iterated" has a largely unrelated value
//curently my view is that the whole thing is good, specifically because all of the folowing requirements are met:
//	the flow of information to output is one word per call
//	the flow of information from the mixing pool to the tables is at least one word per call - arguablly more, because "c" is used to index in to the indirection table
//	the flow of information from the tables to the mixing pool is one word per call
//	the flow of information from the mixing pool to output is one word per call, and largely independent of the flow to the mixing pool
//	the flow of information from the tables to output is one word per call
//	passes statistical tests even with table sizes are very small
//using only output, once the above requirements are met, I think the only way it breaks is if either there's a very obvious and simple flaw, or if the mixing isn't actually chaotic, both of which would have good odds of showing up in thorough testing
//	...I suppose the mixing being too slow could also be a problem?  but since it does fine with table sizes set to 1, that's probably okay?
//to be clear, I don't think those are the only possible set of sufficient requirements - for instance, I think output could be taken entirely from the mixing pool or entirely from the tables IFF both the flow from mixing to tables and the flow from tables to mixing were each clearly wider than the output rate

//#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word iterated = iteration_table  [i % ITERATION_SIZE];\
	Word indirect = indirection_table[c % INDIRECTION_SIZE];\
	indirection_table[c % INDIRECTION_SIZE] = iterated ^ a;\
	iteration_table  [i % ITERATION_SIZE  ] = indirect;\
	Word old = a + i++;\
	a = b + iterated;\
	b = c ^ indirect;\
	c = old + rotate ## BITS (c, SHIFT_AMOUNT);\
	return b;
//adequate statistically, good speed
//wider connection between mixing pool (a,b,c) and indirection pool (the two arrays) means it's now safe to use a simpler output function
//...I think, anyway
//"old" now being a function of (a,i) instead of (a,b) means worse mixing in the pool, I think


//#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word iterated = iteration_table  [i % ITERATION_SIZE] ^ i;\
	Word indirect = indirection_table[c % INDIRECTION_SIZE] + a;\
	indirection_table[c % INDIRECTION_SIZE] = iterated;\
	iteration_table  [i % ITERATION_SIZE  ] = indirect;\
	Word old = a^b;\
	a = b + indirect;\
	b = c + iterated;\
	c = old + rotate ## BITS (c, SHIFT_AMOUNT);\
	i++;\
	return b;
//not good enough statistically (8x1+1: 128 GB, 8x2+2: 32 MB, 8x4+4: 8 GB, 8x8+8: 32 GB), adequate speed
//looks like the above one but better mixing in the pool?


//#define EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT) \
	Word iterated = iteration_table  [i % ITERATION_SIZE] ^ a;\
	Word indirect = indirection_table[c % INDIRECTION_SIZE] + i;\
	indirection_table[c % INDIRECTION_SIZE] = iterated;\
	iteration_table  [i % ITERATION_SIZE  ] = indirect;\
	Word old = a + b;\
	a = b + iterated;\
	b = c + indirect;\
	c = old ^ rotate ## BITS (c, SHIFT_AMOUNT);\
	i++; return old;
//decent statistically, decent speed
//wider connection between mixing pool (a,b,c) and indirection pool (the two arrays) means it's now safe to use a simpler output function
//...I think, anyway
//but directly returning "old" means that if they guess (a,b,c) correctly once, they can figure out all other (a,b,c) sets from looking at the output

#define _EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, ITER, IND, DO_RET) \
	Word iterated = ITER;\
	Word indirect = IND;\
	IND = iterated + a;\
	ITER = indirect;\
	Word old = a ^ b; \
	a = b + i++; \
	b = c + indirect; \
	c = old + rotate ## BITS(c, SHIFT_AMOUNT); \
	if (DO_RET) return b ^ iterated;


/*
	RNG             Quality     PractRand   TestU01     RaBiGeTe    Dieharder  Dieharder
	                subscores   standard    *Crush      Extended*1  -a         custom
*	efiix4x(2+2)     0/0/0/0    512 MB      2/?/?       pass(b)
*	efiix4x(1+4)     0/0/0/0    128 MB      3/?/?       128 Mb
*	efiix4x(4+1)     1/0/0/0    8 GB        0/0/12      pass
*	efiix4x(2+4)     0/0/1/0    4 GB        2/?/?       pass?
*	efiix4x(4+2)     2/0/1/0    16 GB       pass        pass?
*	efiix4x(4+4)     1/1/1/1    256 GB      0/1/2       pass?
*	efiix4x(1+8)     1/1/1/1    1 GB        0/14/?      512 Mb
*	efiix4x(8+1)     5/1/1/1    > 1 TB      > 1 TB      pass        pass?
*	efiix4x(2+8)     ?/1/1/1    ?           ?           0/13/?      pass?
*	efiix4x(1+16)    ?/3/2/1    ?           ?           0/1/6       pass
---
*	efiix8x(1+1)     4/1/1/1    > 8 TB      > 4 TB      pass        pass?
*	efiix8x(1+2)     5/1/1/1    > 2 TB      > 1 TB      pass        pass?
*	efiix8x(2+1)     ?/1/1/1                
*	efiix8x(2+2)     ?/2/2/2    

*	efiix4x(2+2)     0/0/0/0    2 GB*       256 MB      2/?/?       pass
*	efiix4x(1+4)     0/0/0/0    1 GB*       64 MB       3/?/?       128 Mb
*	efiix4x(4+1)     1/0/0/0    16 GB(pe)*  4 GB        0/0/12      pass
*	efiix4x(2+4)     0/0/1/0    16 GB*      4 GB        2/?/?       pass?
*	efiix4x(4+2)     2/0/1/0    256 GB*     16 GB       pass        pass?
*	efiix4x(4+4)     1/1/1/1    2 TB*       256 GB      0/1/2       pass?
*	efiix8x(1+1)     4/1/1/1    > 8 TB*     > 4 TB      pass        pass?


It looks like quality is suffering with 4 bit words... 
	one possibly is how few shift constants are possible... maybe they all suck

Some oddities related to empirical quality vs parameterization visible there at 4 bit
	impossible to tell if they remain at 8 bit


*/


/*	on the algorithm:
		examining small windows of the output is never sufficient to learn anything significant, 
			due to relatively large FIFO "iteration_table"
		examining small windows at a stride of ITERATION_SIZE apart might be sufficient, but...
			if the indirection pattern is not predicted/guessed that is meaningless
			if the indirection pattern is predicted/guessed, 
				if it doesn't reach back ITERATION_SIZE elements
					that will simply add more unknowns
				if it does reach back ITERATION_SIZE elements
					you still need to establish the contents of (a,b,c)
					and besides that's an enormous amount of indirection pattern to figure out
		figuring out any significant portion of the indirection pattern should be very very difficult
	consider:
		if we guess (a,b,c,i) and "indirect" correctly (a LOT of correct guesses)
		then we can use all that to calculate "iterated" from the output
		and then use all that plus a guess that the indirection index is identical
			(which isn't exactly a guess since that's a function of other values we guessed, but it doesn't happen all that often)
		to figure out all the same information at the next position
		it might look like we have it almost cracked - we can predict important future (partial) states with relatively high probabilities
		but that's useless if we can't TELL if our guesses were right
		1. 
			which we can't do from a short output window
				because we need the output to figure out "iterated" values feeding in
			so we keep guessing that the indirection index remains constant, ITERATIONS_SIZE times in a row
			THEN we can start to confirm our guesses
			but that's already on the order of 2**(W*5+ITERATION_SIZE*L) operations to get a success
			W is word size, in bits; L is log2(INDIRECTION_SIZE)
		2. 
			which we can't do from a single short output window
			...so we have to skip forward roughly ITERATION_SIZE positions to a 2nd window
			and do the whole thing over again (but we don't have to guess "i" again)
			and use the old predicted values of "indirect" fed in to the iteration table
			to check the new predicted values of "iterated"
			and then we can START confirming our guesses
			at region1 we guessed (a,b,c) and "i" and "indirect"
			at region2 we guessed (a,b,c) and "indirect"
		3. 
			we check if the "iterated" values observed actually match the constant-indirection-index required
			but even if they appear to, that doesn't do us much good, as it's still far more likely to be a false positive than a true positive
			unless they appear to for a very long sequence... but the chances of such a region even existing in a cyphertext length <2**64 is very small
	overall, I think the 32 bit word variant ought to be good for 128 bit security, and the 64 bit word variant ought to be good for 160 bit security
	I don't have much confidence in the 8 or 16 bit variants though
*/

#define EFIIX_SEED_SIMPLE_128(BITS, seed_low, seed_high) \
	for (unsigned long x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = x; \
	for (unsigned long x = 0; x < ITERATION_SIZE; x++) iteration_table[x] = x; \
	a = Uint ## BITS(seed_low); \
	b = Uint ## BITS(seed_high); \
	c = 0; i = 1; \
	for (unsigned long x = 1; x < 64 / BITS; x++) { \
		for (unsigned long y = 0; y < 16; y++) raw ## BITS(); \
		a ^= Uint ## BITS(seed_low >> (BITS * x)); \
		b ^= Uint ## BITS(seed_high >> (BITS * x)); \
	} \
	for (unsigned long x = 0; x < ITERATION_SIZE + INDIRECTION_SIZE + 16; x++) raw ## BITS();

#define EFIIX_SEED_SIMPLE(BITS, _a, _b, _c) \
	for (unsigned long x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = x; \
	for (unsigned long x = 0; x < ITERATION_SIZE; x++) iteration_table[x] = x; \
	a = _a; b = _b; c = _c; i = 1; \
	for (unsigned long x = 0; x < ITERATION_SIZE + INDIRECTION_SIZE + 16; x++) raw ## BITS();

#define EFIIX_SEED_VECTOR_RAW(BITS, source, length) \
	a = 0xDEADBEEF; b = 0; c = a; i = 0; \
	for (unsigned long x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = x; \
	if (length < ITERATION_SIZE / 2) {\
		for (unsigned long x = 0; x < length; x++) iteration_table[x] = source[x]; \
		for (unsigned long x = length; x < ITERATION_SIZE / 2; x++) iteration_table[x] = x; \
		for (unsigned long x = 0; x < length; x++) iteration_table[x + ITERATION_SIZE / 2] = source[x]; \
		for (unsigned long x = length + ITERATION_SIZE / 2; x < ITERATION_SIZE; x++) iteration_table[x] = x; \
		for (unsigned long x = 0; x < length * 2; x++) raw ## BITS(); \
	}\
	else {\
		for (unsigned long x = 0; x < ITERATION_SIZE / 2; x++) iteration_table[x] = source[x]; \
		for (unsigned long x = 0; x < ITERATION_SIZE / 2; x++) iteration_table[x + ITERATION_SIZE / 2] = source[x]; \
		for (unsigned long x = 0; x < ITERATION_SIZE; x++) raw ## BITS(); \
		length -= ITERATION_SIZE / 2;\
		source += ITERATION_SIZE / 2;\
		while (length >= ITERATION_SIZE / 2) {\
			for (unsigned long x = 0; x < ITERATION_SIZE / 2; x++) iteration_table[x] ^= source[x]; \
			for (unsigned long x = 0; x < ITERATION_SIZE / 2; x++) iteration_table[x + ITERATION_SIZE / 2] ^= source[x]; \
			for (unsigned long x = 0; x < ITERATION_SIZE; x++) raw ## BITS(); \
			length -= ITERATION_SIZE / 2;\
			source += ITERATION_SIZE / 2;\
		}\
		if (length) {\
			for (unsigned long x = 0; x < length; x++) iteration_table[x] ^= source[x]; \
			for (unsigned long x = 0; x < length; x++) iteration_table[x + ITERATION_SIZE / 2] ^= source[x]; \
			for (unsigned long x = 0; x < length << 1; x++) raw ## BITS(); \
		}\
	}\
	for (unsigned long x = 0; x < ITERATION_SIZE + 16; x++) raw ## BITS();

#define EFIIX_SEED_VECTOR_PUMP64(BITS, _seed_vector, _seed_length) \
	a = b = c = i = 1;\
	for (int w = 0; w < ITERATION_SIZE; w++) iteration_table[w] = w;\
	for (int w = 0; w < INDIRECTION_SIZE; w++) indirection_table[w] = w;\
	for (int w = 0; w < _seed_length; w++) {\
		if (BITS == 64) { \
			a ^= _seed_vector[w]; \
			raw ## BITS(); raw ## BITS(); \
		}\
		else for (int bp = 0; bp < 64; bp += 2*BITS) {\
			a ^= Word(_seed_vector[w] >> bp); \
			b ^= Word(_seed_vector[w] >> (bp + BITS)); \
			raw ## BITS(); raw ## BITS(); raw ## BITS(); \
		}\
		a ^= _seed_length;\
	}\
	for (int w = 0; w < ITERATION_SIZE; w++) raw8();\
	for (int w = 0; w < INDIRECTION_SIZE * 2; w++) raw8();

#define EFIIX_SEED_BRUTE_FORCE(BITS, _s1, _s2) \
	{PractRand::RNGs::Raw::arbee seeder; \
	Uint64 s1 = _s1, s2 = _s2;\
	seeder.seed(s1, s2); \
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) indirection_table[w] = Word(seeder.raw64()); \
	i = Word(seeder.raw64()); \
	for (unsigned long w = 0; w < ITERATION_SIZE; w++) iteration_table[(w + i) % ITERATION_SIZE] = Word(seeder.raw64()); \
	a = Word(seeder.raw64()); \
	b = Word(seeder.raw64()); \
	c = Word(seeder.raw64()); \
	for (unsigned long w = 0; w < 64; w++) raw ## BITS(); \
	seeder.raw64(); s1 += seeder.raw64(); s2 += seeder.raw64(); \
	s1 ^= a; s2 ^= b; \
	if (BITS < 32) {raw ## BITS(); raw ## BITS(); s1 = rotate64(s1, BITS) ^ a; s2 = rotate64(s2, BITS) ^ b;} \
	if (BITS < 16) {raw ## BITS(); raw ## BITS(); s1 = rotate64(s1, BITS) ^ a; s2 = rotate64(s2, BITS) ^ b;} \
	raw ## BITS(); raw ## BITS(); \
	seeder.seed(s1, s2); \
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) indirection_table[w] ^= Word(seeder.raw64()); \
	a ^= Word(seeder.raw64()); \
	b ^= Word(seeder.raw64()); \
	c ^= Word(seeder.raw64()); \
	for (unsigned long w = 0; w < ITERATION_SIZE + 16; w++) raw ## BITS(); }

#define EFIIX_SEED_FOLDING(BITS, SHIFT_AMOUNT, _s1, _s2, _s3) \
	a = _s1; b = _s2; c = _s3;\
	for (int w = 0; w < INDIRECTION_SIZE; w++) indirection_table[w] = w; \
	for (unsigned long w = 0; w < 1 + 2; w++) { EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, 1, INDIRECTION_SIZE) } \
	if (ITERATION_SIZE >= 2) { for (unsigned long w = 0; w < 1 + 2; w++) { EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, 1, INDIRECTION_SIZE) } } \
	if (ITERATION_SIZE >= 4) { for (unsigned long w = 0; w < 1 + 2; w++) { EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, 4, INDIRECTION_SIZE) } } \
	if (ITERATION_SIZE >= 8) { for (unsigned long w = 0; w < 1 + 2; w++) { EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, 8, INDIRECTION_SIZE) } } \
	if (ITERATION_SIZE >= 16) { for (unsigned long w = 0; w < 1 + 2; w++) { EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, 16, INDIRECTION_SIZE) } } \
	if (ITERATION_SIZE >= 32) { for (unsigned long w = 0; w < 1 + 2; w++) { EFIIX_ALGORITHM(BITS, SHIFT_AMOUNT, 32, INDIRECTION_SIZE) } } \
	;

Uint8 PractRand::RNGs::Raw::efiix8min::raw8() {
	//_EFIIX_ALGORITHM(8, 3, iteration_table[i%ITERATION_SIZE], indirection_table[c%INDIRECTION_SIZE], 1)
	_EFIIX_ALGORITHM(8, 3, iteration_table, indirection_table, 1)
}
void PractRand::RNGs::Raw::efiix8min::seed(Uint64 seed_low, Uint64 seed_high) {
	// only the lowest 48 bits of the 128 bit seed are used
	a = Uint8(seed_low >> 0);
	b = Uint8(seed_low >> 8);
	c = Uint8(seed_low >> 16);
	i = Uint8(seed_low >> 24);
	iteration_table = Uint8(seed_low >> 32);
	indirection_table = Uint8(seed_low >> 40);
	for (int x = 0; x < 12; x++) raw8();
	/*
		redo tests
	*/
}
void PractRand::RNGs::Raw::efiix8min::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
	walker->handle(i);
	walker->handle(iteration_table);
	walker->handle(indirection_table);

}


PractRand::RNGs::Raw::efiix8x48::~efiix8x48() { std::memset(this, 0, sizeof(*this)); }
PractRand::RNGs::Raw::efiix16x48::~efiix16x48() { std::memset(this, 0, sizeof(*this)); }
PractRand::RNGs::Raw::efiix32x48::~efiix32x48() { std::memset(this, 0, sizeof(*this)); }
PractRand::RNGs::Raw::efiix64x48::~efiix64x48() { std::memset(this, 0, sizeof(*this)); }

Uint8 PractRand::RNGs::Raw::efiix8x48::raw8() {
	typedef Word check_efiix_array_sizes[(ITERATION_SIZE & (ITERATION_SIZE-1)) || (INDIRECTION_SIZE & (INDIRECTION_SIZE-1)) ? -1 : 1];
	EFIIX_ALGORITHM(8, 3)
}
void PractRand::RNGs::Raw::efiix8x48::add_entropy8(Uint8 value) {
	a += value;
	raw8();
}
void PractRand::RNGs::Raw::efiix8x48::seed(Uint64 seed_low, Uint64 seed_high) {
	EFIIX_SEED_SIMPLE_128(8, seed_low, seed_high)
}
/*void PractRand::RNGs::Raw::efiix8x48::seed(PractRand::RNGs::vRNG *source_rng) {
	a = source_rng->raw8();
	b = source_rng->raw8();
	c = source_rng->raw8();
	i = source_rng->raw8();
	// number of possible seeded states is kept much smaller than the number of valid states
	// in order to make it extremely unlikely that any bad seeds exist
	// as opposed to how it would be otherwise, in which finding a bad seed would be *extremely* difficult, but a few would likely exist
	for (int x = 0; x < ITERATION_SIZE; x++) if (!(x & 3)) iteration_table[x] = source_rng->raw8(); else iteration_table[x] = iteration_table[x & ~3];
	for (int x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = source_rng->raw8();

	//ought to be a secure seeding if source_rng->get_flags() includes FLAG::CRYPTOGRAPHIC_SECURITY, so...
	for (int x = 0; x < ITERATION_SIZE; x++) raw8();
}*/
void PractRand::RNGs::Raw::efiix8x48::walk_state(StateWalkingObject *walker) {
	for (unsigned long w = 0; w < ITERATION_SIZE; w++) walker->handle(iteration_table[w]);
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) walker->handle(indirection_table[w]);
	walker->handle(i);
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
}


Uint16 PractRand::RNGs::Raw::efiix16x48::raw16() {
	typedef Word check_efiix_array_sizes[(ITERATION_SIZE & (ITERATION_SIZE-1)) || (INDIRECTION_SIZE & (INDIRECTION_SIZE-1)) ? -1 : 1];//making sure they're powers of 2
	EFIIX_ALGORITHM(16, 7)
}
void PractRand::RNGs::Raw::efiix16x48::seed(Uint64 seed_low, Uint64 seed_high) {
	EFIIX_SEED_SIMPLE_128(16, seed_low, seed_high)
	//EFIIX_SEED_BRUTE_FORCE(16, _s1, _s2)
}
/*void PractRand::RNGs::Raw::efiix16x48::seed(PractRand::RNGs::vRNG *source_rng) {
	a = source_rng->raw16();
	b = source_rng->raw16();
	c = source_rng->raw16();
	i = source_rng->raw16();
	// number of possible seeded states is kept much smaller than the number of valid states
	// in order to make it extremely unlikely that any bad seeds exist
	// as opposed to how it would be otherwise, in which finding a bad seed would be *extremely* difficult, but a few would likely exist
	for (int x = 0; x < ITERATION_SIZE; x++) if (!(x & 3)) iteration_table[x] = source_rng->raw16(); else iteration_table[x] = iteration_table[x & ~3];
	for (int x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = source_rng->raw16();

	//ought to be a secure seeding if source_rng->get_flags() includes FLAG::CRYPTOGRAPHIC_SECURITY, so...
	for (int x = 0; x < ITERATION_SIZE; x++) raw16();
}*/
void PractRand::RNGs::Raw::efiix16x48::walk_state(StateWalkingObject *walker) {
	for (unsigned long w = 0; w < ITERATION_SIZE; w++) walker->handle(iteration_table[w]);
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) walker->handle(indirection_table[w]);
	walker->handle(i);
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
}


Uint32 PractRand::RNGs::Raw::efiix32x48::raw32() {
	typedef Word check_efiix_array_sizes[(ITERATION_SIZE & (ITERATION_SIZE-1)) || (INDIRECTION_SIZE & (INDIRECTION_SIZE-1)) ? -1 : 1];
	EFIIX_ALGORITHM(32, 13)
}
static void mix4x32(Uint32 &a, Uint32 &b, Uint32 &c, Uint32 &d) {
	b ^= a + (a << 13);
	c += b ^ (b >> 5);
	d ^= c + b;
	a += d ^ ((c << 8) | (c >> 24));
	b ^= a + (a << 11);
	c += b ^ (b >> 9);
	d ^= c + b;
	a += d ^ ((c << 8) | (c >> 24));
}
void PractRand::RNGs::Raw::efiix32x48::seed(Uint64 seed_low, Uint64 seed_high) {
	EFIIX_SEED_SIMPLE_128(32, seed_low, seed_high)
	//EFIIX_SEED_BRUTE_FORCE(32, _s1, _s2)
}
/*void PractRand::RNGs::Raw::efiix32x48::seed(PractRand::RNGs::vRNG *source_rng) {
	a = source_rng->raw32();
	b = source_rng->raw32();
	c = source_rng->raw32();
	i = source_rng->raw32();
	// number of possible seeded states is kept much smaller than the number of valid states
	// in order to make it extremely unlikely that any bad seeds exist
	// as opposed to how it would be otherwise, in which finding a bad seed would be *extremely* difficult, but a few would likely exist
	for (int x = 0; x < ITERATION_SIZE; x++) if (!(x & 3)) iteration_table[x] = source_rng->raw32(); else iteration_table[x] = iteration_table[x & ~3];
	for (int x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = source_rng->raw32();

	//ought to be a secure seeding if source_rng->get_flags() includes FLAG::CRYPTOGRAPHIC_SECURITY, so...
	for (int x = 0; x < ITERATION_SIZE; x++) raw32();
}*/
//void PractRand::RNGs::Raw::efiix32x48::seed(const Word *seeds, int num_seeds, int seeding_quality) {
//}
void PractRand::RNGs::Raw::efiix32x48::walk_state(StateWalkingObject *walker) {
	for (unsigned long w = 0; w < ITERATION_SIZE; w++) walker->handle(iteration_table[w]);
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) walker->handle(indirection_table[w]);
	walker->handle(i);
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
}


Uint64 PractRand::RNGs::Raw::efiix64x48::raw64() {
	typedef Word check_efiix_array_sizes[(ITERATION_SIZE & (ITERATION_SIZE-1)) || (INDIRECTION_SIZE & (INDIRECTION_SIZE-1)) ? -1 : 1];
	EFIIX_ALGORITHM(64, 25)
}
void PractRand::RNGs::Raw::efiix64x48::seed(Uint64 seed_low, Uint64 seed_high) {
	EFIIX_SEED_SIMPLE_128(64, seed_low, seed_high)
	//EFIIX_SEED_BRUTE_FORCE(64, _s1, _s2)

	/*PractRand::RNGs::Raw::arbee seeder;
	Uint64 s1 = _s1, s2 = _s2, s3 = _s3, s4 = _s4;
	seeder.seed(s1, s2, s3, s4);
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) indirection_table[w] = Word(seeder.raw64());
		i = Word(seeder.raw64());
	for (unsigned long w = 0; w < ITERATION_SIZE; w++) iteration_table[(w + i) % ITERATION_SIZE] = Word(seeder.raw64());
		a = Word(seeder.raw64());
		b = Word(seeder.raw64());
		c = Word(seeder.raw64());
	for (unsigned long w = 0; w < 64; w++) raw64();//
	seeder.raw64(); s1 += seeder.raw64(); s2 += seeder.raw64(); s3 += seeder.raw64();
	seeder.seed(s1^a, s2^b, s3^c, ~s4);
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) indirection_table[w] ^= Word(seeder.raw64());
	for (unsigned long w = 0; w < ITERATION_SIZE + 16; w++) raw64();//
	*/
}
/*void PractRand::RNGs::Raw::efiix64x48::seed(PractRand::RNGs::vRNG *source_rng) {
	a = source_rng->raw32();
	b = source_rng->raw32();
	c = source_rng->raw32();
	i = source_rng->raw32();
	// number of possible seeded states is kept much smaller than the number of valid states
	// in order to make it extremely unlikely that any bad seeds exist
	// as opposed to how it would be otherwise, in which finding a bad seed would be *extremely* difficult, but a few would likely exist
	for (int x = 0; x < ITERATION_SIZE; x++) if (!(x & 3)) iteration_table[x] = source_rng->raw32(); else iteration_table[x] = iteration_table[x & ~3];
	for (int x = 0; x < INDIRECTION_SIZE; x++) indirection_table[x] = source_rng->raw32();

	//ought to be a secure seeding if source_rng->get_flags() includes FLAG::CRYPTOGRAPHIC_SECURITY, so...
	for (int x = 0; x < ITERATION_SIZE; x++) raw64();
}*/
void PractRand::RNGs::Raw::efiix64x48::walk_state(StateWalkingObject *walker) {
	for (unsigned long w = 0; w < ITERATION_SIZE; w++) walker->handle(iteration_table[w]);
	for (unsigned long w = 0; w < INDIRECTION_SIZE; w++) walker->handle(indirection_table[w]);
	walker->handle(i);
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
}

