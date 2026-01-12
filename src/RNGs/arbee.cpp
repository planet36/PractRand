#include <string>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>
#include "PractRand/endian.h"

#include "PractRand/RNGs/arbee.h"

using namespace PractRand;
using namespace PractRand::Internals;

//raw:
void PractRand::RNGs::Raw::arbee::reset_entropy() {
	/*
	a = 1;
	b = 2;
	c = 3;
	d = 4;
	i = 0;
	for (int x = 0; x < 8; x++) raw64();
	//*/

	//skipping the math to just write the conclusion:

	a = 6297939664007638733ull;
	b = 10847301951251065963ull;
	c = 8760839667624891872ull;
	d = 11710998636947437130ull;
	i = 8;
}
Uint64 PractRand::RNGs::Raw::arbee::raw64() {
	// very loosely based upon the 3-rotate variation of the jsf64 algorithm

	// the original:
	/*
	Uint64 e = a - rotate(b,7);
	a = b ^ rotate(c, 13);
	b = c + rotate(d, 37);
	c = e + d;
	d = e + a;
	return d;//*/

	// modified slightly
	//  1. a counter was added to impose a minimum cycle length.  It probably won't ever matter in practice, but people like that, and it didn't hurt too much, and I wanted to do something to make sure zero states were fine anyway.  
	//  2. I changed the subtraction to an addition instead.  That's better for LEA on x86, an in my experience subtraction isn't likely to help this sort of thing.  Really this change was mostly pointless though.  
	//  3. shift constants were adjusted to maximize performance in -ttseed64 mode, as that most closely approximates what we want for add_entropy calls.  In order to make this feasible, I reduced the number of output skipped in seeding down to 5 for a first pass, then 6 once I had a number of candidates.  I didn't test all possibilities exhaustively though, just maybe a few hundred of the more likely-seeming triples.  
	//			the resulting best shift triple was {21,39,28}, there was a tie for 2nd place between {13,5,34} and {43,24,29}, and a rough tie for 4th between {45,11,37} and {39,16,35}
	//*
	Uint64 e = a + rotate(b, 21);//21
	a = b ^ rotate(c, 39);//39
	b = c + rotate(d, 28);//28
	c = e + d + i++;
	d = e + a;
	return d;//*/
}
void PractRand::RNGs::Raw::arbee::seed(Uint64 seed1, Uint64 seed2) {
	a = c = seed1;
	b = d = seed2;
	i = 0;
	for (int x = 0; x < 8; x++) raw64();//8
	/*
		(# of outputs discarded) vs (log2 of # of seeds needed to detect interseed correlation)
			2: 9
			3: 13
			4: 18
			5: 26
			6: >38
			outdated, redo
	*/
}
void PractRand::RNGs::Raw::arbee::seed(vRNG *rng) {
	Uint64 seed1 = rng->raw64();
	Uint64 seed2 = rng->raw64();
	seed(seed1, seed2);
}
void PractRand::RNGs::Raw::arbee::autoseed() {
	StateWalkingObject *walker = Internals::get_autoseeder(this);
	this->walk_state(walker);
	delete walker;
}
void PractRand::RNGs::Raw::arbee::walk_state(StateWalkingObject *walker) {
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
	walker->handle(d);
	walker->handle(i);
	if (walker->is_seeder()) i = 0;
}
void PractRand::RNGs::Raw::arbee::add_entropy_N(const void *_data, size_t length) {
	const Uint8 *data = (const Uint8*) _data;
	while (length >= 8) {
		Uint64 in = Uint64(*(data++));
		in |= Uint64(*(data++)) << 8;
		in |= Uint64(*(data++)) << 16;
		in |= Uint64(*(data++)) << 24;
		in |= Uint64(*(data++)) << 32;
		in |= Uint64(*(data++)) << 40;
		in |= Uint64(*(data++)) << 48;
		in |= Uint64(*(data++)) << 56;
		add_entropy64(in);
		length -= 8;
	}
	if (length >= 4) {
		Uint32 in = Uint32(*(data++));
		in |= Uint32(*(data++)) << 8;
		in |= Uint32(*(data++)) << 16;
		in |= Uint32(*(data++)) << 24;
		add_entropy32(in);
		length -= 4;
	}
	if (length >= 2) {
		Uint16 in = Uint16(*(data++));
		in |= Uint16(*(data++)) << 8;
		add_entropy16(in);
		length -= 2;
	}
	if (length) add_entropy8(*data);
}
void PractRand::RNGs::Raw::arbee::add_entropy8(Uint8 value) {
	add_entropy16(value);
}
void PractRand::RNGs::Raw::arbee::add_entropy16(Uint16 value) {
	b ^= value;
	raw64();
}
void PractRand::RNGs::Raw::arbee::add_entropy32(Uint32 value) {
	b ^= value;
	d += value;
	raw64();
}
void PractRand::RNGs::Raw::arbee::add_entropy64(Uint64 value) {
	b ^= value;
	d += value;
	raw64();
	raw64();
}
void PractRand::RNGs::Raw::arbee::flush_buffers() {
	for (int x = 0; x < 7; x++) raw64();//7
}


//polymorphic:
Uint64 PractRand::RNGs::Polymorphic::arbee::get_flags() const {
	return FLAGS;
}
std::string PractRand::RNGs::Polymorphic::arbee::get_name() const {
	return std::string("arbee");
}
Uint8  PractRand::RNGs::Polymorphic::arbee::raw8 () {
	return Uint8(implementation.raw64());
}
Uint16 PractRand::RNGs::Polymorphic::arbee::raw16() {
	return Uint16(implementation.raw64());
}
Uint32 PractRand::RNGs::Polymorphic::arbee::raw32() {
	return Uint32(implementation.raw64());
}
Uint64 PractRand::RNGs::Polymorphic::arbee::raw64() {
	return Uint64(implementation.raw64());
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy_N(const void *data, size_t length) {
	implementation.add_entropy_N(data, length);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy8 (Uint8  value) {
	implementation.add_entropy8 (value);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy16(Uint16 value) {
	implementation.add_entropy16(value);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy32(Uint32 value) {
	implementation.add_entropy32(value);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy64(Uint64 value) {
	implementation.add_entropy64(value);
}
void PractRand::RNGs::Polymorphic::arbee::seed(Uint64 seed_low, Uint64 seed_high) {
	implementation.seed(seed_low, seed_high);
}
void PractRand::RNGs::Polymorphic::arbee::seed(vRNG *rng) {
	implementation.seed(rng);
}
void PractRand::RNGs::Polymorphic::arbee::reset_entropy() {
	implementation.reset_entropy();
}
void PractRand::RNGs::Polymorphic::arbee::walk_state(StateWalkingObject *walker) {
	implementation.walk_state(walker);
}
void PractRand::RNGs::Polymorphic::arbee::flush_buffers() {
	implementation.flush_buffers();
}
