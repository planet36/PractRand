#include "PractRand/RNGs/arbee.h"
#include "PractRand/endian.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <bit>
#include <cstdint>
#include <string>

using namespace PractRand;
using namespace PractRand::Internals;

//raw:
void PractRand::RNGs::Raw::arbee::reset_entropy() {
	//the default state is arrived at by setting all values to 1, then calling mix()
	a = 9873171087373218264ULL;
	b = 10599573592049074392ULL;
	c = 16865209178899817893ULL;
	d = 5013818595375203225ULL;
	i = 13;
}
uint64_t PractRand::RNGs::Raw::arbee::raw64() {
	uint64_t e = a + std::rotl(b,45);
	a = b ^ std::rotl(c, 13);
	b = c + std::rotl(d, 37);
	c = e + d + i++;
	d = e + a;
	return d;
}
void PractRand::RNGs::Raw::arbee::mix() {
	for (int x = 0; x < 12; x++) raw64();
	/*
		(# of outputs discarded) vs (log2 of # of seeds needed to detect interseed correlation)
			4 - 10-11
			5 - 15-16 5
			6 - 22-23 7
			7 - 34-35 12
	*/
}
void PractRand::RNGs::Raw::arbee::seed(uint64_t seed1, uint64_t seed2, uint64_t seed3, uint64_t seed4) {
	a = seed1;
	b = seed2;
	c = seed3;
	d = seed4;
	i = 1;
	mix();
}
void PractRand::RNGs::Raw::arbee::seed(uint64_t s) {
	a = s;
	b = 1;
	c = 2;
	d = 3;
	i = 1;
	mix();
}
void PractRand::RNGs::Raw::arbee::seed(vRNG* rng) {
	a = rng->raw64();
	b = rng->raw64();
	c = rng->raw64();
	d = rng->raw64();
	i = 13;
}
void PractRand::RNGs::Raw::arbee::walk_state(StateWalkingObject* walker) {
	walker->handle(a);
	walker->handle(b);
	walker->handle(c);
	walker->handle(d);
	walker->handle(i);
	if (walker->is_seeder()) i = 1;
}
void PractRand::RNGs::Raw::arbee::add_entropy_N(const void* _data, size_t length) {
	const auto* data = static_cast<const uint8_t*>(_data);
#if 0 // little-endian only
	//maybe add an ifdef to enable misaligned reads at compile time where appropriate?
/*	if (!(7 & reinterpret_cast<uint64_t>(data))) {
	//if (!(((unsigned long)data) & 7)) {
		while (length >= 8) {
			add_entropy64(*(uint64_t*)data);
			data += 8;
			length -= 8;
		}
	}
	else*/
#endif
	while (length >= 8) {
		auto in = uint64_t(*(data++));
		in |= uint64_t(*(data++)) << 8;
		in |= uint64_t(*(data++)) << 16;
		in |= uint64_t(*(data++)) << 24;
		in |= uint64_t(*(data++)) << 32;
		in |= uint64_t(*(data++)) << 40;
		in |= uint64_t(*(data++)) << 48;
		in |= uint64_t(*(data++)) << 56;
		add_entropy64(in);
		length -= 8;
	}
	while (length >= 2) {
		uint16_t in = uint64_t(*(data++));
		in |= uint16_t(*(data++)) << 8;
		add_entropy16(in);
		length -= 2;
	}
	if (length) add_entropy8(*data);
}
void PractRand::RNGs::Raw::arbee::add_entropy8(uint8_t value) {
	add_entropy16(value);
}
void PractRand::RNGs::Raw::arbee::add_entropy16(uint16_t value) {
	d ^= value;
	b += value;
	raw64();
	b += value;
}
void PractRand::RNGs::Raw::arbee::add_entropy32(uint32_t value) {
	d ^= value;
	b += value;
	raw64();
	raw64();
	b += value;
}
void PractRand::RNGs::Raw::arbee::add_entropy64(uint64_t value) {
	d ^= value;
	b += value;
	raw64();
	raw64();
	raw64();
	raw64();
	b += value;
}


//polymorphic:
uint64_t PractRand::RNGs::Polymorphic::arbee::get_flags() const {
	return FLAGS;
}
std::string PractRand::RNGs::Polymorphic::arbee::get_name() const {
	return {"arbee"};
}
uint8_t  PractRand::RNGs::Polymorphic::arbee::raw8 () {
	return uint8_t(implementation.raw64());
}
uint16_t PractRand::RNGs::Polymorphic::arbee::raw16() {
	return uint16_t(implementation.raw64());
}
uint32_t PractRand::RNGs::Polymorphic::arbee::raw32() {
	return uint32_t(implementation.raw64());
}
uint64_t PractRand::RNGs::Polymorphic::arbee::raw64() {
	return implementation.raw64();
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy_N(const void* data, size_t length) {
	implementation.add_entropy_N(data, length);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy8 (uint8_t  value) {
	implementation.add_entropy8 (value);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy16(uint16_t value) {
	implementation.add_entropy16(value);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy32(uint32_t value) {
	implementation.add_entropy32(value);
}
void PractRand::RNGs::Polymorphic::arbee::add_entropy64(uint64_t value) {
	implementation.add_entropy64(value);
}
void PractRand::RNGs::Polymorphic::arbee::seed(uint64_t s) {
	implementation.seed(s);
}
void PractRand::RNGs::Polymorphic::arbee::seed(vRNG* rng) {
	implementation.seed(rng);
}
void PractRand::RNGs::Polymorphic::arbee::seed(uint64_t s1, uint64_t s2, uint64_t s3, uint64_t s4) {
	implementation.seed(s1, s2, s3, s4);
}
void PractRand::RNGs::Polymorphic::arbee::reset_entropy() {
	implementation.reset_entropy();
}
void PractRand::RNGs::Polymorphic::arbee::walk_state(StateWalkingObject* walker) {
	implementation.walk_state(walker);
}
void PractRand::RNGs::Polymorphic::arbee::flush_buffers() {
	implementation.flush_buffers();
}
