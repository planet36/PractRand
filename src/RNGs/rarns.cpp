#include "PractRand/RNGs/rarns16.h"
#include "PractRand/RNGs/rarns32.h"
#include "PractRand/RNGs/rarns64.h"
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"
#include <bit>
#include <string>

using namespace PractRand;
using namespace PractRand::Internals;

//polymorphic:
PRACTRAND_POLYMORPHIC_RNG_BASICS_C64(rarns64)
void PractRand::RNGs::Polymorphic::rarns64::seed(Uint64 s) { implementation.seed(s); }
//void PractRand::RNGs::Polymorphic::rarns64::seed_fast(Uint64 s) { implementation.seed_fast(s); }
//void PractRand::RNGs::Polymorphic::rarns64::seed(Uint64 s1, Uint64 s2, Uint64 s3) { implementation.seed(s1, s2, s3); }
std::string PractRand::RNGs::Polymorphic::rarns64::get_name() const { return "rarns64"; }

PRACTRAND_POLYMORPHIC_RNG_BASICS_C32(rarns32)
void PractRand::RNGs::Polymorphic::rarns32::seed(Uint64 s) { implementation.seed(s); }
//void PractRand::RNGs::Polymorphic::rarns32::seed_fast(Uint64 s) { implementation.seed_fast(s); }
//void PractRand::RNGs::Polymorphic::rarns32::seed(Uint32 s1, Uint32 s2, Uint32 s3) { implementation.seed(s1, s2, s3); }
std::string PractRand::RNGs::Polymorphic::rarns32::get_name() const { return "rarns32"; }

PRACTRAND_POLYMORPHIC_RNG_BASICS_C16(rarns16)
void PractRand::RNGs::Polymorphic::rarns16::seed(Uint64 s) { implementation.seed(s); }
//void PractRand::RNGs::Polymorphic::rarns16::seed_fast(Uint64 s) { implementation.seed_fast(s); }
//void PractRand::RNGs::Polymorphic::rarns16::seed(Uint16 s1, Uint16 s2, Uint16 s3) { implementation.seed(s1, s2, s3); }
std::string PractRand::RNGs::Polymorphic::rarns16::get_name() const { return "rarns16"; }

//raw:
static char global_s1 = 3, global_s2 = 7, global_s3 = 8;
void set_shift_values(int shift1, int shift2, int shift3) {
	global_s1 = shift1;
	global_s2 = shift2;
	global_s3 = shift3;
}
Uint16 PractRand::RNGs::Raw::rarns16::raw16() {
	//constexpr int S1 = 2;
	//constexpr int S2 = 1;
	//constexpr int S3 = 1;
	constexpr int OUTROT = 5;
	Uint16 rv = std::rotl(static_cast<uint16_t>(xs1 + xs3), OUTROT);
	Uint16 old = xs1 >> global_s1;
	xs3 ^= xs1;
	xs1 ^= xs2;
	xs2 ^= xs3;
	xs1 ^= old;
	xs3 = std::rotl(xs3, global_s2);
	xs1 = std::rotl(xs1, global_s3);
	rv += xs1;
	return rv;
}
void PractRand::RNGs::Raw::rarns16::seed(Uint64 s) {
	xs1 = Uint16(s);
	xs2 = Uint16(s >> 16);
	xs3 = Uint16(s >> 32);
	if (!(xs1 | xs2 | xs3)) {
		xs1 = Uint16(s >> 48);
		xs3 = ~xs3;
		//chosen so that 0 through 2**48-2 (inclusive) are all distinct seeds
	}
}
void PractRand::RNGs::Raw::rarns16::seed(Uint16 s1, Uint16 s2, Uint16 s3) {
	xs1 = s2; xs2 = s2; xs3 = s3;
	if (!(s1 | s2 | s3)) xs1 = 1;
}
void PractRand::RNGs::Raw::rarns16::walk_state(StateWalkingObject *walker) {
	walker->handle(xs1);
	walker->handle(xs2);
	walker->handle(xs3);
	if (!(xs1 | xs2 | xs3)) xs1 = 1;
}

Uint32 PractRand::RNGs::Raw::rarns32::raw32() {
	constexpr int S1 = 13;
	constexpr int S2 = 3;
	constexpr int S3 = 5;
	constexpr int OUTROT = 10;
	Uint32 rv = std::rotl(xs1 + xs3, OUTROT);
	Uint32 old = xs1 >> S1;
	xs3 ^= xs1;
	xs1 ^= xs2;
	xs2 ^= xs3;
	xs1 ^= old;
	xs3 = std::rotl(xs3, S2);
	xs1 = std::rotl(xs1, S3);
	rv += xs1;
	return rv;
}
void PractRand::RNGs::Raw::rarns32::seed(Uint64 s) {
	xs1 = Uint32(s >> 0);
	xs2 = Uint32(s >> 32);
	xs3 = ~(xs1 + xs2);

	for (int i = 0; i < 5; i++) raw32();
}
void PractRand::RNGs::Raw::rarns32::seed(Uint32 s1, Uint32 s2, Uint32 s3) {
	xs1 = s2; xs2 = s2; xs3 = s3;
	if (!(s1 | s2 | s3)) xs1 = 1;
}
void PractRand::RNGs::Raw::rarns32::walk_state(StateWalkingObject *walker) {
	walker->handle(xs1);
	walker->handle(xs2);
	walker->handle(xs3);
	if (!(xs1 | xs2 | xs3)) xs1 = 1;
}
Uint64 PractRand::RNGs::Raw::rarns64::raw64() {
	constexpr int S1 = 13;
	constexpr int S2 = 3;
	constexpr int S3 = 5;
	constexpr int OUTROT = 21;
	Uint64 rv = std::rotl(xs1 + xs3, OUTROT);
	Uint64 old = xs1 >> S1;
	xs3 ^= xs1;
	xs1 ^= xs2;
	xs2 ^= xs3;
	xs1 ^= old;
	xs3 = std::rotl(xs3, S2);
	xs1 = std::rotl(xs1, S3);
	rv += xs1;
	return rv;
}
void PractRand::RNGs::Raw::rarns64::seed(Uint64 s) {
	xs1 = s;
	xs2 = xs3 = 1;
	for (int i = 0; i < 8; i++) raw64();
}
void PractRand::RNGs::Raw::rarns64::seed(Uint64 s1, Uint64 s2) {
	xs1 = s1;
	xs2 = s2;
	xs3 = 1;
	for (int i = 0; i < 8; i++) raw64();
}
void PractRand::RNGs::Raw::rarns64::walk_state(StateWalkingObject *walker) {
	walker->handle(xs1);
	walker->handle(xs2);
	walker->handle(xs3);
	if (!(xs1 | xs2 | xs3)) xs1 = 1;
}


