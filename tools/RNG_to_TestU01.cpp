#include <cstdio>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <string>
extern "C" {
#include "TestU01.h"
}
#include "PractRand_full.h"
#include "PractRand/RNGs/all.h"
#include "PractRand/RNGs/other/mt19937.h"
#include "PractRand/RNGs/other/transform.h"
#include "PractRand/RNGs/other/simple.h"
#include "PractRand/RNGs/other/lcgish.h"
#include "PractRand/RNGs/other/oddball.h"
#include "PractRand/RNGs/other/cbuf.h"
#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/other/special.h"
#include "RNG_from_name.h"


//PractRand RNG to TestU01 adaptor
static PractRand::Uint8 reverse_bits8[256];
class PractRand_unif01_buffered_adaptor {
	//the buffering is done for 64 bit RNGs
	//since it doesn't hurt non-64-bit RNGs much I only bothered with the buffering implementation
	int used;
	PractRand::RNGs::PolymorphicRNG *real_rng;
	static void init_reverse_bits() {
		for (unsigned int i = 0; i < 256; i++) {
			PractRand::Uint8 reversed = 0;
			PractRand::Uint8 tmp = i;
			for (int b = 0; b < 8; b++) {
				if (tmp & (1 << b)) reversed |= 1 << (b ^ 7);
			}
			reverse_bits8[i] = reversed;
		}
	}
	static PractRand::Uint32 bit_reverse32(PractRand::Uint32 in) {
		return (PractRand::Uint32(reverse_bits8[(in >> 0) & 255]) << 24) | (PractRand::Uint32(reverse_bits8[(in >> 8) & 255]) << 16) | (PractRand::Uint32(reverse_bits8[(in >> 16) & 255]) << 8) | (PractRand::Uint32(reverse_bits8[(in >> 24) & 255]) << 0);
	}
	static PractRand::Uint32 byte_reverse32(PractRand::Uint32 in) {
		return ((in & 255) << 24) | ((in & 0xff00) << 8) | ((in & 0xff0000) >> 8) | ((in & 0xff000000) >> 24);
	}

	enum {BLOCKS=1};
	enum {N = BLOCKS * PractRand::Tests::TestBlock::SIZE / 4};
	PractRand::Tests::TestBlock blocks[BLOCKS];

	void refill() {
		used = 0;
		blocks[0].fill(real_rng, BLOCKS);
	}
	unsigned long buffered_get32() {
		if (used < N) return blocks[0].as32[used++];
		refill();
		return blocks[0].as32[used++];
	}
	static unsigned long _unif01_buffered_GetBits(void *param, void *rng_) {
		PractRand_unif01_buffered_adaptor *rng = (PractRand_unif01_buffered_adaptor*)rng_;
		return rng->buffered_get32();
	}
	static unsigned long _unif01_buffered_bitreversed_GetBits(void *param, void *rng_) {
		PractRand_unif01_buffered_adaptor *rng = (PractRand_unif01_buffered_adaptor*)rng_;
		return bit_reverse32(rng->buffered_get32());
	}
	static unsigned long _unif01_buffered_bytereversed_GetBits(void *param, void *rng_) {
		PractRand_unif01_buffered_adaptor *rng = (PractRand_unif01_buffered_adaptor*)rng_;
		return byte_reverse32(rng->buffered_get32());
	}
	static double _unif01_buffered_GetU01(void *param, void *rng_) { return _unif01_buffered_GetBits(param, rng_) / 4294967296.0; }
	static double _unif01_buffered_bitreversed_GetU01(void *param, void *rng_) { return _unif01_buffered_bitreversed_GetBits(param, rng_) / 4294967296.0; }
	static double _unif01_buffered_bytereversed_GetU01(void *param, void *rng_) { return _unif01_buffered_bytereversed_GetBits(param, rng_) / 4294967296.0; }
	static void _unif01_Write(void *rng) { printf("(printing RNG states is not supported)\n"); }
public:
	void bind( unif01_Gen *u, unif01_Gen *bit_reversed = NULL, unif01_Gen *byte_reversed = NULL) {
		u->GetBits = _unif01_buffered_GetBits;
		u->GetU01 = _unif01_buffered_GetU01;
		u->Write = _unif01_Write;
		u->name = strdup(real_rng->get_name().c_str());
		u->state = this;
		u->param = NULL;
		if (bit_reversed) {
			init_reverse_bits();
			bit_reversed->GetBits = _unif01_buffered_bitreversed_GetBits;
			bit_reversed->GetU01 = _unif01_buffered_bitreversed_GetU01;
			bit_reversed->Write = _unif01_Write;
			bit_reversed->name = strdup((real_rng->get_name() + " [BIT REVERSED]").c_str());
			bit_reversed->state = this;
			bit_reversed->param = NULL;
		}
		if (byte_reversed) {
			byte_reversed->GetBits = _unif01_buffered_bytereversed_GetBits;
			byte_reversed->GetU01 = _unif01_buffered_bytereversed_GetU01;
			byte_reversed->Write = _unif01_Write;
			byte_reversed->name = strdup((real_rng->get_name() + " [BYTE REVERSED]").c_str());
			byte_reversed->state = this;
			byte_reversed->param = NULL;
		}
	}
	PractRand_unif01_buffered_adaptor(PractRand::RNGs::vRNG *rng) : real_rng(rng), used(N) {}
};
void do_crush(PractRand::RNGs::vRNG *rng, int quality, int n = 1, int reversal_mode = 0) {
	swrite_Basic = FALSE;
	swrite_Parameters = FALSE;
	swrite_Collectors = FALSE;
	swrite_Counters = FALSE;
	swrite_Classes = FALSE;

	unif01_Gen gen, bit_reverse, byte_reverse;
	PractRand_unif01_buffered_adaptor *binding = new PractRand_unif01_buffered_adaptor( rng );
	binding->bind(&gen, &bit_reverse, &byte_reverse);
	if (reversal_mode == 1) {
		gen = bit_reverse;
	}
	else if (reversal_mode == 2) {
		gen = byte_reverse;
	}
	int repeats[512];
	for (int i = 0; i < 512; i++) repeats[i] = n;
	if (quality == 0) bbattery_RepeatSmallCrush(&gen, repeats);
	else if (quality == 1) bbattery_RepeatCrush(&gen, repeats);
	else if (quality == 2) bbattery_RepeatBigCrush(&gen, repeats);
	else if (quality >= 10) bbattery_RepeatRabbit(&gen, 1ull << quality, repeats);
	else throw "invalid quality number";
	delete binding;
}

bool interpret_seed(const std::string &seedstr, PractRand::Uint64 &seed) {
	//would prefer strtol, but that is insufficiently portable when it has to handle 64 bit values
	PractRand::Uint64 value = 0;
	int position = 0;
	if (seedstr.length() >= 3 && seedstr[0] == '0' && seedstr[1] == 'x') position = 2;
	while (position < seedstr.length()) {
		int c = seedstr[position++];
		if (value >> 60) return false;//too long
		value *= 16;
		if (c >= '0' && c <= '9') value += (c-'0');
		else if (c >= 'a' && c <= 'f') value += (c-'a')+10;
		else if (c >= 'A' && c <= 'F') value += (c-'A')+10;
		else return false;//invalid character
	}
	seed = value;
	return true;
}

int main(int argc, char **argv) {
/*
arguments:
	arg0 is the executable name
	arg1 is the RNG to use
	arg2 is 0 for SmallCrush, 1 for Crush, 2 for BigCrush, 10+ for Rabbit
	//arg3 is a bitfield - bits 0 & 1 dictate reversal, 
		00: normal, 01: bit-reversed, 10: byte-reversed, 11: unused
		bit 2 dictates repeated testing (0: 1 time, 1: 4 times)
	arg3 is the 64 bit seed in hexadecimal - if omitted, a random 32 bit value will be chosen instead
*/
	if (argc < 3) {
		std::printf("usage: %s rng_name crush_level [seed] [modifiers]\n", argv[0]);
		std::printf("   rng_name: any recognized PractRand RNG name, e.g. jsf32\n");
		std::printf("   crush_level: 0 or SmallCrush, 1 or Crush, or 2 or BigCrush\n");
		std::printf("   seed: either a 64 bit value in hexadecimal, or the word \"unspecified\" which will cause a 32 bit value to be randomly chosen\n");
		std::printf("   modifiers: not yet implemented\n");
		std::exit(0);
	}
	PractRand::initialize_PractRand();
	RNG_Factories::register_recommended_RNGs();
	RNG_Factories::register_nonrecommended_RNGs();
	RNG_Factories::register_input_RNGs();
	std::string errmsg;
	PractRand::RNGs::vRNG *rng = RNG_Factories::create_rng(argv[1], &errmsg);
	if (!rng) {
		if (errmsg.empty()) std::fprintf(stderr, "unrecognized RNG name \"%s\".  aborting.\n", argv[1]);
		else std::fprintf(stderr, "%s\n", errmsg.c_str());
		std::exit(1);
	};
	long crush_level;
	if (std::strlen(argv[2]) == 1) crush_level = argv[2][0] - '0';
	else if (!std::strcmp(argv[2], "SmallCrush")) crush_level = 0;
	else if (!std::strcmp(argv[2], "Crush")) crush_level = 1;
	else if (!std::strcmp(argv[2], "BigCrush")) crush_level = 2;
	else if (std::strlen(argv[2]) == 8 && !std::strncmp(argv[2], "Rabbit", 6) && argv[2][6] >= '0' && argv[2][6] <= '9' && argv[2][7] >= '0' && argv[2][7] <= '9') crush_level = 10 + (argv[2][6] - '0') * 10 + (argv[2][7] - '0');
	else {
		std::printf("failed to parse crush_level\n");
		std::printf("aborting\n");
		std::exit(1);
	}
	
	PractRand::Uint64 seed;
	bool unspecified_seed = true;
	if (argc == 4 && strcmp(argv[3], "unspecified")) {
		unspecified_seed = !interpret_seed(argv[3], seed);
		if (unspecified_seed) {
			std::printf("failed to parse seed value\n");
			std::printf("aborting\n");
			std::exit(1);
		}
	}
	if (unspecified_seed) {
		PractRand::RNGs::Polymorphic::arbee good_rng(PractRand::SEED_AUTO);
		seed = good_rng.raw32();
	}
	long repetitions = 1;
	long backwards = 0;
	rng->seed(seed);
	std::printf("starting do_crush\nrng = [%s], crush_level = %ld, repetitions = %ld, backwards = %ld, seed = 0x%llX\n", rng->get_name().c_str(), crush_level, repetitions, backwards, seed);
	do_crush(rng, crush_level, repetitions, backwards);
	return 0;
}
