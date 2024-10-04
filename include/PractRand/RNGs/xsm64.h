#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/xsm.cpp
			class xsm64 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 64;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION | FLAG::SUPPORTS_FASTFORWARD | FLAG::OUTPUT_IS_HASHED;
			protected:
				Uint64 lcg_low, lcg_high, lcg_adder_low, lcg_adder_high;//2**127 cycles of length 2**128
				void step_backwards();
				void step_forwards();
			public:
				Uint64 raw64();
				void seed(Uint64 s) { seed(s, 0); } // no two seeds on the same cycle
				void seed(Uint64 seed_low, Uint64 seed_high); // no two seeds within 2**127 of each other on the same cycle
				void seed(vRNG *seeder_rng); // no two distinct seeded states within 2**95 of each other on the same cycle (2**160 distinct seeded states possible)
				void walk_state(StateWalkingObject *walker);
				void seek_forward(Uint64 how_far_low64, Uint64 how_far_high64);
				void seek_backward(Uint64 how_far_low64, Uint64 how_far_high64);
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class xsm64 final : public vRNG64 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(xsm64)
				void seed(Uint64 s) override; // no two seeds on the same cycle
				void seed(Uint64 seed_low, Uint64 seed_high); // no two seeds within 2**127 of each other on the same cycle
				void seed(vRNG *seeder_rng) override; // no two distinct seeded states within 2**95 of each other on the same cycle (2**160 distinct seeded states possible)
				void seek_forward128(Uint64 how_far_low64, Uint64 how_far_high64) override;
				void seek_backward128(Uint64 how_far_low64, Uint64 how_far_high64) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(xsm64)
}
