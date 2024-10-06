#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/xsm.cpp
			class xsm32 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION | FLAG::SUPPORTS_FASTFORWARD | FLAG::OUTPUT_IS_HASHED;
			protected:
				Uint32 lcg_low, lcg_high, lcg_adder_low, lcg_adder_high;
				void step_backwards();
				void step_forwards();
			public:
				Uint32 raw32();
				void seed(Uint64 s); // no two seeds within 2**63 of each other on the same cycle
				void seed(vRNG *seeder_rng); // no two seeds within 2**48 of each other on the same cycle (2**79 possible seeded states)
				void walk_state(StateWalkingObject *walker);
				void seek_forward (Uint64 how_far);
				void seek_backward(Uint64 how_far);
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class xsm32 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(xsm32)
				void seed(Uint64 s) override; // no two seeds within 2**63 of each other on the same cycle
				void seed(vRNG *seeder_rng) override; // no two seeds within 2**48 of each other on the same cycle (2**79 possible seeded states)
				void seek_forward128(Uint64 how_far_low64, Uint64 how_far_high64) override;
				void seek_backward128(Uint64 how_far_low64, Uint64 how_far_high64) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(xsm32)
}
