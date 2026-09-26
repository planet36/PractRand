#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/sfc.cpp
			class sfc16 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 16;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED;
			protected:
				uint16_t a, b, c, counter;
			public:
				uint16_t raw16();
				void seed(uint64_t s);
				void seed_fast(uint64_t s);
				void seed(uint16_t s1, uint16_t s2, uint16_t s3);
				void walk_state(StateWalkingObject* walker);
			};
		}

		namespace Polymorphic {
			class sfc16 final : public vRNG16 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(sfc16)
				void seed(uint64_t s) override;
				void seed_fast(uint64_t s) override;
				void seed(uint16_t s1, uint16_t s2, uint16_t s3);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(sfc16)
}
