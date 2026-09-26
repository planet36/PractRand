#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/sfc.cpp
			class sfc32 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED;
			protected:
				uint32_t a, b, c, counter;
			public:
				uint32_t raw32();
				void seed(uint64_t s);
				void seed_fast(uint64_t s);
				void seed(uint32_t s1, uint32_t s2, uint32_t s3);
				void walk_state(StateWalkingObject* walker);
			};
		}

		namespace Polymorphic {
			class sfc32 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(sfc32)
				void seed(uint64_t s) override;
				void seed_fast(uint64_t s) override;
				void seed(uint32_t s1, uint32_t s2, uint32_t s3);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(sfc32)
}
