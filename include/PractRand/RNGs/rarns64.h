#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/rarns.cpp
			class rarns64 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 64;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED;
			protected:
				uint64_t xs1, xs2, xs3;
			public:
				uint64_t raw64();
				void seed(uint64_t s);
				void seed(uint64_t s1, uint64_t s2);
				void walk_state(StateWalkingObject* walker);
			};
		}

		namespace Polymorphic {
			class rarns64 final : public vRNG64 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(rarns64)
				void seed(uint64_t s) override;
				void seed(uint64_t s1, uint16_t s2);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(rarns64)
}
