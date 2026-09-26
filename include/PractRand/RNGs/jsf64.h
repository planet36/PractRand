#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/jsf.cpp
			class jsf64 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 64;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::ENDIAN_SAFE;
			protected:
				uint64_t a, b, c, d;
			public:
				uint64_t raw64();
				void seed(uint64_t s);
				void seed_fast(uint64_t s);
				void walk_state(StateWalkingObject* walker);
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class jsf64 final : public vRNG64 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(jsf64)
				void seed(uint64_t s) override;
				void seed_fast(uint64_t s) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(jsf64)
}
