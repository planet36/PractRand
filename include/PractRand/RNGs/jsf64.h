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
				Uint64 a, b, c, d;
			public:
				Uint64 raw64();
				void seed(Uint64 s);
				void seed_fast(Uint64 s);
				void walk_state(StateWalkingObject *walker);
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class jsf64 final : public vRNG64 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(jsf64)
				void seed(Uint64 s) override;
				void seed_fast(Uint64 s) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(jsf64)
}
