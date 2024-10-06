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
				Uint16 a, b, c, counter;
			public:
				Uint16 raw16();
				void seed(Uint64 s);
				void seed_fast(Uint64 s);
				void seed(Uint16 s1, Uint16 s2, Uint16 s3);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class sfc16 final : public vRNG16 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(sfc16)
				void seed(Uint64 s) override;
				void seed_fast(Uint64 s) override;
				void seed(Uint16 s1, Uint16 s2, Uint16 s3);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(sfc16)
}
