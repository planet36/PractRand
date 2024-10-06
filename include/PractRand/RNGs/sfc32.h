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
				Uint32 a, b, c, counter;
			public:
				Uint32 raw32();
				void seed(Uint64 s);
				void seed_fast(Uint64 s);
				void seed(Uint32 s1, Uint32 s2, Uint32 s3);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class sfc32 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(sfc32)
				void seed(Uint64 s) override;
				void seed_fast(Uint64 s) override;
				void seed(Uint32 s1, Uint32 s2, Uint32 s3);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(sfc32)
}
