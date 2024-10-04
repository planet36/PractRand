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
				Uint64 xs1, xs2, xs3;
			public:
				Uint64 raw64();
				void seed(Uint64 s);
				void seed(Uint64 s1, Uint64 s2);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class rarns64 final : public vRNG64 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(rarns64)
				void seed(Uint64 s) override;
				void seed(Uint64 s1, Uint16 s2);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(rarns64)
}
