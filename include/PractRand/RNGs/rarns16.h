#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/rarns.cpp
			class rarns16 {
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 16,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED
				};
			protected:
			public:
				Uint16 xs1, xs2, xs3;
			public:
				Uint16 raw16();
				void seed(Uint64 s);
				void seed(Uint16 s1, Uint16 s2, Uint16 s3);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class rarns16 final : public vRNG16 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(rarns16)
				void seed(Uint64 s) override;
				void seed(Uint16 s1, Uint16 s2, Uint16 s3);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(rarns16)
}
