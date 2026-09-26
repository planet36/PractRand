#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/rarns.cpp
			class rarns32 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED;
			protected:
				uint32_t xs1, xs2, xs3;
			public:
				uint32_t raw32();
				void seed(uint64_t s);
				void seed(uint32_t s1, uint32_t s2, uint32_t s3);// { xs1 = s2; xs2 = s2; xs3 = s3; if (!(s1 | s2 | s3)) xs1 = 1; }
				void walk_state(StateWalkingObject* walker);
			};
		}

		namespace Polymorphic {
			class rarns32 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(rarns32)
				void seed(uint64_t s) override;
				void seed(uint32_t s1, uint32_t s2, uint32_t s3);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(rarns32)
}
