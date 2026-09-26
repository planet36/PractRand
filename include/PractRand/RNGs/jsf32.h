#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/jsf.cpp
			class jsf32 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::ENDIAN_SAFE;
			protected:
				uint32_t a, b, c, d;
			public:
				uint32_t raw32();
				void seed(uint64_t s);
				void seed_fast(uint64_t s);
				void seed(vRNG* seeder_rng);
				void seed(uint32_t seed1, uint32_t seed2, uint32_t seed3, uint32_t seed4);//custom seeding
				void walk_state(StateWalkingObject* walker);
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class jsf32 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(jsf32)
				void seed(uint64_t s) override;
				void seed_fast(uint64_t s) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(jsf32)
}
