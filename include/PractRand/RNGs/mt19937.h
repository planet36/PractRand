#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			class mt19937 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::OUTPUT_IS_BUFFERED | FLAG::OUTPUT_IS_HASHED | FLAG::ENDIAN_SAFE;
			protected:
				static constexpr int ARRAY_SIZE = 624;
				static constexpr int OFFSET = 397;
				uint32_t state[ARRAY_SIZE];
				uint32_t used;
				void _advance_state();
			public:
				void flush_buffers() {used = ARRAY_SIZE;}
				uint32_t raw32();
				void walk_state(StateWalkingObject* walker);

				//seeds < 2**32 use the standard MT19937 seeding algorithm
				//seeds >= 2**32 use a nonstandard MT19937 seeding algorithm
				void seed(uint64_t s);
				void seed(const uint32_t s[], int seed_length);//alternate seeding algorithm added to MT in 2002

				uint32_t untempered_raw32() {
					if ( used >= ARRAY_SIZE ) {
						_advance_state();
						return state[used++];
					}
					else { return state[used++]; }
				}
				static void self_test();
			};
		}

		namespace Polymorphic {
			class mt19937 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(mt19937)
				mt19937 (uint32_t s[], int seed_length) {seed(s, seed_length);}
				void seed(uint64_t s) override;
				void seed(uint32_t s[], int seed_length);//alternate seeding algorithm added to MT in 2002
				void flush_buffers() override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(mt19937)
}
