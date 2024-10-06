#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			class isaac32x256 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::CRYPTOGRAPHIC_SECURITY | FLAG::OUTPUT_IS_BUFFERED | FLAG::ENDIAN_SAFE;
			protected:
				static constexpr int SIZE_L2 = 8;
				static constexpr int SIZE = 1 << SIZE_L2;
				Uint32 results[SIZE];
				Uint32 used;
				Uint32 state[SIZE];
				Uint32 a, b, c;
				void _advance_state();
				void _seed(bool flag = true);
			public:
				~isaac32x256();
				void flush_buffers() {used = SIZE;}
				Uint32 raw32() {//LOCKED, do not change
					//note: this walks the buffer in the same direction as the buffer is filled
					//  whereas (some of) Bob Jenkins original code walked the buffer backwards
					if ( used >= SIZE ) _advance_state();
					return results[used++];
				}
				void seed(Uint64 s);
				void seed(Uint32 s[256]);
				void seed(vRNG *seeder_rng);
				void walk_state(StateWalkingObject *walker);
				static void self_test();
			};
		}

		namespace Polymorphic {
			class isaac32x256 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(isaac32x256)
				void seed(Uint64 s) override;
				void seed(vRNG *seeder_rng) override;
				void flush_buffers() override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(isaac32x256)
}
