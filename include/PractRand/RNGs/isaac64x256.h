#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			class isaac64x256 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 64;
				static constexpr int FLAGS = FLAG::CRYPTOGRAPHIC_SECURITY | FLAG::OUTPUT_IS_BUFFERED | FLAG::ENDIAN_SAFE;
			protected:
				static constexpr int SIZE_L2 = 8;
				static constexpr int SIZE = 1 << SIZE_L2;
				uint64_t results[SIZE];
				uint64_t state[SIZE];
				uint64_t a, b, c;
				uint32_t used;
				void _advance_state();
				void _seed(bool flag=true);
			public:
				isaac64x256() = default;
				~isaac64x256();
				isaac64x256(const isaac64x256&) = delete;
				isaac64x256& operator=(const isaac64x256&) = delete;
				isaac64x256(isaac64x256&&) = delete;
				isaac64x256& operator=(isaac64x256&&) = delete;
				void flush_buffers() {used = SIZE;}
				uint64_t raw64() {//LOCKED, do not change
					//note: this walks the buffer in the same direction as the buffer is filled
					//  whereas (some of) Bob Jenkins original code walked the buffer backwards
					if ( used >= SIZE ) _advance_state();
					return results[used++];
				}
				void seed(uint64_t s);
				void seed(const uint64_t s[256]);
				void seed(vRNG* seeder_rng);
				void walk_state(StateWalkingObject* walker);
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class isaac64x256 final : public vRNG64 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(isaac64x256)
				void seed(uint64_t s) override;
				void seed(vRNG* seeder_rng) override;
				void flush_buffers() override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(isaac64x256)
}
