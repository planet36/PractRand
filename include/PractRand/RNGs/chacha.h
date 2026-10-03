#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			//implemented in RNGs/chacha.cpp
			class chacha {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::OUTPUT_IS_BUFFERED | FLAG::OUTPUT_IS_HASHED | FLAG::ENDIAN_SAFE | FLAG::CRYPTOGRAPHIC_SECURITY;
			protected:
				uint32_t outbuf[16]{};
				uint32_t state[12]{};
				uint32_t used{};
				uint32_t position_overflow{};
				uint8_t rounds{20};
				bool extend_cycle{};//true allows carries from the position field to overflow into the upper word of the IV
				bool short_seed{};
				uint8_t padding[5]{};//just to make the size a round number

				void _advance_1();
				//void _reverse_1();
				void _set_position(uint64_t low, uint64_t high);
				void _get_position(uint64_t& low, uint64_t& high) const;

				void _core();
				uint32_t _refill_and_raw32();
			public:
				chacha() = default;
				~chacha();
				uint32_t raw32() {
					if (used < 16) return outbuf[used++];
					return _refill_and_raw32();
				}
				void seed(uint64_t s);
				void seed( const uint32_t seed_and_iv[10], bool extend_cycle_ = false );
				void seed_short( const uint32_t seed_and_iv[6], bool extend_cycle_ = false );
				void walk_state(StateWalkingObject* walker);
				void seek_forward (uint64_t how_far_low, uint64_t how_far_high);
				void seek_backward(uint64_t how_far_low, uint64_t how_far_high);

				//normally rounds is 8, 12, or 20, but it can be anywhere from 1 to 255
				//the default is 20
				//the author recommends 8 for weak crypto or non-crypto, 12 for moderate crypto, and 20 for strong crypto
				//I'd say that you can go as low as 4 or 5 rounds without adversely effecting the output for non-cryptographic purposes
				//when using an odd number of rounds the final transposition gets skipped - this may or may not match other implementations
				//therefore it is recommended that only even numbers of rounds be used
				void set_rounds(int rounds_);
				[[nodiscard]] int get_rounds() const {return rounds;}

				static void self_test();
			};
		}
		namespace Polymorphic {
			class chacha final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(chacha)
				explicit chacha(uint32_t seed_and_iv[10], bool extend_cycle_ = false) {seed(seed_and_iv, extend_cycle_);}
				void seed(uint64_t s) override;
				void seed(uint32_t seed_and_iv[10], bool extend_cycle_ = false);
				void seed_short(uint32_t seed_and_iv[6], bool extend_cycle_ = false);
				void seek_forward128 (uint64_t how_far_low64, uint64_t how_far_high64) override;
				void seek_backward128(uint64_t how_far_low64, uint64_t how_far_high64) override;

				//normally rounds is 8, 12, or 20, but lower and higher values are also possible
				//default is 20
				//for non-crypto applications, 4 rounds is sufficient to qualify for a 3 star quality rating, 6 rounds for a 5 star quality rating
				//for crypto applications, 8 rounds is sufficient to qualify for a 1 star crypto rating, 12 for a 3 star crypto rating, 20 for a 4 star crypto rating
				void set_rounds(int rounds_);
				[[nodiscard]] int get_rounds() const;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(chacha)
}
