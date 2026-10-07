#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			class hc256 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::CRYPTOGRAPHIC_SECURITY | FLAG::OUTPUT_IS_HASHED | FLAG::OUTPUT_IS_BUFFERED | FLAG::ENDIAN_SAFE;
			protected:
				static constexpr int OUTPUT_BUFFER_SIZE = 64;//should be a multiple of 16
				uint32_t outbuf[OUTPUT_BUFFER_SIZE];
				uint32_t used;
				uint32_t X[16], Y[16];
				uint32_t P[1024], Q[1024];
				uint16_t counter;
				void _do_batch();
				[[nodiscard, gnu::always_inline]] inline uint32_t _h1(uint32_t x) const;
				[[nodiscard, gnu::always_inline]] inline uint32_t _h2(uint32_t x) const;
				[[gnu::always_inline]] inline void _step_A(uint32_t& u, uint32_t v, uint32_t& a, uint32_t b, uint32_t c, uint32_t d, uint32_t& m) const;
				[[gnu::always_inline]] inline void _step_B(uint32_t& u, uint32_t v, uint32_t& a, uint32_t b, uint32_t c, uint32_t d, uint32_t& m) const;
				[[gnu::always_inline]] inline void _feedback_1(uint32_t& u, uint32_t v, uint32_t b, uint32_t c) const;
				[[gnu::always_inline]] inline void _feedback_2(uint32_t& u, uint32_t v, uint32_t b, uint32_t c) const;
			public:
				hc256() = default;
				~hc256();
				hc256(const hc256&) = delete;
				hc256& operator=(const hc256&) = delete;
				hc256(hc256&&) = delete;
				hc256& operator=(hc256&&) = delete;
				void flush_buffers() {used = OUTPUT_BUFFER_SIZE;}
				uint32_t raw32() {//LOCKED, do not change
					if (used < OUTPUT_BUFFER_SIZE) return outbuf[used++];
					_do_batch();
					return outbuf[used++];
				}
				void walk_state(StateWalkingObject* walker);

				//The standard seeding algorithm for HC-256 uses a sequence
				//  of 16 numbers to seed the state.  The first 8 of those
				//  numbers are called the key and the last 8 are called the
				//  initialization vector.  Each number in the sequence is a
				//  32 bit value.  Seeding is very slow.
				void seed(const uint32_t key_and_iv[16]);
				void seed(uint64_t s);
				void seed(vRNG* seeder_rng);
				static void self_test();
			};
		}

		namespace Polymorphic {
			class hc256 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(hc256)
				void seed(uint64_t s) override;
				void seed(uint32_t key_and_iv[16]);
				void seed(vRNG* seeder_rng) override;
				void flush_buffers() override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(hc256)
}
