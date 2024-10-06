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
				Uint32 outbuf[OUTPUT_BUFFER_SIZE];
				Uint32 used;
				Uint32 X[16], Y[16];
				Uint32 P[1024], Q[1024];
				Uint16 counter;
				void _do_batch();
			public:
				~hc256();
				void flush_buffers() {used = OUTPUT_BUFFER_SIZE;}
				Uint32 raw32() {//LOCKED, do not change
					if (used < OUTPUT_BUFFER_SIZE) return outbuf[used++];
					_do_batch();
					return outbuf[used++];
				}
				void walk_state(StateWalkingObject *walker);

				//The standard seeding algorithm for HC-256 uses a sequence
				//  of 16 numbers to seed the state.  The first 8 of those
				//  numbers are called the key and the last 8 are called the
				//  initialization vector.  Each number in the sequence is a
				//  32 bit value.  Seeding is very slow.
				void seed(Uint32 key_and_iv[16]);
				void seed(Uint64 s);
				void seed(vRNG *seeder_rng);
				static void self_test();
			};
		}

		namespace Polymorphic {
			class hc256 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(hc256)
				void seed(Uint64 s) override;
				void seed(Uint32 seed_and_iv[16]);
				void seed(vRNG *seeder_rng) override;
				void flush_buffers() override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(hc256)
}
