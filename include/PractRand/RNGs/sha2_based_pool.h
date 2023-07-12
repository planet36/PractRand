#pragma once

#include "PractRand/rng_basics.h"

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			class sha2_based_pool final : public vRNG8 {
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL,
					OUTPUT_BITS = 8,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::SUPPORTS_ENTROPY_ACCUMULATION | FLAG::CRYPTOGRAPHIC_SECURITY
				};
				enum {STATE_SIZE = 128-24};
				enum {INPUT_BUFFER_SIZE=128, OUTPUT_BUFFER_SIZE=64};
				Uint8 state[STATE_SIZE];
				Uint8 input_buffer[128];
				Uint8 output_buffer[64];
				Uint16 input_buffer_left, output_buffer_left, state_phase;

				sha2_based_pool(Uint64 s) {seed(s);}
				sha2_based_pool(vRNG *seeder) {seed(seeder);}
				sha2_based_pool(SEED_AUTO_TYPE ) {autoseed();}
				sha2_based_pool(SEED_NONE_TYPE ) {reset_state();}
				sha2_based_pool() {reset_state();}
				~sha2_based_pool();

				std::string get_name() const override;
				Uint64 get_flags() const override;

				Uint8  raw8 () override;
				void seed(Uint64 s) override;
				void reset_state();
				using vRNG::seed;
				void walk_state(StateWalkingObject *walker) override;
				void reset_entropy() override {reset_state();}
				void add_entropy8 (Uint8  value) override;
				void add_entropy16(Uint16 value) override;
				void add_entropy32(Uint32 value) override;
				void add_entropy64(Uint64 value) override;
				void flush_buffers() override;
//				static void self_test();
			protected:
				void empty_input_buffer();
				void refill_output_buffer();
			};
		}
	}
}
