#pragma once

#include "PractRand/rng_basics.h"

namespace PractRand::RNGs::Polymorphic {
			class sha2_based_pool final : public vRNG8 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;
				static constexpr int OUTPUT_BITS = 8;
				static constexpr int FLAGS = FLAG::ENDIAN_SAFE | FLAG::SUPPORTS_ENTROPY_ACCUMULATION | FLAG::CRYPTOGRAPHIC_SECURITY;
				static constexpr int STATE_SIZE = 128 - 24;
				static constexpr int INPUT_BUFFER_SIZE = 128;
				static constexpr int OUTPUT_BUFFER_SIZE = 64;
				uint8_t state[STATE_SIZE]{};
				uint8_t input_buffer[128]{};
				uint8_t output_buffer[64]{};
				uint16_t input_buffer_left{}, output_buffer_left{}, state_phase{};

				explicit sha2_based_pool(uint64_t s) {seed(s);}
				explicit sha2_based_pool(vRNG* seeder) {seed(seeder);}
				explicit sha2_based_pool(SEED_AUTO_TYPE ) {autoseed();}
				explicit sha2_based_pool(SEED_NONE_TYPE ) {reset_state();}
				sha2_based_pool() {reset_state();}
				~sha2_based_pool() override;
				sha2_based_pool(const sha2_based_pool&) = delete;
				sha2_based_pool& operator=(const sha2_based_pool&) = delete;
				sha2_based_pool(sha2_based_pool&&) = delete;
				sha2_based_pool& operator=(sha2_based_pool&&) = delete;

				[[nodiscard]] std::string get_name() const override;
				[[nodiscard]] uint64_t get_flags() const override;

				uint8_t  raw8 () override;
				void seed(uint64_t s) override;
				void reset_state();
				using vRNG::seed;
				void walk_state(StateWalkingObject* walker) override;
				void reset_entropy() override {reset_state();}
				void add_entropy8 (uint8_t  value) override;
				void add_entropy16(uint16_t value) override;
				void add_entropy32(uint32_t value) override;
				void add_entropy64(uint64_t value) override;
				void flush_buffers() override;
//				static void self_test();
			protected:
				void empty_input_buffer();
				void refill_output_buffer();
			};
}
