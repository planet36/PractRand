#pragma once

#include "PractRand/config.h"

namespace PractRand::Crypto {
		class SHA2_512 {
			using Word = uint64_t;
			Word state[8]{};
			uint64_t length{};
			union InputBlock {
				Word as_word[16];
				uint8_t as_byte[16 * sizeof(Word)];
			};
			InputBlock input_buffer{};
			unsigned long leftover_input_bytes{};
			void process_block();
			void process_final_block();
			static Word endianness_word(Word);
			void endianness_input();
			void endianness_state();
		public:
			static constexpr int RESULT_LENGTH = 64;
			void reset();
			void handle_input(const uint8_t* input, unsigned long length);
			void finish(uint8_t destination[RESULT_LENGTH]);
			SHA2_512() {reset();}
		};
}//PractRand
