#pragma once

#include "PractRand/test_helpers.h"
#include "PractRand/tests.h"

namespace PractRand::Tests {
		class CoupGap final : public TestBaseclass {
			//static constexpr int MAX_OLDEST_AGE = 256 * 16;
			static constexpr int MAX_CURRENT_AGE = 256 * 2;

			//low8 = current_sym, high8 = oldest_sym
			FixedSizeCount<uint8_t, 256*256> count_syms_by_oldest_sym;

			//low8 = oldest_sym, high8 = current_gap
//			FixedSizeCount<uint8_t, 256*MAX_CURRENT_AGE> count_gaps_by_oldest_sym;

//			uint32_t last_sym_pos[256];
			unsigned long oldest_sym{};
			unsigned long youngest_sym{};
			uint8_t next_younger_sym[256]{};
			bool sym_has_appeared[256]{};
			int symbols_ready{};

			uint32_t autofail{};
//			uint64_t blocks;
		public:
			CoupGap() = default;

			void init( RNGs::vRNG* known_good ) override;
			[[nodiscard]] std::string get_name() const override;
			void get_results ( std::vector<TestResult>& results ) override;
			void test_blocks(TestBlock* data, int numblocks) override;
		};
}//PractRand
