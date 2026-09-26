#pragma once

#include "PractRand/test_helpers.h"
#include "PractRand/tests.h"

namespace PractRand::Tests {
		//class BitMatrix;
		class BRank final : public TestBaseclass {
		public:
			explicit BRank (
				uint32_t rate_hl2_ // 2 * log2(time units per KB)
			);
			void init( PractRand::RNGs::vRNG* known_good ) override;
			[[nodiscard]] std::string get_name() const override;
			void get_results ( std::vector<TestResult>& results ) override;
			void deinit() override;

			void test_blocks(TestBlock* data, int numblocks) override;
		protected:
			uint32_t rate_hl2;
			uint64_t rate;

			void pick_next_size();
			void finish_matrix();

			uint64_t saved_time{};

			//partially complete matrix:
			BitMatrix* in_progress{nullptr};
			uint32_t blocks_in_progress{};
			int size_index{};//which PerSize is active atm?

			//stats:
			class PerSize {
			public:
				//PerSize() {}
				uint32_t size{};
				uint64_t time_per{};
				uint64_t total{};
				static constexpr int NUM_COUNTS = 10;
				static constexpr int MAX_OUTLIERS = 100;
				uint64_t counts[NUM_COUNTS]{};

				uint64_t outliers_overflow{};
				std::vector<uint32_t> outliers;
				void reset();
			};
			std::vector<PerSize> ps;

		};
}//PractRand
