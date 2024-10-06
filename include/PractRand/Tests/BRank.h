#pragma once

namespace PractRand::Tests {
		//class BitMatrix;
		class BRank final : public TestBaseclass {
		public:
			BRank (
				Uint32 rate_hl2_ // 2 * log2(time units per KB)
			);
			void init( PractRand::RNGs::vRNG *known_good ) override;
			std::string get_name() const override;
			void get_results ( std::vector<TestResult> &results ) override;
			void deinit() override;

			void test_blocks(TestBlock *data, int numblocks) override;
		protected:
			Uint32 rate_hl2;
			Uint64 rate;

			virtual void pick_next_size();
			virtual void finish_matrix();

			Uint64 saved_time;

			//partially complete matrix:
			BitMatrix *in_progress{nullptr};
			Uint32 blocks_in_progress;
			int size_index;//which PerSize is active atm?

			//stats:
			class PerSize {
			public:
				//PerSize() {}
				Uint32 size;
				Uint64 time_per;
				Uint64 total;
				static constexpr int NUM_COUNTS = 10;
				static constexpr int MAX_OUTLIERS = 100;
				Uint64 counts[NUM_COUNTS];

				Uint64 outliers_overflow;
				std::vector<Uint32> outliers;
				void reset();
			};
			std::vector<PerSize> ps;

		};
}//PractRand
