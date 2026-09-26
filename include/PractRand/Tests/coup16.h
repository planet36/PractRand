#pragma once

#include "PractRand/test_helpers.h"
#include "PractRand/tests.h"

namespace PractRand::Tests {
		class Coup16 final : public TestBaseclass {
		protected:
			static constexpr int S = 65536 / 32;
			uint32_t flags[S]{};
			FixedSizeCount<uint8_t, 65536> counts;
			//to do: also measure sequential correlation, or do something like BCFN on whether or not each sample is >= the expected value
		public:
			Coup16() = default;
			void init(PractRand::RNGs::vRNG* known_good) override;
			[[nodiscard]] std::string get_name() const override;
			void get_results(std::vector<TestResult>& results) override;
			void test_blocks(TestBlock* data, int numblocks) override;
		};
}//PractRand
