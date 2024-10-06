#pragma once

#include "PractRand/test_helpers.h"
#include "PractRand/tests.h"

namespace PractRand::Tests {
		/*
			BCFN now has a 3rd parameter.  Set to false it acts like the BCFN in versions 0.90 and earlier.
				Set to false it acts in a broadly similar manner, but is faster.  I don't *think* it's any worse, but am not sure.
			Also, BCFN now has some variants:
				BCFN_FF - slightly slower, but perhaps slightly more powerful?
					looks like it doesn't add much
				//BCFN_MT - MUCH slower, but hopefully more powerful?
					not really implemented yet
		*/
		class BCFN final : public TestBaseclass {
		protected:
			static constexpr int LEVELS = 32;
			VariableSizeCount<Uint16> counts[LEVELS];
			long mask[LEVELS];
			long cur[LEVELS];
			bool even[LEVELS];     // <4 is never used?
			long leftovers[LEVELS];// <4 is never used? <10 not really needed?
			long warmup[LEVELS];   // <4 is never used?
			long tbits;
			long unitsL2;
			bool unbalanced;
			Uint64 blocks;
		public:
			BCFN(int unitsL2_ = 0, int tbits_ = 10, bool unbalanced_ = true);
			void init( PractRand::RNGs::vRNG *known_good ) override;
			void deinit( ) override;
			std::string get_name() const override;
			//virtual double get_result();
			//virtual double result_to_pvalue ( Uint64 blocks, double r );
			void get_results ( std::vector<TestResult> &results ) override;

			void test_blocks(TestBlock *data, int numblocks) override;
			void handle_high_levels_balanced   ( long level, long bits );
			void handle_high_levels_unbalanced ( long level, long bits );
		};
		class BCFN_FF final : public TestBaseclass {
		protected:
			static constexpr int LEVELS = 32;
			VariableSizeCount<Uint16> counts[LEVELS];
			//long mask[LEVELS];
			unsigned long mask;
			long cur[LEVELS];
			bool even[LEVELS];     // <4 is never used?
			long leftovers[LEVELS];// <4 is never used? <10 not really needed?
			long warmup[LEVELS];   // <4 is never used?

			static constexpr int COUNTS2_SIZE = 256;//review the optimizations near the end of test_blocks() before reducing this value
			Uint64 counts2[LEVELS][COUNTS2_SIZE];
			std::vector<Sint32> extreme_counts2[LEVELS];// <= 1 in a hundred million events
			long shifts[LEVELS];

			long tbits;
			long unitsL2;

			bool unbalanced;
			Uint64 blocks;
		public:
			BCFN_FF(int unitsL2_ = 0, int tbits_ = 10, bool unbalanced_ = true);
			void init( PractRand::RNGs::vRNG *known_good ) override;
			void deinit( ) override;
			std::string get_name() const override;
			//virtual double get_result();
			//virtual double result_to_pvalue ( Uint64 blocks, double r );
			void get_results ( std::vector<TestResult> &results ) override;

			void test_blocks(TestBlock *data, int numblocks) override;
			void handle_high_levels ( int level, int bits );
		};
}//PractRand
