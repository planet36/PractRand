#pragma once

namespace PractRand::Tests {
		/*
			to do:
				consider more complex forms, like a BCFN-equivalent or somesuch
		*/
		class mod3_simple : public TestBaseclass {
		public:
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;
		protected:
			typedef Uint32 Word;
			static constexpr int WORD_BITS = 8 * sizeof(Word);
			static constexpr int EXP = 11;//the number of words used in an overlapping sample; 9 matches what gjrand does, but I think 10 is good given cache sizes
			static constexpr int K = ((EXP & 1) ? 3 : 1) * ((EXP & 2) ? 9 : 1) * ((EXP & 4) ? 81 : 1) * ((EXP & 8) ? 6561 : 1);
			static constexpr int KRNDUPA = K | (K >> 1) | (K >> 2) | (K >> 3) | (K >> 4) | (K >> 5) | (K >> 6) | (K >> 7);
			static constexpr int KRNDUPB = (KRNDUPA | (KRNDUPA >> 8) | (KRNDUPA >> 16) | (KRNDUPA >> 24)) + 1;
			static constexpr int BITS = int(1.5849625007211561814537389439478 * EXP + 1);
			static constexpr int P2 = 1 << BITS;
			FixedSizeCount<Uint8, P2> counts;
			unsigned long index;
			void update_index(Word byte);
		};
		class mod3n : public TestBaseclass {
		public:
			mod3n(int block_fraction_);//0 is all, 1 is half, 2 is a quarter, 3 is an 8th, etc
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;
		protected:
			int block_fraction;
			int block_scale;
			Sint64 block_phase;
			Sint64 total_blocks_on;
			static constexpr int EXP = 9;//the number of words used in an overlapping sample; 9 matches what gjrand does, I'd like 10 but in the leveled version with all the extra cache used it may not be worth it
			static constexpr int K = ((EXP & 1) ? 3 : 1) * ((EXP & 2) ? 9 : 1) * ((EXP & 4) ? 81 : 1) * ((EXP & 8) ? 6561 : 1);
			static constexpr int KRNDUPA = K | (K >> 1) | (K >> 2) | (K >> 3) | (K >> 4) | (K >> 5) | (K >> 6) | (K >> 7);
			static constexpr int KRNDUPB = (KRNDUPA | (KRNDUPA >> 8) | (KRNDUPA >> 16) | (KRNDUPA >> 24)) + 1;
			static constexpr int BITS = int(1.5849625007211561814537389439478 * EXP + 1);
			static constexpr int P2 = 1 << BITS;
			static constexpr int LEVELS = 16;
			static constexpr bool PACKED_INDEX = true;
			struct PerLevel {
				unsigned long index;
				Uint8 remainder;
				Uint8 warmup;
				bool odd;
				FixedSizeCount<Uint16, PACKED_INDEX ? K : P2> counts;
			};
			PerLevel levels[LEVELS];
			static unsigned long update_index(unsigned long index, Uint8 remainder);
			void handle_level(int level, Uint8 remainder);
		};
}//PractRand
