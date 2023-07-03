namespace PractRand {
	namespace Tests {
		class Pat5 final : public TestBaseclass {
		public:
			Pat5();
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;

		protected:
			typedef Uint32 Word;
			static constexpr int WORD_BITS = sizeof(Word) * 8;
			static constexpr int ZERO_FILTER_BITS = 20;//should not exceed WORD_BITS - PATTERN_INDEX_BITS - (1 << PRIMARY_WORD_DISTANCE_BITS) + 1
			static constexpr int PATTERN_INDEX_BITS = 3;
			static constexpr int PRIMARY_WORD_DISTANCE_BITS = 2;
			static constexpr int PRIMARY_WORD_DISTANCE_EXTRA_BITS = 0;
			static constexpr int NUM_SECONDARY_WORDS = 1;
			static constexpr int SECONDARY_WORD_DISTANCE_BITS = 6;
			static constexpr int SECONDARY_WORD_DISTANCE_EXTRA_BITS = 0;
			static constexpr int NUM_TERTIARY_WORDS = 2;
			static constexpr int TERTIARY_WORD_DISTANCE_BITS = 5;
			static constexpr int TERTIARY_WORD_DISTANCE_EXTRA_BITS = 1;
			static constexpr int NUM_QUATERNARY_WORDS = 2;
			static constexpr int QUATERNARY_WORD_DISTANCE_BITS = 3;
			static constexpr int QUATERNARY_WORD_DISTANCE_EXTRA_BITS = 3;
			static constexpr int TABLE_SIZE_L2 = PATTERN_INDEX_BITS + PRIMARY_WORD_DISTANCE_BITS + SECONDARY_WORD_DISTANCE_BITS + TERTIARY_WORD_DISTANCE_BITS + QUATERNARY_WORD_DISTANCE_BITS;
			static constexpr int PATTERN_WIDTH = 1 + 2 * (NUM_SECONDARY_WORDS + NUM_TERTIARY_WORDS + NUM_QUATERNARY_WORDS);
			class Pattern {
			public:
				Word base_pattern[PATTERN_WIDTH];
				Word _padding;
				Sint64 total_count;
			};
			//PractRand::RNGs::Raw::arbee internal_rng;
			//state:
			Pattern patterns[1 << PATTERN_INDEX_BITS];
			FixedSizeCount<Uint8, 1 << TABLE_SIZE_L2> counts;
			//Uint64 lifespan;
			//internal helpers:
			int transform_bitcount_primary(int bit_count) const;
			int transform_bitcount_secondary(int bit_count) const;
			int transform_bitcount_tertiary(int bit_count) const;
		};
	}//Tests
}//PractRand
