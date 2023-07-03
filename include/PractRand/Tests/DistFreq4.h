namespace PractRand {
	namespace Tests {
		class DistFreq4 final : public TestBaseclass {
		public:
			DistFreq4(int blocks_per_) : blocks_per(blocks_per_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;

		protected:
			enum {
				ALIGNMENT1 = 4,//the alignments must be powers of 2
				ALIGNMENT2 = 4,
				// TODO: try it with a 3rd window
				SIZE1 = 4,
				SIZE2 = 4,
				POSITIONS1_L2 = 0,//	4,8:38		0,10:38		0,8(8):>39		0(8),8(8):>39		0,11:39
				POSITIONS2_L2 = 10,
				TOTAL_INDEX_BITS = SIZE1 + SIZE2 + POSITIONS1_L2 + POSITIONS2_L2
			};
			int blocks_till_next;
			int blocks_per;
			FixedSizeCount<Uint8, 1 << TOTAL_INDEX_BITS> counts;
		};
		class TripleFreq final : public TestBaseclass {
		public:
			TripleFreq(int passes_at_once_, int blocks_per_pass_) : blocks_per_pass(blocks_per_pass_), passes_at_once(passes_at_once_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;

		protected:
			static constexpr int BASE_ALIGNMENT = 64;//don't change
			static constexpr int WINDOW_ALIGNMENT = 4;
			static constexpr int SIZE1 = 4;//sizes must be multiples of WINDOW_ALIGNMENT
			static constexpr int SIZE2 = 4;
			static constexpr int SIZE3 = 4;
			static constexpr int POSITIONS2_L2 = 6;
			static constexpr int POSITIONS3_L2 = 6;
			static constexpr int TOTAL_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + POSITIONS2_L2 + POSITIONS3_L2;
			static constexpr int REGION_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + POSITIONS3_L2;
			static constexpr int PASSES_PER_REGION = 1 << 14; //
			static constexpr int NUMBER_OF_REGIONS = 1 << POSITIONS2_L2;
			int regions_tested;
			int passes_till_next_region;
			int blocks_till_next_pass;
			int blocks_per_pass;
			int passes_at_once;
			FixedSizeCount<Uint16, 1 << TOTAL_INDEX_BITS> counts;
			// reordered order:                (pos2 aka region), (pos3), (window1), (window2), (window3)
			// index order, from high to low:  (pos2 aka region), (window1), (window2), (pos3), (window3)
		};
		class TripleMirrorFreq final : public TestBaseclass {
		public:
			TripleMirrorFreq(int passes_at_once_, int blocks_per_pass_) : blocks_per_pass(blocks_per_pass_), passes_at_once(passes_at_once_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;
			virtual int get_blocks_to_repeat() const override;

		protected:
			static constexpr int BASE_ALIGN_L2 = 2;
			static constexpr int POSITION_ALIGN_L2 = 2;
			static constexpr int BLOCK_STEP = 16;
			static constexpr int SIZE1 = 3;//units of bits
			static constexpr int SIZE2 = 3;
			static constexpr int SIZE3 = 3;
			static constexpr int POSITIONS_L2 = 6;//can't exceed (10-SAMPLE_ALIGN_L2) atm
			static constexpr int TOTAL_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + POSITIONS_L2;
			int blocks_till_next_pass;
			int blocks_per_pass;
			int passes_at_once;
			FixedSizeCount<Uint16, 1 << TOTAL_INDEX_BITS> counts;
			// index order, from high to low:  (pos), (window1), (window2), (window3)
		};
		class TripleMirrorFreqN : public TestBaseclass {
		public:
			TripleMirrorFreqN(int minimum_level_) : blocks_per_pass(1 << minimum_level_), minimum_level(minimum_level_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;
			//virtual int get_blocks_to_repeat() const;

		protected:
			static constexpr int MAX_LEVELS = 16;
			static constexpr int ALIGN_L2 = 2;
			static constexpr int ALIGN = 1 << ALIGN_L2;
			static constexpr int SIZE1 = 3;//units of bits, can't exceed ALIGN
			static constexpr int SIZE2 = 3;
			static constexpr int SIZE3 = 3;
			static constexpr int POSITIONS_L2 = 4;//can't exceed (6-ALIGN_L2) atm
			static constexpr int TOTAL_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + POSITIONS_L2;
			Uint64 saved_blocks[MAX_LEVELS * 2];
			char level_state[MAX_LEVELS];
			//states: 
			//0: no blocks saved
			//1: 1 block saved
			//2: 2 blocks saved, in order
			//3: 2 blocks saved, reverse order
			char level_polarity[MAX_LEVELS];// 0: 
			int blocks_till_next_pass;
			int blocks_per_pass;
			int minimum_level;
			FixedSizeCount<Uint16, MAX_LEVELS << TOTAL_INDEX_BITS> counts;
			// index order, from high to low:  (pos), (window1), (window2), (window3)
		};
		class TripleMirrorCoup final : public TestBaseclass {
		public:
			TripleMirrorCoup(int passes_at_once_, int blocks_per_pass_) : blocks_per_pass(blocks_per_pass_), passes_at_once(passes_at_once_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;
			virtual int get_blocks_to_repeat() const override;

		protected:
			static constexpr int BASE_ALIGN_L2 = 2;
			static constexpr int POSITION_ALIGN_L2 = 2;
			static constexpr int BLOCK_STEP = 16;
			static constexpr int SIZE1 = 3;//units of bits
			static constexpr int SIZE2 = 3;
			static constexpr int SIZE3 = 3;
			static constexpr int POSITIONS_L2 = 6;//can't exceed (10-SAMPLE_ALIGN_L2) atm
			static constexpr int TOTAL_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + POSITIONS_L2;
			static constexpr int COUP_BUCKETS = 256;
			Uint64 pass_number;
			int blocks_till_next_pass;
			int blocks_per_pass;
			int passes_at_once;
			FixedSizeCount<Uint16, 1 << TOTAL_INDEX_BITS> counts;
			Uint64 coup_masks[1 << (TOTAL_INDEX_BITS >> 6)];
			Uint64 coup_counts[COUP_BUCKETS];
			Uint64 coup_last[1 << POSITIONS_L2];
			Uint64 coup_collected;
			// index order, from high to low:  (pos), (window1), (window2), (window3)
		};
		class QuadFreq : public TestBaseclass {//not yet implemented
		public:
			QuadFreq(int passes_at_once_, int blocks_per_pass_) : blocks_per_pass(blocks_per_pass_), passes_at_once(passes_at_once_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;

		protected:
			static constexpr int BASE_ALIGNMENT = 64;//don't change
			static constexpr int WINDOW_ALIGNMENT = 3;
			static constexpr int SIZE1 = 3;//sizes must be multiples of WINDOW_ALIGNMENT
			static constexpr int SIZE2 = 3;
			static constexpr int SIZE3 = 3;
			static constexpr int SIZE4 = 3;
			static constexpr int POSITIONS2_L2 = 5;
			static constexpr int POSITIONS3_L2 = 5;
			static constexpr int POSITIONS4_L2 = 5;
			static constexpr int TOTAL_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + SIZE4 + POSITIONS2_L2 + POSITIONS3_L2 + POSITIONS4_L2;
			static constexpr int REGION_INDEX_BITS = SIZE1 + SIZE2 + SIZE3 + POSITIONS4_L2;
			static constexpr int PASSES_PER_REGION = 1 << 14; //
			static constexpr int NUMBER_OF_REGIONS = 1 << POSITIONS2_L2;
			int regions_tested;
			int passes_till_next_region;
			int blocks_till_next_pass;
			int blocks_per_pass;
			int passes_at_once;
			FixedSizeCount<Uint16, 1 << TOTAL_INDEX_BITS> counts;
			// reordered order:                (pos2 aka region), (pos3), (window1), (window2), (window3)
			// index order, from high to low:  (pos2 aka region), (window1), (window2), (pos3), (window3)
		};
		class LPerm16 final : public TestBaseclass {//not yet implemented
		public:
			LPerm16(int word_bits_, int passes_at_once_ = 0, int blocks_per_pass_ = 1) : word_bits(word_bits_), blocks_per_pass(blocks_per_pass_), passes_at_once(passes_at_once_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;

			virtual void test_blocks(TestBlock *data, int numblocks) override;

		protected:
			enum {
				LPERM_BUCKETS = 1 << 15 // do not change
			};
			int word_bits;// 8, 16, 32, or 64 bits; 16 is recommended
			int blocks_till_next_pass;
			int blocks_per_pass;
			int passes_at_once;
			FixedSizeCount<Uint16, LPERM_BUCKETS> lperm_counts;//consider: also doing longer range tests on the comparisons done here
		};
	}//Tests
}//PractRand
