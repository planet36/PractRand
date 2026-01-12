namespace PractRand {
	namespace Tests {
		class Coup16 : public TestBaseclass {
		protected:
			enum { S = 65536 / 32 };
			Uint32 flags[S];
			FixedSizeCount<Uint8, 65536> counts;
			//to do: also measure sequential correlation, or do something like BCFN on whether or not each sample is >= the expected value
		public:
			Coup16() {}
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void test_blocks(TestBlock *data, int numblocks);
		};

		class Coup32 : public TestBaseclass {
		protected:
			typedef Uint64 BitmapWord;//must be either 64 bit or 32 bit
			enum {
				BITMAP_WORD_BITS_L2 = sizeof(BitmapWord) == sizeof(Uint64) ? 6 : (sizeof(BitmapWord) == sizeof(Uint32) ? 5 : -1),
				COUPON_BITS = 32, // hardwired to 32, don't change
			};
			BitmapWord bitmap[1ull << (COUPON_BITS - BITMAP_WORD_BITS_L2)];
			Uint64 set_start_position;
			Uint64 filled_words;
			Uint64 num_complete_sets;
			double sum_log2_set_size;
			double sum_sqr_log2_set_size;
		public:
			Coup32() {}
			void handle_set_completion(Uint64 position);
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void test_blocks(TestBlock *data, int numblocks);
		};

	}//Tests
}//PractRand
