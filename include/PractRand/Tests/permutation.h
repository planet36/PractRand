namespace PractRand {
	namespace Tests {

		class LPerm16 : public TestBaseclass {
		public:
			LPerm16(int word_bits_, int passes_at_once_ = 0, int blocks_per_pass_ = 1) : word_bits(word_bits_), passes_at_once(passes_at_once_), blocks_per_pass(blocks_per_pass_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);

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
		class LPerm16A : public TestBaseclass {
		public:
			LPerm16A(int word_bits_, int passes_at_once_ = 0, int blocks_per_pass_ = 1) : word_bits(word_bits_), passes_at_once(passes_at_once_), blocks_per_pass(blocks_per_pass_) {}
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);

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


		/*class permutation_simple : public TestBaseclass {
		public:
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			typedef Uint32 Word;
			enum { // don't change ANYTHING, this is more internal constants than settings, the code is not flexible
				WORD_BITS = 8 * sizeof(Word),
				WORDS_PER_BLOCK = 5,
				POSSSIBLE_ORDERS = 120,
			};
			FixedSizeCount<Uint8, POSSSIBLE_ORDERS> counts;
			unsigned long index;
			void update_index(Word byte);
		};*/
	}//Tests
}//PractRand
