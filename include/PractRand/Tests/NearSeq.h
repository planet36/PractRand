namespace PractRand {
	namespace Tests {
		class NearSeq : public TestBaseclass {
		protected:
			typedef Uint64 Word;
			enum {
				/*
					parameterizations of interest:
						blocks*bits=core		thresholds				evaluation
					word			extra bits			derived
					32	14x8=112	128			2/3		87.5, 35 M		meh
					64	9x7=63		128			1/2		1.2 K, 268 M	???
					64	12x5=60		128			1/2		1, 129 K		
					64	10x12=120	128			3/4		13 K, 227 M		
					32?	10x4=40		64			0/1		110, 1 B		
					64	8x8=64		128			2/3		13, 20 K		
					32	8x8=64		64?			1/2		20 K, 13 B		
					32	8x4=32		64			0/1		43, 17 M		
				*/

				WORD_BITS = sizeof(Word)* 8,
				WORD_BITS_L2 = (WORD_BITS == 8) ? 3 : ((WORD_BITS == 16) ? 4 : ((WORD_BITS == 32) ? 5 : ((WORD_BITS == 64) ? 6 : -1))),
				NUM_BUCKETS_L2 = 9,
				NUM_BUCKETS = 1 << NUM_BUCKETS_L2,
				BITS_PER_BLOCK = 7,
				BLOCKS_PER_CORE = NUM_BUCKETS_L2, 
				CORE_SEQUENCE_BITS = BLOCKS_PER_CORE * BITS_PER_BLOCK,

				CORE_WORDS = (CORE_SEQUENCE_BITS + WORD_BITS - 1) / WORD_BITS, 
				EXTRA_WORDS = 2,
				SEQUENCE_WORDS = CORE_WORDS + EXTRA_WORDS,
				SEQUENCE_BITS = SEQUENCE_WORDS * WORD_BITS,
				SEQUENCE_WORD_OFFSET = (EXTRA_WORDS + 1) / 2,
				MAX_ERRORS_PER_BLOCK = 2,
				GOOD_ERRORS_PER_BLOCK = 1,

				MAX_CORE_DISTANCES = CORE_SEQUENCE_BITS / 2,//distances of this much or higher get combined with (one less than this)
				//MAX_CORE_DISTANCES = MAX_ERRORS_PER_BLOCK * BLOCKS_PER_CORE,//distances of this much or higher get combined with (one less than this)
				//RARES_THRESHOLD = 8,//core distances less than or equal to this are required for rares
				//RARES_SIZE_L2 = 16,
			};
			struct Bucket {
				Uint64 sequence[SEQUENCE_WORDS];
			};
			Bucket buckets[NUM_BUCKETS];
			Uint64 core_distances[MAX_CORE_DISTANCES];
			Uint64 sum_extra_distances[MAX_CORE_DISTANCES];//sum of all overall distances at a given core distances... divide by the same index on core_distance to get the average
			//Uint32 rare_index;
			//Uint32 rare_warmup;
			//FixedSizeCount<Uint16, 1 << RARES_SIZE_L2> count_rares;
			//FixedSizeCount<Uint16, NUM_BUCKETS * EXTRA_WORDS * 16> count_region;
			//void handle_rare(Uint8 one);
			Uint8 *lookup_table;//used by core_to_index
			Uint8 *lookup_table2;//used by is_core_good
			int core_to_index(const Word *core) const;//returns -1 on invalid core
			int is_core_good(const Word *core) const;
			int get_core_distance(const Word *core, int bucket_index) const;
			int get_extra_distance(const Word *core, int bucket_index) const;
		public:
			NearSeq();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};



		class NearSeq2 : public TestBaseclass {
		protected:
			typedef Uint64 Word;
			enum {
				/*
				parameterizations of interest:
				blocks*bits=core		thresholds				evaluation
				word			extra bits			derived
				32	14x8=112	128			2/3		87.5, 35 M		meh
				64	9x7=63		128			1/2		1.2 K, 268 M	???
				64	12x5=60		128			1/2		1, 129 K
				64	10x12=120	128			3/4		13 K, 227 M
				32?	10x4=40		64			0/1		110, 1 B
				64	8x8=64		128			2/3		13, 20 K
				32	8x8=64		64?			1/2		20 K, 13 B
				32	8x4=32		64			0/1		43, 17 M

				64	8x16=128	128?		5/6/7	65, 549, 263K
				64	4x32=128	128?		10-13	49, 465, 6784, 159K
				32	8x32=256	64?			12-14	62, 2446, 216K
				64	4x64=256	128?		23-28	47, 217, 1236, 8740, 77K, 868K

				*/

				// stride is one Word
				WORD_BITS = 8*sizeof(Word),
				WORD_BITS_L2 = (WORD_BITS == 8) ? 3 : ((WORD_BITS == 16) ? 4 : ((WORD_BITS == 32) ? 5 : ((WORD_BITS == 64) ? 6 : -1))),
				NUM_BUCKETS_L2 = 4,
				NUM_BUCKETS = 1 << NUM_BUCKETS_L2,
				BITS_PER_BLOCK = 64,

				BLOCKS_PER_CORE = NUM_BUCKETS_L2,
				CORE_SEQUENCE_BITS = BLOCKS_PER_CORE * BITS_PER_BLOCK,

				CORE_WORDS = (CORE_SEQUENCE_BITS + WORD_BITS - 1) / WORD_BITS,
				EXTRA1_FULL_WORDS = 2,//should be an even number, otherwise call it an invalid parameterization
				EXTRA2_FULL_WORDS = 2,//should be an even number, otherwise call it an invalid parameterization
				EXTRA1_PARTIAL_WORD_BITS = CORE_WORDS * WORD_BITS - CORE_SEQUENCE_BITS, //should be less than 1 word, otherwise call it an invalid parameterization
				EXTRA1_BITS = EXTRA1_FULL_WORDS * WORD_BITS + EXTRA1_PARTIAL_WORD_BITS,
				//EXTRA2_BYTES = EXTRA_FULL_WORDS * 
				SEQUENCE_WORDS = CORE_WORDS + EXTRA1_FULL_WORDS,
				SEQUENCE_BITS = SEQUENCE_WORDS * WORD_BITS,
				SEQUENCE_WORD_OFFSET = EXTRA1_FULL_WORDS / 2,
				MAX_HDIST_PER_BLOCK = 27,
				MAX_TOTAL_HDIST = MAX_HDIST_PER_BLOCK * BLOCKS_PER_CORE,
				//MAX_CONCEIVABLE_HDIST = BLOCKS_PER_CORE * (BITS_PER_BLOCK / 2),
				HDIST_BINS = 16,
				HDIST_BIN_SHIFT = 2,

				MAX_LOOKUP_L2 = 12,//otherwise the tables get too big
				CHECK_VALIDITY_EARLY = 0,//if true, checks for a bad core after each block for the first word, then once per word thereafter
			};
			struct Bucket {
				//Word ideal_sequence[SEQUENCE_WORDS];
				Uint64 core_hdist[MAX_TOTAL_HDIST + 1];
				Uint64 extra_counts[HDIST_BINS][EXTRA1_BITS];
				//Uint64 extra_data[HDIST_BINS][EXTRA_BYTES]
				void reset();
			};
			Bucket buckets[NUM_BUCKETS];

			//Uint64 _total_cores;
			//Uint64 _total_valid_cores;
			//Uint64 _total_invalid_cores;

			double prob_of_valid_block;
			double prob_of_valid_core;
			Uint64 base_chances[MAX_HDIST_PER_BLOCK + 1];
			Uint64 total_base_chances;
			double base_probs[MAX_HDIST_PER_BLOCK + 1];
			double normalized_base_probs[MAX_HDIST_PER_BLOCK + 1];
			double full_probs[MAX_TOTAL_HDIST + 1];
			double hdist_bin_probs[HDIST_BINS];
			double hdist_scores[MAX_TOTAL_HDIST + 1];

			Sint8 *lookup_table1;//bit 7: valid or invalid value for a core block, bit 0: high or low value for core block
			Uint8 *lookup_table2;//hamming distance from idealized value for core block
			Uint8 *lookup_table3;//hdist -> bin

			Sint8 lookup1(Word value) const {
				using namespace PractRand::Internals;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
				if (BITS_PER_BLOCK < WORD_BITS) value &= Word((1ull << BITS_PER_BLOCK) - 1);
#pragma GCC diagnostic pop
				if (BITS_PER_BLOCK <= MAX_LOOKUP_L2) return lookup_table1[value];
				else if (BITS_PER_BLOCK <= 16) return lookup_table1[count_ones16(value)];
				else if (BITS_PER_BLOCK <= 32) return lookup_table1[count_ones32(value)];
				else return lookup_table1[count_ones64(value)];
			}
			Sint32 _lookup1(Word value) const {
				using namespace PractRand::Internals;
				if (BITS_PER_BLOCK <= MAX_LOOKUP_L2) return lookup_table1[value];
				else if (BITS_PER_BLOCK <= 16) return lookup_table1[count_ones16(value)];
				else if (BITS_PER_BLOCK <= 32) return lookup_table1[count_ones32(value)];
				else return lookup_table1[count_ones64(value)];
			}
			Uint32 _lookup2(Word value) const {
				using namespace PractRand::Internals;
				if (BITS_PER_BLOCK <= MAX_LOOKUP_L2) return lookup_table2[value];
				else if (BITS_PER_BLOCK <= 16) return lookup_table2[count_ones16(value)];
				else if (BITS_PER_BLOCK <= 32) return lookup_table2[count_ones32(value)];
				else return lookup_table2[count_ones64(value)];
			}
			void analyze_block(Word block_value, long &bucket, int bucket_bit, long &hdist) const {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
				if (BITS_PER_BLOCK < WORD_BITS) block_value &= Word((1ull << BITS_PER_BLOCK) - 1);
#pragma GCC diagnostic pop
				bucket |= (_lookup1(block_value) & 1) << bucket_bit;
				hdist += _lookup2(block_value);
			}
			int get_hdist_bin(int hdist) const;
			bool is_word_bad(Word word) const;
			bool is_core_bad(const Word *core) const;
			void core_analysis(const Word *core, int &index, int &ham) const;//only call on valid cores
			void count_ones_distribution(Word bits, Uint64 *counts, int num = WORD_BITS);
		public:
			NearSeq2();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		class NearSeq3 : public TestBaseclass {
				/*
					NEWEST PLAN:
						look for aligned blocks of 64,128,256,512,etc bits that have unusually low hamming weights
						unusual defined as 1/1K, 1/1M, 1/1B, etc
							but there is no value to (much) less than the best available rarity, is there?
						for the last ~64 such blocks, record the preceding block and following block
							actually, might as well just stop once the max is reached?
						record & test the hamming distance distribution among preceding & following blocks
							not a chi-squared (or equivalent) test - the events we're looking for will be too low likelyhood
								do something like what we did for BRank
							theoretically, using the same block more than twice biases the results, and we use them a LOT more than twice eventually
								but in practice I think we can ignore that, so long as the number of times used is less than the dimensionality, which has a minimum of 64
						ideally also checking the hamming distance distribution among center blocks, but what it should be is too hard to figure out

					DONE:
						this idea was significantly more actionable than prior ones
						still need to support larger block sizes - it's currently capping out at 512 bytes, could up that to 8KB without too much work
						but there's not a lot of point
						could also look at the number of cores meeting the tightest threshold too?

					*/
			typedef Uint64 Word;

			enum {
				BITS_PER_WORD = sizeof(Word) * 8,
				NUM_SIZES = 7,
				BUNCH_SIZE = 2048,
				DUMP_ALL_ON_THRESHOLD_CHANGE = 0,
				THRESHOLD_TIGHTENING_FRACTION = 500,//out of 1,000
			};
			//Uint32 current_thresholds[NUM_SIZES];
			//Uint32 partial_counts[NUM_SIZES];
			struct ScoringData {
				Uint32 threshold;
				Uint32 total_pairs;
				Uint16 count;
				Uint16 flushed;
				double score_hw_pre;
				double score_hw_post;
				double score_hw_core;
				double score_hw_both;
				double score_pd_pre;
				double score_pd_post;
				//double score_pd_core;
				double score_pd_both;
			};
			struct PerSizeData {
				Sint32 threshold;//inclusive
				Sint32 count;
				Sint32 max_count;//should be some function of block size... possibly quite large for block size 1
				Uint32 num_bits;
				Uint32 partial_hw;
				Uint8 size_index;
				bool even;
				//bool warmup_needed;//not used for anything smaller than 1 KB

				std::vector<ScoringData> old_thresholds;

				std::vector<Word> raw_core;
				std::vector<Word> raw_pre_;
				std::vector<Word> raw_post;
				std::vector<Sint32> hw_core;
				std::vector<Sint32> hw_pre_;
				std::vector<Sint32> hw_post;

				//double running_score_core;
				void tighten_threshold();
				void handle_sample(long weight, Word *position);
				long do_scoring(ScoringData &data, bool filter);
			};
			static double _basic_hw_scoring_function(std::vector<double> &pdf, std::vector<double> &cdf, long hw);
			PerSizeData per_size_data[NUM_SIZES];
			void handle_excess_sizes(long si, long weight, Word *ptr);

		public:
			NearSeq3();
			~NearSeq3();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			//int get_blocks_to_repeat() const;

			virtual void test_blocks(TestBlock *data, int numblocks);
		};

		/*class NLHWE : public TestBaseclass {
			typedef Uint64 Word;
			typedef Uint8 Category;
			typedef Uint16 Position;
			enum { WINDOW_SIZE = 32768, WORDS_PER_BLOCK_L2 = 1, NUM_CATS = 3 + sizeof(Word)+WORDS_PER_BLOCK_L2 + 0 };
			long warmup;
			//Word history_window[WINDOW_SIZE];
			Category history_cats[WINDOW_SIZE * ];
			Position prior_occurance[WINDOW_SIZE];
			Uint64 count_cats[NUM_CATS];
			Position last_occurance[NUM_CATS];

			Uint64 total_in_window;
			Category block_to_category(Word *block);
		public:
			NLHWE();
			//~NLHWE();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};*/

		class LHWSCN : public TestBaseclass {//low-hamming-weight sequence correlation, N-level
			typedef Uint64 Word;//do not change, it's kinda hardwired to 64 bits words
			typedef Uint16 Classification;
			typedef Uint16 Index;
			//enum { WINDOW_SIZE = 32768, WORDS_PER_BLOCK_L2 = 1, NUM_CATS = 3 + sizeof(Word)+WORDS_PER_BLOCK_L2 + 0 };
			enum {
				WORD_BITS_L2 = 6,//do not change
				BLOCK_BITS_L2 = TestBlock::SIZE_L2 + 3,
				SAMPLES_KEPT = 1024, 
				MAX_MINIMUM_LEVEL = BLOCK_BITS_L2 - WORD_BITS_L2,// 1 KB
				MAX_LEVELS = 8 + MAX_MINIMUM_LEVEL
			};
			PractRand::RNGs::LightWeight::arbee internal_rng; // round corner-cases randomly to preserve distribution
			//int lowest_level;
			struct LevelData {
				struct SampleData {
					Uint64 weight;
					Classification cat;
					Index index;
				};
				SampleData samples[SAMPLES_KEPT];
				long num_samples;
				long parity;
				long carry_data[6 + 1 + MAX_LEVELS];//replace with some kind of dynamic allocation?
				Word *storage;

				long threshold_low, threshold_high;
				std::vector<long> count_at_weight;//number of commons at each hamming weight from theshold_low to threshold_high, inclusive on both ends
				std::map<long, long> outliers;

				//blocks have hamming wieght HW such that threshold_low <= HW <= threshold_high
				//anything with hamming weight below threshold_low gets redirected to outliers ; that's probably a failure though, since any such should be REALLY rare
				//if samples overflow then threshold_high is reduced until there is at least 1/16th of SAMPLES_KEPT is available ; if that would drop threshold_high below threshold_low then autofail is set
			};
			LevelData level_data[MAX_LEVELS];
			int minimum_level;
			bool autofail;

			//void make_room(unsigned long level);
			static void sample_to_counts(unsigned long level, Word *sample, long counts[]);//returns hamming weight of entire region
			//static long large_sample_to_counts(unsigned long level, TestBlock *blocks, unsigned long blocks_available, long counts[]);//returns number of testblocks used
			static Classification counts_to_classification(unsigned long level, long counts[]);
			void handle_sample(unsigned long level, Word *sample);
			//void handle_next_level(unsigned long level, Word *sample, )
		public:
			LHWSCN(unsigned long _minimum_level);
			~LHWSCN();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			//virtual std::string get_name() const;
			//virtual void get_results(std::vector<TestResult> &results);
			//virtual int get_blocks_to_repeat() const;//number needed varies depending upon MAX_LEVELS, unless I add some kind of internal buffering or accumulation

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
	}//Tests
}//PractRand
