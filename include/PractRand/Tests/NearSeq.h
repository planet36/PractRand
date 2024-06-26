#pragma once

#include <bit>

namespace PractRand::Tests {
		class NearSeq final : public TestBaseclass {
		protected:
			typedef Uint64 Word;
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

			static constexpr int WORD_BITS = sizeof(Word)* 8;
			static constexpr int WORD_BITS_L2 = (WORD_BITS == 8) ? 3 : ((WORD_BITS == 16) ? 4 : ((WORD_BITS == 32) ? 5 : ((WORD_BITS == 64) ? 6 : -1)));
			static constexpr int NUM_BUCKETS_L2 = 9;
			static constexpr int NUM_BUCKETS = 1 << NUM_BUCKETS_L2;
			static constexpr int BITS_PER_BLOCK = 7;
			static constexpr int BLOCKS_PER_CORE = NUM_BUCKETS_L2;
			static constexpr int CORE_SEQUENCE_BITS = BLOCKS_PER_CORE * BITS_PER_BLOCK;

			static constexpr int CORE_WORDS = (CORE_SEQUENCE_BITS + WORD_BITS - 1) / WORD_BITS;
			static constexpr int EXTRA_WORDS = 2;
			static constexpr int SEQUENCE_WORDS = CORE_WORDS + EXTRA_WORDS;
			static constexpr int SEQUENCE_BITS = SEQUENCE_WORDS * WORD_BITS;
			static constexpr int SEQUENCE_WORD_OFFSET = (EXTRA_WORDS + 1) / 2;
			static constexpr int MAX_ERRORS_PER_BLOCK = 2;
			static constexpr int GOOD_ERRORS_PER_BLOCK = 1;

			static constexpr int MAX_CORE_DISTANCES = CORE_SEQUENCE_BITS / 2;//distances of this much or higher get combined with (one less than this)
			//static constexpr int AX_CORE_DISTANCES = MAX_ERRORS_PER_BLOCK * BLOCKS_PER_CORE;//distances of this much or higher get combined with (one less than this)
			//static constexpr int ARES_THRESHOLD = 8;//core distances less than or equal to this are required for rares
			//static constexpr int ARES_SIZE_L2 = 16;
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
			void init(PractRand::RNGs::vRNG *known_good) override;
			void deinit() override;
			std::string get_name() const override;
			void get_results(std::vector<TestResult> &results) override;

			void test_blocks(TestBlock *data, int numblocks) override;
		};
		class NearSeq2 final : public TestBaseclass {
		protected:
			typedef Uint32 Word;
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

			static constexpr int WORD_BITS = sizeof(Word)* 8;
			static constexpr int WORD_BITS_L2 = (WORD_BITS == 8) ? 3 : ((WORD_BITS == 16) ? 4 : ((WORD_BITS == 32) ? 5 : ((WORD_BITS == 64) ? 6 : -1)));
			static constexpr int NUM_BUCKETS_L2 = 8;
			static constexpr int NUM_BUCKETS = 1 << NUM_BUCKETS_L2;
			static constexpr int BITS_PER_BLOCK = 32;

			static constexpr int BLOCKS_PER_CORE = NUM_BUCKETS_L2;
			static constexpr int CORE_SEQUENCE_BITS = BLOCKS_PER_CORE * BITS_PER_BLOCK;

			static constexpr int CORE_WORDS = (CORE_SEQUENCE_BITS + WORD_BITS - 1) / WORD_BITS;
			static constexpr int EXTRA_FULL_WORDS = 2;//should be an even number, otherwise call it an invalid parameterization
			static constexpr int EXTRA_PARTIAL_WORD_BITS = CORE_WORDS * WORD_BITS - CORE_SEQUENCE_BITS; //should be less than 1 word, otherwise call it an invalid parameterization
			static constexpr int EXTRA_BITS = EXTRA_FULL_WORDS * WORD_BITS + EXTRA_PARTIAL_WORD_BITS;
			static constexpr int SEQUENCE_WORDS = CORE_WORDS + EXTRA_FULL_WORDS;
			static constexpr int SEQUENCE_BITS = SEQUENCE_WORDS * WORD_BITS;
			static constexpr int SEQUENCE_WORD_OFFSET = EXTRA_FULL_WORDS / 2;
			static constexpr int MAX_HDIST_PER_BLOCK = 12;
			static constexpr int MAX_TOTAL_HDIST = MAX_HDIST_PER_BLOCK * BLOCKS_PER_CORE;
			//static constexpr int MAX_CONCEIVABLE_HDIST = BLOCKS_PER_CORE * (BITS_PER_BLOCK / 2);
			static constexpr int HDIST_BINS = 16;

			static constexpr int MAX_LOOKUP_L2 = 12;//otherwise the tables get too big
			static constexpr int CHECK_VALIDITY_EARLY = 0;//if true, checks for a bad core after each block for the first word, then once per word thereafter
			struct Bucket {
				//Word ideal_sequence[SEQUENCE_WORDS];
				Uint64 core_hdist[MAX_TOTAL_HDIST + 1];
				Uint64 extra_counts[HDIST_BINS][EXTRA_BITS];
				void reset();
			};
			Bucket buckets[NUM_BUCKETS];
			Uint64 _total_cores;
			Uint64 _total_invalid_cores;

			Sint8 *lookup_table1;//bit 7: valid or invalid value for a core block, bit 0: high or low value for core block
			Uint8 *lookup_table2;//hamming distance from idealized value for core block

			Sint8 lookup1(Word value) const {
				if constexpr (BITS_PER_BLOCK < WORD_BITS) value &= (1UL << BITS_PER_BLOCK) - 1;
				if constexpr (BITS_PER_BLOCK <= MAX_LOOKUP_L2) return lookup_table1[value];
				else if constexpr (BITS_PER_BLOCK <= 16) return lookup_table1[std::popcount(value)];
				else if constexpr (BITS_PER_BLOCK <= 32) return lookup_table1[std::popcount(value)];
				else return lookup_table1[std::popcount(value)];
			}
			Sint32 _lookup1(Word value) const {
				if constexpr (BITS_PER_BLOCK <= MAX_LOOKUP_L2) return lookup_table1[value];
				else if constexpr (BITS_PER_BLOCK <= 16) return lookup_table1[std::popcount(value)];
				else if constexpr (BITS_PER_BLOCK <= 32) return lookup_table1[std::popcount(value)];
				else return lookup_table1[std::popcount(value)];
			}
			Uint32 _lookup2(Word value) const {
				if constexpr (BITS_PER_BLOCK <= MAX_LOOKUP_L2) return lookup_table2[value];
				else if constexpr (BITS_PER_BLOCK <= 16) return lookup_table2[std::popcount(value)];
				else if constexpr (BITS_PER_BLOCK <= 32) return lookup_table2[std::popcount(value)];
				else return lookup_table2[std::popcount(value)];
			}
			void analyze_block(Word block_value, long &bucket, int bucket_bit, long &hdist) const {
				block_value &= (1UL << BITS_PER_BLOCK) - 1;
				bucket |= (_lookup1(block_value) & 1) << bucket_bit;
				hdist += _lookup2(block_value);
			}
			int get_hdist_bin(int hdist) const;
			bool is_core_bad(const Word *core) const;
			void core_analysis(const Word *core, int &index, int &ham) const;//only call on valid cores
			void count_bits_distribution(Word bits, Uint64 *counts, int num = WORD_BITS);
		public:
			NearSeq2();
			void init(PractRand::RNGs::vRNG *known_good) override;
			void deinit() override;
			std::string get_name() const override;
			void get_results(std::vector<TestResult> &results) override;

			void test_blocks(TestBlock *data, int numblocks) override;
		};
}//PractRand
