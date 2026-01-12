namespace PractRand {
	namespace Internals { namespace XBG_helpers { class BitMatrix; } }
	namespace Tests {
		//class BitMatrix;


		class BRank : public TestBaseclass {
		public:
			BRank(
				Uint32 rate_hl2_ // 2 * log2(time units per KB)
				);
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void deinit();

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			//configuration:
			Uint32 rate_hl2;
			Uint64 rate;

			//partially complete matrix:
			void finish_matrix();
			Internals::XBG_helpers::BitMatrix *in_progress;
			Uint64 blocks_in_progress;
			Uint64 blocks_needed;

			//scoring:
			enum { BASE_RANKS = 16 };
			double base_probs[BASE_RANKS];
			double base_scores[BASE_RANKS];
			double get_score(long relative_rank);
			Uint64 rank_counts[BASE_RANKS];
			double rank_expected;
			double rank_deviation;
			Uint64 total_count;//total number of matrices so far
			double total_score;//
			Uint64 total_rank;//sum of binary ranks of all matrices so far
			Sint64 first_fail_size;

			//waiting stuff:
			Uint64 time_for_size(Uint64 matrix_size) const;
			Uint64 calculate_size_step(Uint64 matrix_size) const;
			Sint64 time_till_next_matrix;
			Uint64 current_size, size_step, last_size;
		};
		class BRank_old : public TestBaseclass {
		public:
			BRank_old(
				Uint32 rate_hl2_ // 2 * log2(time units per KB)
			);
			virtual void init( PractRand::RNGs::vRNG *known_good );
			virtual std::string get_name() const;
			virtual void get_results ( std::vector<TestResult> &results );
			virtual void deinit();

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			Uint32 rate_hl2;
			Uint64 rate;
			
			void pick_next_size();
			void finish_matrix();

			Uint64 saved_time;

			//partially complete matrix:
			Internals::XBG_helpers::BitMatrix *in_progress;
			Uint32 blocks_in_progress;
			int size_index;//which PerSize is active atm?

			//stats:
			class PerSize {
			public:
				//PerSize() {}
				Uint32 size;
				Uint64 time_per;
				Uint64 total;
				enum { NUM_COUNTS = 10, MAX_OUTLIERS = 100 };
				Uint64 counts[NUM_COUNTS];

				Uint64 outliers_overflow;
				std::vector<Uint32> outliers;
				void reset();
			};
			std::vector<PerSize> ps;
		};
	}//Tests
}//PractRand
