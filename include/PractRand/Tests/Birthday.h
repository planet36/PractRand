namespace PractRand {
	namespace Tests {

/*
Birthday tests
	originally I had little respect for the Birthday Test, as it looks kind of awful, and most implementations (on 32 bit integers) are kind of awful, and scaling up quality seemed difficult
		however, on finding a paper (I think it was by the TestU01 authors) suggesting that it works best at lambda~=1, I became more interested
			the results in that paper looked kind of awesome, and testing bore that out somewhat when using large sort buffers
			unfortunately, large sort buffers are... problematic
				PractRand perfers to assume CPU time (and to a lesser extent, cache footprint) are the constraining factor in testing
					when users don't have fast enough CPUs, they can simply run tests longer ; when quality isn't high enough, they can simply run tests MUCH longer
				when the constraining factor is memory, scaling tests becomes impossible
					well, much more problematic anyway
				also, PractRand makes some choices that don't go well with memory constraints
					for instance, the "foldings" means that PractRand typically has 3 to 4 copies of each test in memory simultaneously
	lambda=1 (expected number of duplicate deltas of about 1.0... or at least keeping it between 0.5 and 4.0) means roughly a buffer size proportional to the cube root of the maximum word value
	something like a 2^11 word buffer for 32 bit words, 2^21 word buffer for 64 bit words, 2^43 word buffer for 128 bit words, etc
		which is usually impractical for 128 bit words - that would be 128 terabytes assuming an in-place sorting algorithm was used
	intermediate word sizes are fine - just mask out some bits ; no need to fix up the alignments or anything, ignoring some bits is fine
	quality of test seems to scale quite decently with sort-buffer size, so long as word sizes are adjusted to keep lambda at least faintly near 1
	the actual test involves sorting the buffers, looking at the differences between adjacent sorted elements, 
		and then sorting the differences to check for duplicates
	BirthdaySystematic adds the minor innovation of permitting preliminary results to be accessed even before the buffer is first filled
		this is done by masking out extra bits, so the word size is effectively smaller
		additionally, this means that the buffer is partially sorted afterwards, so ideally the later sortings should be faster
		and this means that the buffer has plenty of unused space in it, which the deltas can be stored in temporarily
		provided that preliminary results aren't needed after the buffer is more than half full, but before it's completed
	additional innovations are required
		since masking out bits and then looking for duplicates seems wasteful, I thought I'd leave as many bits in as possible
			and looks at the spacings between deltas instead of the duplicates of deltas (log and exp of delta between deltas)
			that's what I did in BirthdayAlt, but unfortunately, this did not seem to work well at all
			could try it by looking at the distrubtion of the number of bits that match between two adjacent sorted deltas
			but I was disappointed enough with BirthdayAlternative that I haven't felt up to it
		look for duplicates between different sets of deltas
			not entirely clear that this would be helpful, even if it does kinda approximate a larger sort buffer
			and the sizes of the sets of deltas are too large for this to be any more practical that larger sort buffers without alternative means of storing the deltas
			I considered storing only deltas that had shown up as duplicates, and looking for more matches of those
				but in preliminary testing it didn't look very good
			average delta value should be (2^word_bits / buffer_size)
				that puts it at something like 2^21 for 32 bit words, 2^43 for 64 bit words, and 2^85 for 128 bit words
				if the distribution is fairly tight then, instead of making a list of deltas and sorting them to look for duplicates
					we could make a bitmask of them (maybe plus a list of outliers)
						though it would have to be a VERY narrow distribution at 2^85 or even 2^43
						probably not practical for useful word sizes under normal circumstances
					we could make a hashtable of them... actually... that sounds promising
						index table by, say, bits 0-15 of deltas, fill table with bits 16-79 of deltas
							assuming expected deltas are large, and no weird structures in the fine details...
							it would probably be worse than simply using a large sort buffer
							but it might scale better when looking for matching deltas between multiple sort buffers
							or 
						number of times when overwriting the old value with the new changes nothing is the result
		filter input - only buffer lowest 0.1% of values (those that would end up at the start of the sorted buffer, so we can actually calculate the deltas we would have seen correctly... some of them anyway)
			possibly even increase the filtering over time, after we reach lambda=1 a few times or something
				that would let us mask out fewer bits over time while keeping the same buffer size
			would it actually help anything though?
				we can adjust bits discarded to keep lambda=1
				adjust filters to keep the buffer size constant as fewer bits are discarded on the low end
				but fundamentally, is quality tied to the buffer size, the number of deltas we can compare?
			possibly run at multiple filtering levels... actually that does sound good, provided that they're not too close
				say, one at 2^8, one at 2^16, one at 2^24, and a final one at 2^32
					going as high as 2^32 means that there would be useful new results coming in no matter how far
				an unfiltered version may or may not be useful
					skipping that would improve speed a lot, but unfiltered is the most useful at least early on
					and we're not actually sure filtered versions are any good yet
				could handle the higher filter versions only once the lower filter versions are sorted
					that has the advantage of being faster, and minor disadvantage of delaying higher level results
					... but it's probably kind of pointless, given that higher filterings are so rarely hit anyway
					would probably only be worthwhile if the filtering levels weren't very different
					but in that case, we'd run out of memory from having too many sort buffers
		p-value granulaity is a problem.  At lambda=1, the granularity is very poor, but it doesn't really improve that much until vastly higher lambda values.  
			but requiring high lambda generally means reducing sensitivity (relative to data tested)
			on long enough sample sizes it will work out fine, assuming memory is limited, but that's not the usual case

later:
	I may need a variant for the core test set - certain PRNGs are broken easily by BDay tests and poorly by other things
		for instance, some NP2CLCGs like np2clcg2_53 are failed by SmallCrush, yet do fairly well on PractRand
	However, given the need to keep memory usage under control, it would end up being a weak, specialized, test.  
	Probably it would have to skip a lot of input.  Comparable to BRank, TMF, or mod3.  
	It might be best to switch back and forth between a 32 bit BDay and a 64 bit BDay.  
	A single set for 32-bit is... ~8 KB, while a single set for 64-bit is... ~16 MB.  I *could* do a 48-bit intermediate, with ~512 KB set size.  
	I couldn't really go past 64-bit while keeping memory usage under control, unless it was by filtering.  Maybe 80 bit in a pinch, or 96 bit with a super-large dataset.  
		Could I do that while keeping CPU usage minimal?  They're naturally low on CPU usage, but not *that* low, right?  I should check.  
	I could do something like 1 32-bit BDay set per 512KB, 1 48-bit per 32 MB, 1 64-bit per 2 GB, that would work
	but... how about a fixed size buffer of 128-bit samples, gradually increasing bits used as it fills up, and once full gradually increasing filtering
		should be at least 2^21 samples, 32 MB, preferably larger, 2 or 4 times that
		never actually using more than 64 bits, the 128-bit samples are only for filtering
		if that's slow-ish, then don't test all data, only every N kilobytes, where N is a parameter to the test
	progressive filtering:
		there's no filtering until the buffer fills up all the way the first time
		pick the largest power of 2 smaller than the largest value in the buffer
		that's the new filter cap - anything >= that gets discarded
		immediately after tightening the filter, the buffer will drop to about half full (autofail if it drops way too little or way too much)
*/
		class Birthday32 : public TestBaseclass {
			enum {
				BUFFER_SIZE_L2 = 12, // must be at least 8... no, 12, and no more than 30
				BUFFER_SIZE = 1 << BUFFER_SIZE_L2,
				MAX_DUPLICATES = 32
			};
			Uint32 buffer[1 << BUFFER_SIZE_L2];
			Uint64 counts[MAX_DUPLICATES];
			int num_buffered;
			void flush_buffer();
		public:
			Birthday32();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		class Birthday64 : public TestBaseclass {
			enum {
				BUFFER_SIZE_L2 = 23, // must be at least 7
				BUFFER_SIZE = 1 << BUFFER_SIZE_L2,
				MAX_DUPLICATES = 64,
				SORT_HELPER_BITS = 10
			};
			Uint64 buffer[1 << BUFFER_SIZE_L2];
			Uint64 counts[MAX_DUPLICATES];
			static void _histogram_in_place_sort64(Uint64 *base, long length, long bits_already, Uint32 freq_counts[1 << SORT_HELPER_BITS]);
			static void _histogram_in_place_sort64(Uint64 *base, long length);
			static void _histogram_sort64(Uint64 *base, long length, long bits_already, Uint32 freq_counts[1 << SORT_HELPER_BITS]);
			static void _histogram_sort64(Uint64 *base, long length);
			static void _weird_sort64(Uint64 *base, long length);
			int num_buffered;
			void flush_buffer();
		public:
			Birthday64();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		namespace BirthdayHelpers {
			enum { SORT_HELPER_BITS = 8 };
			struct i128 {
				Uint64 low;
				Uint64 high;
				bool operator==(const i128 &other) const {
					return high == other.high && low == other.low;
				}
				bool operator<(const i128 &other) const {
					if (high < other.high) return true;
					if (high > other.high) return false;
					return low < other.low;
				}
				i128 operator-(const i128 &other) const {
					i128 rv;
					rv.high = high - other.high;
					rv.low = low - other.low;
					rv.high -= other.low > low ? 1 : 0;
					return rv;
				}
			};

			// this is the fastest in-place sort I've tried so far
			// in-place is important if I want to sort something huge without assuming I can allocate a comparable amount of memory to help with the sorting
			void histogram_in_place_sort128(i128 *base, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void histogram_in_place_sort128(i128 *base, Uint64 length, long bits_already = 0);

			// ..but sometimes I'm sorting smaller buffers, with pre-allocated regions to sort into..., so maybe another interface for that would help
			void histogram_sort_and_copy(i128 *base, i128 *dest, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void histogram_sort_and_copy(i128 *base, i128 *dest, Uint64 length, long bits_already = 0);
			//possibly faster algorithm for the same interface?
			void radix_sort_and_copy(i128 *base, i128 *dest, Uint64 length, long bits_already = 0);

			void _sorted_deltas_of_sorted_values(i128 *base, long length_L2, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void _sorted_deltas_of_sorted_values(i128 *base, long length_L2);

			void histogram_in_place_sort64(Uint64 *buffer, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void histogram_in_place_sort64(Uint64 *buffer, Uint64 length, long bigs_alreaady = 0);
			void histogram_sort_and_copy64(Uint64 *buffer, Uint64 *dest, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			void histogram_sort_and_copy64(Uint64 *buffer, Uint64 *dest, Uint64 length, long bigs_alreaady = 0);
			void radix_sort_and_copy64(Uint64 *buffer, Uint64 *dest, Uint64 length, long bits_already = 0);

			bool count_low16_duplicates(Uint64 *buffer, Uint64 length, Uint64 &count);// untested and not used by anything yet ; returns true if every value was less than 2**16

			//void _sorted_deltas_of_sorted_values(Uint64 *buffer, long length_L2, Uint64 freq_counts[1 << SORT_HELPER_BITS]);
			//void _sorted_deltas_of_sorted_values(Uint64 *buffer, long length_L2);

			//void _histogram_in_place_sort64(Uint64 *base, long length, long bits_already, Uint32 freq_counts[1 << SORT_HELPER_BITS]);
			//void _histogram_in_place_sort64(Uint64 *base, long length);
			//void _histogram_sort64(Uint64 *base, long length, long bits_already, Uint32 freq_counts[1 << SORT_HELPER_BITS]);
			//void _histogram_sort64(Uint64 *base, long length);
			//void _weird_sort64(Uint64 *base, long length);

		};
		class BirthdayLambda1 : public TestBaseclass {
		protected:
			//optimized for lambda=1, few runs, as described in "On the performance of birthday spacings tests with certain families of random number generators" (L'ecuyer & Simard, 2001)
			typedef BirthdayHelpers::i128 i128;
			enum {
				SORT_HELPER_BITS = BirthdayHelpers::SORT_HELPER_BITS,
				DO_LARGEST_SPACING = 1,
			};
			bool autofail;
			Uint64 sort_helper_counts[1 << SORT_HELPER_BITS];
			//i128 buffer[1 << BUFFER_SIZE_L2];//can't have arrays this large inside a class due to object file or executable file format constraints
			//std::vector<i128> buffer;// ... and the STL vector implementation I'm using throws some kind of exception if it exceeds about 4 GB or so
			i128 *buffer;
			Uint64 num_buffered;
			virtual Uint64 flush_buffer();
			double duplicates;
			double expected_duplicates;
			double longest_spacing;
			int buffer_size_L2;
			int bits_to_use;
		public:
			BirthdayLambda1(int buffer_size_L2_ = 26);
			virtual ~BirthdayLambda1();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		/*
		class BirthdaySystematic64 {
			// does Birthday Spacings test on 64 bit integers
			// if a result is requested mid-buffer, the buffer will be sorted and processed but the contents left reusable and the results treated as ephemeral
			// generally, the buffer size is aimed to have a lambda value of about 1 when using all 64 bits
			// multiple tests are combined by simply adding their observed duplicated and expected duplicates
			// currently undecided on whether or not early use of the buffer will suppress some bits or not
			enum {
				BUFSIZE_L2 = 22, 
				BUFSIZE = 1 << BUFSIZE_L2
			};
			Uint64 buffer[BUFSIZE];
			Uint64 elements_buffered;
			Uint64 num_sorted; // this many elements in the buffer are already sorted, starting at the beginning, potentially allowing optimization to the final sorting of its contents
			Uint64 observed_duplicates;
			double expected_duplicates;
			Uint64 evaluate_buffer();
		public:
			;
		};*/
		class BirthdaySystematic128 : public BirthdayLambda1 {
			// similar to BirthdayLambda1 above
			// but if a result is requested before the first sample is ready, it will return a result for a partial buffer
			// and attempts to have everything optimized for the possibility of that partial-buffer case
			virtual Uint64 flush_buffer();
			static Uint64 get_target_num_at_bufsize(int bufsize_L2_);
			int already_sorted;//if this is half of (1ull << bufsize_L2) then incomplete_duplicates should hold 

			double score;//for scoring method 2
			static double evaluate_score(double lambda, Uint64 duplicates);

			void do_incomplete_buffer();
			double incomplete_duplicates;
			double incomplete_expected_duplicates;
		public:
			BirthdaySystematic128(int max_bufsize_L2_ = 28);
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void test_blocks(TestBlock *data, int numblocks);
		};
		class BirthdayAlt : public TestBaseclass {
			//as for BirthdayLambda1, but: 
			// keep all bits regardless of buffer size, just count an adjusting range of near deltas as if they were exact matches (or score them based upon how exact they are?)
			// try filtering the initial samples range, as if it was a small part of a larger sort buffer
			typedef BirthdayHelpers::i128 i128;
			enum {
				SORT_HELPER_BITS = BirthdayHelpers::SORT_HELPER_BITS,
				//FILTER_BITS = 0,
			};
			//i128 buffer[1 << BUFFER_SIZE_L2];
			//std::vector<i128> buffer;
			i128 *buffer;
			int num_buffered;
			int buffer_size_L2;
			int filter_bits;
			Uint64 sort_helper_counts[1 << SORT_HELPER_BITS];
			bool autofail;
			void flush_buffer();

			double score_sum_log;
			double score_sum_log2;
			double score_sum_log_sqr;
			Uint64 count;
			static void _lookup_constants(int buffer_size_L2, long double *offset, long double *deviation, long double *samples);
		public:
			BirthdayAlt(int buffer_size_L2_, int filter_bits_ = 0);
			~BirthdayAlt();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			//virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		};//*/
		class BirthdaySystematic2 : public TestBaseclass {
			// Birthday test
			// on short tests, it reduces bits tested (first it refuses to return results at all, then 19 bits, then 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 61, and finally 64)
			// on longer tests it's a 64 bit birthday tests that adds progressively tighter filtering (no filtering, then only if bit65 is 0, then only if bit66 is also 0, etc)
			// minimum usable test length varies depending upon settings, might be around 1 MB
			// this one actually works, in the standard usage scenario for PractRand, it's going in the core set (at a low priority, but still)
			// also, this finally confirms that filtering works as intended, producing solid-looking results:
			//
			//	skip_level=0	
			//	BUFFER_SIZE_L2		23		23		25		25
			//	DISABLE_FILTERING	0		1		0		1
			//	LAMBDA_PER_LEVEL	128		128		128		128
			//	size (in MB)		64		64		256		256
			//  np2clcg2_60			33		33		29		30
			//	np2clcg2_64			38		39		35		35
			//	np2clcg2_67			41		>44		39		>43
			//
			// filtering also speeds things up, eventually - it's about 700 seconds per terabyte assymptotically, plus about 60 seconds of startup overhead, in my testing
			// while unfiltered it's maybe 17 times that assymptotically
			// the speed gains take quite a while to fully take effect though
			enum {
				//BUFFER_SIZE_L2 = 24, // minimum 22, max 25, at the moment ; 22 is a 32 MB buffer, 23 is a 64 MB buffer, 24 is 128 MB, 25 is 256 MB
				//BUFFER_SIZE = 1 << BUFFER_SIZE_L2,
				TEST_BITS_MIN = 18,
				TEST_BITS_STEP = 1, // probably leave at 1?
				VALUES_PER_BLOCK = 64, //max 64, could reduce to 16 or 32 as a substitute for increasing skip_level, but that's not recommended
				LAMBDA_PER_LEVEL = 128,//don't increase above 128, the aproximation as poisson distribution breaks down if it's too high ; higher values are better for p-value granularity and should be slightly worse for actual test sensitivity
				DISABLE_FILTERING = 0
			};

			Uint64 *buffer;
			Uint64 num_buffered;
			Uint64 num_sorted;
			Uint64 buffer_size_L2;
			Uint64 total_flushes;
			Uint64 filter_check_mask;//any value that is non-zero when the upper 64 bits get anded by this gets discarded
			long filter_bits;

			double cur_lambda;
			Uint64 cur_duplicates;
			Uint64 cur_bday_bits;
			double old_lambda;
			Uint64 old_duplicates;

			Uint64 skips_left;
			const Uint64 skip_level;
			bool autofail;//currently there is no actual way to trigger this

			static double calculate_poisson_suspicion(double lambda, Uint64 duplicates);
			void tighten_filter();//called by flush_buffer, when cur_lambda gets too high
			void flush_buffer();//called when buffer is full
			void update_partial_buffer_results();//called only on the very first buffer, before it fills
		public:
			BirthdaySystematic2(Uint8 skip_level, Uint8 buffer_size_L2_);
			~BirthdaySystematic2();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void test_blocks(TestBlock *data, int numblocks);
		};
	/*	class BirthdaySystematic3 : public TestBaseclass {
			//	previous version worked well in a lot of ways, but still had some problems
			//		The big one was that some failures WENT AWAY as the bits got larger
			//			actually, is this really happening, or is it just not handling some lambda levels properly?
			//		A smaller issue was the switchover from bit reduction to filtering was kind of hardwired to happen around 1 full buffer, which created some problems
			//  so... intended changes:
			//		birthday sizes are restricted to a smaller variety - ONLY multiples of 8 bits (16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96)
			//		separate p-values for each bday size
			//		old sizes get revisited periodically like BRank
			//		support bit size reduction even on full blocks
			//		may have to have blocks tested for multiple different sizes too...
			//		?autowindowing of results?
			enum {
				//BUFFER_SIZE_L2 = 24, // minimum 22, max 25, at the moment ; 22 is a 32 MB buffer, 23 is a 64 MB buffer, 24 is 128 MB, 25 is 256 MB
				//BUFFER_SIZE = 1 << BUFFER_SIZE_L2,
				TEST_BITS_MIN = 16,
				TEST_BITS_STEP = 8, 
				VALUES_PER_BLOCK = 64, //max 64, could reduce to 16 or 32 as a substitute for increasing skip_level, but that's not recommended
				MAX_LAMBDA_PER_BLOCK = 512,
				DISABLE_FILTERING = 0
			};

			Uint64 *buffer;
			Uint64 num_buffered;
			Uint64 num_sorted;
			Uint64 buffer_size_L2;
			Uint64 total_flushes;
			Uint64 filter_check_mask;//any value that is non-zero when the upper 64 bits get anded by this gets discarded
			long filter_bits;

			double cur_lambda;
			Uint64 cur_duplicates;
			Uint64 cur_bday_bits;
			double old_lambda;
			Uint64 old_duplicates;

			Uint64 skips_left;
			const Uint64 skip_level;
			bool autofail;//currently there is no actual way to trigger this

			static double calculate_poisson_suspicion(double lambda, Uint64 duplicates);
			void tighten_filter();//called by flush_buffer, when cur_lambda gets too high
			void flush_buffer();//called when buffer is full
			void update_partial_buffer_results();//called only on the very first buffer, before it fills
		public:
			BirthdaySystematic3(Uint8 skip_level, Uint8 buffer_size_L2_);
			~BirthdaySystematic3();
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);
			virtual void test_blocks(TestBlock *data, int numblocks);
		};*/
	}//Tests
}//PractRand
