namespace PractRand {
	namespace Tests {
		class Gap16 : public TestBaseclass {
		public:
			virtual void init( PractRand::RNGs::vRNG *known_good );
			virtual std::string get_name() const;// {return std::string("Gap16");}
			//virtual double get_result();
			//virtual double result_to_pvalue ( Uint64 blocks, double r );
			virtual void get_results ( std::vector<TestResult> &results );

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			enum { 
				SIZE1 = 1<<9, //handles very short lags, tracks sets of (1 << SET1_SHIFT) gaps
				SIZE2 = 1<<14, //handles medium lags, tracks each set of (1 << SET2_SHIFT) gaps
				SIZE3 = 1<<14, //handles very long, very rare lags, tracks each set of (1 << SET3_SHIFT) gaps (VERY rare cases go to extreme_lags instead)
				SET1_SHIFT = 0,//prior to v0.95, settings were: {18:1,17:2,17:3,Uint8} ; now considering {8:0,16:3,16:5,Uint8} and {8:0,14:5,14:6,Uint16}
				SET2_SHIFT = 5,
				SET3_SHIFT = 7
			};
			/*
				idea: replace those three piecewise linear gap distance counts with a pseudo floating point gap count
					9 bits of significand and 4 bits of biased exponent means 1<<13 table size, which is *very* reasonable
					and enough that denormalized small values are only up to 1<<9, while overflowing large values aren't needed until 1<<22 or so, which is plenty
					should be good so long as BSF/BSR are fast enough...
					could even do solely by table... but that's problematic since a single table lookup doesn't give enough dynamic range easily
						maybe on the table section have a denormalized region up to 1<<12
						then bits 13-20 would be the psueudo floating point
						anything 1<<20 and 1<<24 would go to a 3rd path, and anything over 1<<24 would go to extreme_lags
			*/
			FixedSizeCount<Uint16, SIZE1 + SIZE2 + SIZE3> counts;
			std::vector<Uint32> extreme_lags;
			void increment_lag(Uint32 lag);
			bool autofail;
			Uint32 last[65536];
			int warmup;
		};

		class LimitedBigGapPrototype : public TestBaseclass {
		public:
			LimitedBigGapPrototype() { last = NULL; }
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual void deinit();
			virtual std::string get_name() const;// {return std::string("LG64");}
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			enum {
				TOTAL_BITS = 64,//don't change
				TOTAL_BITS_L2 = 6,//log2 of the above
				ZERO_BITS = 4,
				INDEX_BITS = 8, //may not exceed 32 - and even that would take a ton of memory
				CHECK_BITS = TOTAL_BITS - ZERO_BITS - INDEX_BITS,
				CHECK_COUNT_SHIFT = 2,
				MAX_MATCH_LEVEL = CHECK_BITS >> CHECK_COUNT_SHIFT,
				MATCH_LEVELS_L2 = TOTAL_BITS_L2 - CHECK_COUNT_SHIFT,
				MATCH_LEVELS = 1 << MATCH_LEVELS_L2,
				NUM_DISTANCE_BUCKETS = 64,//don't change
				SIZE = 1 << INDEX_BITS,
			};
			struct Entry {
				Uint64 value;
				Uint64 position;
			};
			Entry *last;
			Entry *first;
			Uint64 counts[NUM_DISTANCE_BUCKETS * MATCH_LEVELS];
			long warmup1, warmup2;
			long distance_to_bucket(Uint64 gap_value);
			void bucket_to_distance_range(long distance_bucket, Uint64 &min_distance, Uint64 &max_distance);//inclusive on both sides
			//void handle_gap(Uint64 gap_value);
		};
		/*
			lessons from LimitedGap64Prototype:
		
			Gap tests significantly larger than 16 bits are viable with enough sacrifices.  
			But full 64 bit gaps aren't viable unless you're testing thousands of terabytes.  
			No one set of parameters from the above is good for a wide range of test lengths.  
			Fortunately, TOTAL_BITS (and thus CHECK_BITS) can be implemented for multiple values simultaneously efficiently.  
			I'm currently looking at 3 variants:
				small memory (4 to 8 megabtes) with 8 ZERO_BITS and INDEX_BITS set to 18 or 9, and TOTAL_BITS varying from 32 to 64, 
				medium memory (~128 megabytes), with 10 ZERO_BITS and INDEX_BITS set to 23 or so, and TOTAL_BITS varying from 36 to 64, 
				and large memory (2 gigabytes) with 12 ZERO_BITS and INDEX_BITS set to 27 or so, and TOTAL_BITS varying from 40 to 64.  
			The larger ones aren't as much better as you might expect, but still can reach usable 64 bit gaps faster than the smaller ones, though may be useless on test runs less than a few gigabytes.  
			Currently expected useful test length vs TOTAL_BITS vs test size:
				TOTAL_BITS		min test length for 100 gaps, given small memory		same for medium memory			same for large memory
				36								2**29 bytes								2**30 bytes						n/a
				40								2**32 bytes								2**32 bytes						2**33 bytes
				44								2**35 bytes								2**34 bytes						2**35 bytes
				48								2**39 bytes								2**37 bytes						2**37 bytes
				52								2**43 bytes								2**39 bytes						2**39 bytes
				56								2**47? bytes							2**43 bytes						2**41 bytes
				60								2**51? bytes							2**47 bytes						2**44 bytes
				64								2**55? bytes							2**51? bytes					2**47? bytes

			time has passed, later thoughts:
				1. chi-squared test results on this were kinda weird.  They could work, but it was a bit counterintuitive when they would vs when they wouoldn't.  
				2. how much of what Gap16 catches can be caught by this?
				3. tradeoffs of the check_bits approach compared to the FPF/FPG/FPM approach?  mixing the two?

				I'm gonna say INDEX_BITS * 2 + ZERO_BITS + (1 to 4), with CHECK_BITS not factoring in at this parameter range?

*/
		class LimitedBigGapPrototype2 : public TestBaseclass {
		public:
			LimitedBigGapPrototype2() {}
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual void deinit();
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			enum {
				MINIMUM_ZEROES = 4,
				MAX_EXP = 16,
				SIGNIFICAND_BITS = 8,
				NUM_DISTANCE_BUCKETS = 256,
				MATCH_LEVELS = 32
			};
			struct Entry {
				Uint64 check;
				Uint64 position;
			};
			struct Platter {
				Entry *last;
				Entry *first;
				Uint64 counts[NUM_DISTANCE_BUCKETS * MATCH_LEVELS];
			};
			long warmup1, warmup2;
			long distance_to_bucket(Uint64 gap_value);
			void bucket_to_distance_range(long distance_bucket, Uint64 &min_distance, Uint64 &max_distance);//inclusive on both sides
			//void handle_gap(Uint64 gap_value);
		};

		class GapUniversal1 : public TestBaseclass {
			/*
				Basically just the the test described in "A Universal Statistical Test for Random Bit Generators" by Ueli Maurer, appeared Journal of Cryptography, vol. 5, no. 2, 1992, pp. 89
				The only real changes are those necessary to integrate it in to PractRand, plus some disabled code in init() to calculate the constants it uses.  
				This implementation operates on 16 bit values (L=16) and discards the first 65536*10 of those in warmup (Q=655360).  

				Ultimately the algorithm itself is just a weaker variation of the classic gap test.  Of course the paper isn't just the algorithm but also arguments about what kind of bias if can detect.  
				If I'm reading it correctly, it basically say that this can find any possible bias provided that the RNG is imperfect and is equivalent to a nondeterministic state machine with no more than 65536 states.  
					...the agument might be okay, but that's not really anything I care about, unless it could be scaled up to much drastically larger state sizes (it can't)
				However, there is one detail of interest to me: compared to normal gap tests, this represents the observed gap distribution as a single real number - the mean of the logs of the gap distances.  
					Which makes it weaker than classic gap tests (which typically use a chi-squared test on different ranges of gap distances), 
					but uses far less space for the gap distribution, making possible innovations that would otherwise have been infeasible.  
			*/
		public:
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			static const double expected_log_gap;
			static const double gap_log2_variance;
			Uint64 last[65536];
			enum {WARMUP = 10 * 65536};
			double log_gap_sum;
			Uint64 num_gaps;
			long warmup;
			void handle_gap(Uint64 gap_value);
		};
		class GapUniversal2 : public TestBaseclass {
		/*
			same as GapUniversal1, but this time I added a bunch of optimizations and innovations:
				1. It now stores the position of previously observed values in 32 bits instead of 64 - some additional stuff is needed to correctly handle large gap values because of that, but it's worthwhile
					I stole this optimization from my own pre-existing gap tests variations.  
				2. It now uses multiplication of gap values instead of summing logs of gap values.  This is a big deal for performance, but requires a little extra work to keep everything working properly.  
					A fairly obvious optimization, I'd be surprised if no one else has done it, though I haven't actually seen or heard of anyone else doing it.  
				3. It no longer needs extensive warmup.  I'm not sure how useful it is on short data streams, but it at least produces sane results without discarding any data, and with a minimum useful stream length no worse than the 1 KB that PractRand treats as the minimum possible stream length.  
					Gap tests have long had a problem with requiring discarded data at the start (or end), I've tried several things in the past to fix that, but this is my most complete fix yet.  
					I think I've heard of others trying a number of things too, none of which were much good.  
					So overall I'm pretty happy with this atm, even though I don't expect much utility from it.  

			Two variations to try to improve sensitivity suggest themselves:
				1. the gap-product could be tracked on a per-value basis.  
						Straightforward
						the larger cache footprint would hurt a lot
						the data would have to be combined on shorter runs (elmininating all advantages) for testability reasons
						even on longer runs I think it would only make a difference on a few PRNGs.  
						OVERALL: bad idea
				2. the gap-product could be tracked on a per-low-bits-of-gap basis (or per-high-bits, or per-middle-bits? ... both low and high bits seperately?)
						Gathering the data is easy... not 100% sure about evaluating it
						The cache footprint wouldn't necessarily increase a lot? just 256 or so seperate gap products would probably be enough.  
						the data would have to be combined on shorter runs (elmininating all advantages) for testability reasons
						OVERALL: decent and straightforward

			If I wanted to go farther afield, some other possiblities come to mind:
				3. coupon tests on regions of gap values
						Collecting the data wouldn't be too hard, but figuring out the correct ratios at the end in any non-empiric fashion would be... difficult
						essentially all gap values are less then 2^23, all coupon bitfields combined should fit in to 1 MB, and the frequently accessed portions should fit in to 128 KB
						if we do coupon tests for sets of 256 gap values then the counts of completed sets is reasonable, as is the counts of hits on each region
						and we can collect gap-products for each region too
						and pre-region frequences to do a more direct gap test on
						any disturbance in the distribution of gap values should register on at least one of: the coupon tests, the gap product tests, or the frequency tests
						likely more, but at least one should catch it nearly as quickly as reasonably possible, at least from a gap perspective
						how fast is it to collect all the data?
							test that gap value isn't insane (a quick, always-predicted branch)
							bitshift to find the region, increment the region hit count, multiply the region gap product (with expected inverse?), normalize the gap product if the region hit count is a multiple of 64
							update the bitfield, check if that word of bitfield was completed, if so check the rest of its region, if it's all complete update coupon statistics and reset
							looks pretty decent overall?
						and how easy is it to evaluate the data?
							the region hit counts are easy to evaluate, and the gap products aren't hard either
							but the coupon completion counts... calculating exact predictions there seems impossible
							the good news is that the completion to hits ratio should be identical for all regions
							so if I can get remotely decent empirical data or theoretical estimate for one region, then I can use it for all of them
							try: mask out all but the lowest 8 bits of gap, run ~64 TB of that, resulting ratio should be good enough for ... runs up to 8192 TB with the masking removed?)
							but short runs will need special handling - once the expected number of hits reaches ~10 it assymptotically matches one of the well studied distribution, but when the expected number is less than 2 the ratio gets weird
							so also try: with the masking, do numerous short runs to determine at least roughly what that curve looks like - unfortunately, more than just powers of 2 will be needed since the data must be applicable to all regions, so the hits won't conform to a neat picture
						OVERALL: possibly excessive?

				If integrating in to FPG/FPM, #2 would work well, but the bitfields for #3 would be utterly impossible.  
		*/
		public:
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;// {return std::string("Gap16");}
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			enum { DUMMY_VALUE32 = 0xFFffFFfful };
			static const double expected_inverse;
			static const double gap_log2_variance;
			Uint32 last[65536];
			double gap_product;
			long gap_product_extracted_L2;
			Uint64 num_gaps;
			int warmup;
			bool autofail;
			Uint32 first[65536];
			Uint8 state[65536];
			void normalize();
			void handle_gap(Uint32 gap_value);
		};

		class GapUniversal3 : public TestBaseclass {
			/*
				Like the above GapUniversal2, but with variation #2 from its comments applied.  
			*/
		public:
			virtual void init(PractRand::RNGs::vRNG *known_good);
			virtual std::string get_name() const;// {return std::string("Gap16");}
			virtual void get_results(std::vector<TestResult> &results);

			virtual void test_blocks(TestBlock *data, int numblocks);
		protected:
			enum { DUMMY_VALUE32 = 0xFFffFFfful };
			enum { LOW_BITS = 8 };
			Uint32 last[65536];
			struct LBI {
				double gap_product;
				double expected_gap_inverse;
				long gap_product_extracted_L2;
				Uint64 hits;
			} low_bit_indexed_data[1 << LOW_BITS];
			Uint64 count_pow2;
			double freq_probs[1 << LOW_BITS];
			double mean_lg[1 << LOW_BITS];
			double variances[1 << LOW_BITS];
			double overall_mean;
			double overall_variance;
			void normalize(LBI &lbi);
			void handle_gap(long gap_value);
			int warmup;
			bool autofail;
			Uint32 first[65536];
			Uint8 state[65536];
		};
	}//Tests
}//PractRand
