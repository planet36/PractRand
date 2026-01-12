namespace PractRand {
	namespace Tests {
		class CoupGap : public TestBaseclass {
			//enum {MAX_OLDEST_AGE = 256 * 16};
			enum {MAX_CURRENT_AGE = 256 * 2};

			//low8 = current_sym, high8 = oldest_sym
			FixedSizeCount<Uint8, 256*256> count_syms_by_oldest_sym;

			//low8 = oldest_sym, high8 = current_gap
//			FixedSizeCount<Uint8, 256*MAX_CURRENT_AGE> count_gaps_by_oldest_sym;

//			Uint32 last_sym_pos[256];
			unsigned long oldest_sym;
			unsigned long youngest_sym;
			Uint8 next_younger_sym[256];
			bool sym_has_appeared[256];
			int symbols_ready;

			Uint32 autofail;
			Uint64 blocks;
		public:
			CoupGap() {}

			virtual void init( RNGs::vRNG *known_good );
			virtual std::string get_name() const;
			virtual void get_results ( std::vector<TestResult> &results );
			virtual void test_blocks(TestBlock *data, int numblocks);
		};
/*
thoughts for next test:

for 32 bit values:
	ignore any that don't have 8 zero bits at the top
	keep an array of 2^24 bools keeping track of which other values have shown up
	that's 2 megabytes
	every 32 megabytes of input, swap out that block for a clear one
		we expect to see the first dup at 2^(2+(32+8)/2)= 4 MB
	keep the old block
	if we already had an old block
		bitwise-OR what would be the new old block with the previous one
		and put the result in the next level
		repeat recursively
		if we keep 20 levels, that's something like 42 megabytes
	when we're doing a bitwise-OR, also do other stuff
		popcount the new block
		and also do a bitwise-AND - don't keep the results, just popcount it
	check distributions of each type of popcount against its expected distribution
		at each level
		popcount of total ORed
		popcount of ORed up to last level then ANDed

for 64 bit values:
	ignore any that don't start with 32 zero bits at the top
	keep an array of 2^24 entries
		each has an 8-bit frequency and an 8-bit value
	32 zero bits + 24 entry index bits + 8 value bits = 64 bits
	but... that means the first *partial* duplicate should be found after...
		2^(3+32+(24/2)) = 2^47 bytes = 128 terabytes!
	NO
	is the same basic approach as at 32 bit still viable, at all?
		zero bits:		32	28	24	20	16
		base array size:	512 MB	8 GB	128 GB	2 TB	32 TB
		1st dup dist:		2048 TB	512 TB	128 TB	32 TB	8 TB
		...ridiculous
	how about if I just keep a list of past hits?
		pack the data tightly...
		zero bits:	8	12	16	18	20	24
		storage/TB:	3.5 GB	224 MB	12 MB	3 MB	768 KB	40 KB
		1st dup dist:	512 GB	2 TB	8 TB	16 TB	32 TB	128 TB
		closer to practical
		Z=16 or Z=18 looks best, for what that's worth
		
suppose I go with 64 bit alignment, but 48 bit values?
	doing the same stuff as for 32 bit values above
	zero bits:	16	18	20		(32/32/8)
	base size:	512 MB	128 MB	32 MB		2 MB
	1st dup dist:	32 GB	64 GB	128 GB		4 MB
	region size?:	64 GB	64 GB	64 GB		32 MB
	the Z=18 version looks kinda acceptable?
64 bit alignment, 56 bit values
	zero bits:	16	20	24	26	28	30
	base size:	128 GB	8 GB	512 MB	128 MB	32 MB	8 MB
	1st dup dist:	512 GB	2 TB	8 TB	16 TB	32 TB	64 TB
	the Z=28 version looks best, but is almost unusable

now suppose I go with 48/64 again, but the list of past hits approach?
	zero bits:	8	12	16	18	20	24
	storage/TB:	2.5 GB	160 MB	8 MB	2 MB	512 KB	24 KB
	1st dup dist:	2 GB	8 GB	32 GB	64 GB	128 GB	512 GB
	...could store extra metadata too, that could be useful
	best is probably Z=16, Z=18, or Z=20

now suppose I go with 56/64 with the list of past hits approach?
	zero bits:	8	12	16	18	20	24
	1st dup dist:	32 GB	128 GB	512 GB	1 TB	2 TB	8 TB
	best is probably Z=16, Z=18, or Z=20

now suppose I go with 60/64 with the list of past hits approach?
	zero bits:	8	12	16	18	20	24
	1st dup dist:	8 GB	32 GB	2 TB	4 TB	8 TB	32 TB

for now I'm thinking: the full 32 bit version @ Z=8
	AND, at extreme lengths, 64/64 history based, Z= 16 or 18

*/
	}//Tests
}//PractRand
