#ifndef _practrand_rng_basics_h
#define _practrand_rng_basics_h

#include <string>
#ifndef __PRACTRAND_CONFIG_H__
#include "PractRand/config.h"
#endif //__PRACTRAND_CONFIG_H__

namespace PractRand {
	extern const char *version_str;//like "0.91", for PractRand 0.91

	bool initialize_PractRand(); //returns true normally
	//will return false if it failed to find a good source of entropy
	//  in which case the autoseeding mechanism may have trouble
	//  Usually not catastrophic, but some programs might want to abort if that happens.  
	//NOTE: initialize_PractRand() is NOT threadsafe, it should be called before threads get spun off
	
	void self_test_PractRand();
	void issue_error(const char *msg = 0);//PractRand calls this any time there is an internal error
	void hook_error_handler(void(*callback)(const char *));//this can be used to replace the default behavior of issue_error

	class StateWalkingObject;

	class SEED_AUTO_TYPE {};
	class SEED_NONE_TYPE {};
	extern SEED_AUTO_TYPE SEED_AUTO;
	extern SEED_NONE_TYPE SEED_NONE;

	namespace RNGs {
		class vRNG {
		public:
		//constuctors, destructors, seeding, serialization, & low level state manipulation:
			//vRNG(Uint64 seed1, Uint64 seed2 = 0) {seed(seed1, seed2);}
			//vRNG(vRNG *rng) {seed(rng);}
			//vRNG(SEED_AUTO_TYPE *) {autoseed();}
			//vRNG(SEED_NONE_TYPE *) {}
			virtual ~vRNG();
			virtual void seed(Uint64 seed_low, Uint64 seed_high = 0); // seed from 128 bits is the new standard - 
			// if there are fewer than 2**128 valid seeds, then every 128 bit seed from 0 up to 1 less than the number of valid seeds should act as a distinct valid seed
			virtual void seed_fast(Uint64 seed_value);//fast but low-quality seeding, if that is supported ; otherwise, just a wrapper for the normal seeding function
			virtual void seed(vRNG *rng); // by default, this uses the seeder vRNG to generate a 128 bit integer seed, then passes those 128 bits to the seedee's normal seed() method
			// individual PRNGs can also define algorithm-specific seed functions with different parameter lists, but all should support these methods
			void autoseed();
			long serialize( char *buffer, long buffer_size );//returns serialized size, or zero on failure
			char *serialize( size_t *size );//returns malloced block, or NULL on error, sets *size to size of block
			bool deserialize( const char *buffer, long size );//returns true on success, false on failure
			std::string print_state();//returns RNG state as a comma-delimited sequence of numbers
			virtual void walk_state(StateWalkingObject *) = 0;

#if defined __SIZEOF_INT128__
			// if 128 bit integers are supported, then these might be more convenient equivalents to the above
			void seed128(unsigned __int128 seedval) {seed(Uint64(seedval >> 0), Uint64(seedval >> 64));}
#endif


		//raw random bits
			virtual Uint8  raw8 () = 0;
			virtual Uint16 raw16() = 0;
			virtual Uint32 raw32() = 0;
			virtual Uint64 raw64() = 0;
			//virtual void raw_N(Uint8 *, size_t length) = 0;
#if defined __SIZEOF_INT128__
			// if 128 bit integers are supported, then this might be a convenient wrapper for raw64:
#if defined PRACTRAND_TARGET_IS_LITTLE_ENDIAN
			unsigned __int128 raw128() { Uint64 low = raw64(); Uint64 high = raw64(); return (static_cast<unsigned __int128>(high) << 64) | low; }
#elif defined PRACTRAND_TARGET_IS_BIG_ENDIAN
			unsigned __int128 raw128() { Uint64 low = raw64(); Uint64 high = raw64(); return (static_cast<unsigned __int128>(low) << 64) | high; }
#endif
#endif

		//uniform distributions
			Uint32 rand_i32(Uint32 max); // [0,max)
			Uint32 rand_i32_fast(Uint32 max); // [0,max) - faster but biased for large values of max ; do not use on platforms that lack fast 32x32->64 integer multiplication
			Uint64 rand_i64(Uint64 max); // [0,max)
			//Uint32 rand_i64_fast(Uint32 max); // not yet implemented - [0,max) - faster but biased for large values of max ; do not use on platforms that lack fast 64x64->128 integer multiplication

			float rand_float(); // [0,1)
			double rand_double(); // [0,1)
			float rand_float(float max) { return rand_float() * max; } // [0,max)
			double rand_double(double max) {return rand_double() * max;}// [0,max)

		//non-uniform distributions
			double gaussian();//mean 0.0, stddev 1.0
			double gaussian(double mean, double stddev) { return gaussian() * stddev + mean; }

		//metadata functions
			virtual Uint64 get_flags() const;
			virtual std::string get_name() const = 0;
			virtual int get_native_output_size() const = 0;//generally 8, 16, 32, 64, or -1 (unknown)
			virtual void get_maximum_seed(Uint64 &seed_low, Uint64 &seed_high) const;//returns the highest seed to produce a seeded state distinct from that of all lower seed values
			//if seed_low and seed_high are both set to zero, that means that the highest distinct seed is unknown (should only happen for non-recommended PRNGs)
			//good quality PRNGs should both set both to 0xFFFFFFFFFFFFFFFFull (meaning that all 2**128 possible seeds produce distinct seeded states)
#if defined __SIZEOF_INT128__
			// if 128 bit integers are supported, then these might be more convenient equivalents to the above
			unsigned __int128 get_maximum_seed128() const { Uint64 seed_low, seed_high; get_maximum_seed(seed_low, seed_high); return (static_cast<unsigned __int128>(seed_high) << 64) | seed_low; }
#endif

		//exotic methods (not supported by many implementations - check flags to see if they support it):
		//exotic methods 1: random access
			virtual void seek_forward(Uint64 how_far_low64, Uint64 how_far_high64=0);
			virtual void seek_backward(Uint64 how_far_low64, Uint64 how_far_high64=0);
#if defined __SIZEOF_INT128__
			void seek_forward128(unsigned __int128 how_far) {Uint64 low = Uint64(how_far); Uint64 high = Uint64(how_far >> 64); seek_forward(low, high);}
			void seek_backward128(unsigned __int128 how_far) {Uint64 low = Uint64(how_far); Uint64 high = Uint64(how_far >> 64); seek_backward(low, high);}
#endif

		//exotic methods 2: entropy pooling
			virtual void reset_entropy();//returns an entropy pool to its default state
			virtual void add_entropy8 (Uint8 );
			virtual void add_entropy16(Uint16);
			virtual void add_entropy32(Uint32);
			virtual void add_entropy64(Uint64);
			//note that "add_entropy_N(&byte_buffer[0], 13)" will typically NOT produce the same state transition 
			//  as "add_entropy_N(&byte_buffer[0], 7);add_entropy_N(&byte_buffer[7], 6);"
			virtual void add_entropy_N(const void *, size_t length);//beware of endianness issues
#if defined __SIZEOF_INT128__
#if defined PRACTRAND_TARGET_IS_LITTLE_ENDIAN
			void add_entropy128(unsigned __int128 data) { Uint64 low = Uint64(data); Uint64 high = Uint64(data >> 64); add_entropy64(low); add_entropy64(high); }
#elif defined PRACTRAND_TARGET_IS_BIG_ENDIAN
			void add_entropy128(unsigned __int128 data) { Uint64 low = Uint64(data); Uint64 high = Uint64(data >> 64); add_entropy64(high); add_entropy64(low); }
#endif
#endif

			//add_entropy_automatically returns true if a good amount (>= 128 bits) of entropy was added
			//the milliseconds parameter is the maximum amount of time it is allowed to block while waiting for entropy
			virtual bool add_entropy_automatically(int milliseconds = 0);

			virtual void flush_buffers();// some entropy pooling PRNGs have internal buffers that need to be flushed before inputs can effect outputs - this makes sure that all prior input will impact all future output

		// C++2011 <random> compatibility:
#if defined PRACTRAND_BOOST_COMPATIBILITY
			typedef Uint64 result_type;
			result_type operator()() {return raw64();}
			static const bool has_fixed_value = true;
			static const result_type min_value = 0;
			static const result_type max_value = ~(result_type)0;
			result_type min() const {return min_value;}
			result_type max() const {return max_value;}
#endif
		};
		class vRNG8 : public vRNG {
		public:
			enum {OUTPUT_BITS = 8};
			virtual Uint16 raw16();
			virtual Uint32 raw32();
			virtual Uint64 raw64();
			virtual int get_native_output_size() const;
			virtual void add_entropy16(Uint16);
			virtual void add_entropy32(Uint32);
			virtual void add_entropy64(Uint64);
		};
		class vRNG16 : public vRNG {
		public:
			enum {OUTPUT_BITS = 16};
			virtual Uint8  raw8 ();
			virtual Uint32 raw32();
			virtual Uint64 raw64();
			virtual int get_native_output_size() const;
			virtual void add_entropy8 (Uint8 );
			virtual void add_entropy32(Uint32);
			virtual void add_entropy64(Uint64);
		};
		class vRNG32 : public vRNG {
		public:
			enum {OUTPUT_BITS = 32};
			virtual Uint8  raw8 ();
			virtual Uint16 raw16();
			virtual Uint64 raw64();
			virtual int get_native_output_size() const;
			virtual void add_entropy8 (Uint8 );
			virtual void add_entropy16(Uint16);
			virtual void add_entropy64(Uint64);
		};
		class vRNG64 : public vRNG {
		public:
			enum {OUTPUT_BITS = 64};
			virtual Uint8  raw8 ();
			virtual Uint16 raw16();
			virtual Uint32 raw32();
			virtual int get_native_output_size() const;
			virtual void add_entropy8 (Uint8 );
			virtual void add_entropy16(Uint16);
			virtual void add_entropy32(Uint32);
		};
		namespace OUTPUT_TYPES {enum {
	//		SIMPLE_1 = 0,     //one of 8,16,32,64 as _raw()
			NORMAL_1 = 1,     //one of 8,16,32,64 as raw ## X ()
			NORMAL_ALL = 2,   //all of 8,16,32,64 as raw ## X ()
	//		TEMPLATED_ALL = 3,//all of 8,16,32,64 as raw ## X() AND as _raw<X>()
		};}
//		enum DISTRIBUTIONS_TYPE {
//			DISTRIBUTIONS_TYPE__NONE = 0,
//			DISTRIBUTIONS_TYPE__NORMAL = 1
//		};
//		namespace SEEDING_TYPES { enum {
//			SEEDING_TYPE_INT = 1,//seed(Uint64)
//			SEEDING_TYPE_VRNG = 2//seed(vRNG *)
//		};}
//		enum INTERNAL_STATES_VALID {
//			INTERNAL_STATES_VALID__LOW = 0,//don't manually modify state
//			INTERNAL_STATES_VALID__MED = 1,//manually modification of state not recommended, but unlikely to go horribly wrong
//			INTERNAL_STATES_VALID__HIGH = 2,//may manually modify state; don't expect decent results from low entropy states though
//			INTERNAL_STATES_VALID__ALL = 3//all states pretty much equally valid
//		};
//		enum STATE_TRANSITION_TYPE {
//			UNKNOWN = 0,
//			IRREVERSIBLE_MULTI_CYCLIC = 1,
//			REVERSIBLE_MULTI_CYCLIC = 2,
//			REVERSIBLE_SINGLE_CYCLE = 3,
//			IRREVERSIBLE_SINGLE_CYCLE = 4
//		};
		namespace FLAG { enum {
			SUPPORTS_FASTFORWARD = 1<<0,//also includes rewind
			SUPPORTS_ENTROPY_ACCUMULATION = 1<<1,//supports add_entropy*
			CRYPTOGRAPHIC_SECURITY = 1<<2,
			USES_SPECIFIED = 1<<3,//true if all the other USES_* flags are properly set
			USES_MULTIPLICATION = 1<<4,
			USES_COMPLEX_INSTRUCTIONS = 1<<5,//division, sqrt, exp, log, etc
			USES_VARIABLE_SHIFTS = 1<<6,
			USES_INDIRECTION = 1<<7,
			USES_CYCLIC_BUFFER = 1<<8,
			USES_FLOW_CONTROL = 1<<9,//very simple flow control is not counted
			USES_BIT_SCANS = 1<<10,//bsf & bsr opcodes on x86
			USES_OTHER_WORD_SIZES = 1 << 11,//uses mathematical primitives that do not match the size of its output
			ENDIAN_SAFE = 1<<12,//single flag for output (raw*) and input (add_entropy*)
			OUTPUT_IS_BUFFERED = 1<<13,
			OUTPUT_IS_HASHED = 1<<14,
			STATE_UNAVAILABLE = 1<<15,//don't trust any state-walking operations other than simple seeding (never true on recommended RNGs)
			SEEDING_UNSUPPORTED = 1<<16,//PRNG does not support conventional seeding (example: an RNG that just returns data from standard input)
			NEEDS_GENERIC_SEEDING = 1 << 31,
		};}
		typedef vRNG PolymorphicRNG;
		typedef vRNG8 PolymorphicRNG8;
		typedef vRNG16 PolymorphicRNG16;
		typedef vRNG32 PolymorphicRNG32;
		typedef vRNG64 PolymorphicRNG64;
		namespace Polymorphic {
			using PractRand::RNGs::vRNG;
			using PractRand::RNGs::vRNG8;
			using PractRand::RNGs::vRNG16;
			using PractRand::RNGs::vRNG32;
			using PractRand::RNGs::vRNG64;
			typedef vRNG PolymorphicRNG;
			typedef vRNG8 PolymorphicRNG8;
			typedef vRNG16 PolymorphicRNG16;
			typedef vRNG32 PolymorphicRNG32;
			typedef vRNG64 PolymorphicRNG64;
		}
	}//namespace RNGs
	namespace Tests { union TestBlock; }
}//namespace PractRand
#endif//_practrand_rng_basics_h
