#pragma once

#include <cstdint>
#include <string>

namespace PractRand {
	extern const char* version_str;//like "0.91", for PractRand 0.91

	bool initialize_PractRand(); //returns true normally
	//will return false if it failed to find a good source of entropy
	//  in which case the autoseeding mechanism may have trouble
	//  Usually not catastrophic, but some programs might want to abort if that happens.
	//NOTE: initialize_PractRand() is NOT threadsafe, it should be called before threads get spun off

	void self_test_PractRand();
	void print_err(const char* msg);
	void issue_error(const char* msg = nullptr);//PractRand calls this any time there is an internal error
	void hook_error_handler(void(*callback)(const char*));//this can be used to replace the default behavior of issue_error

	class StateWalkingObject;

	class SEED_AUTO_TYPE {};
	class SEED_NONE_TYPE {};
	extern SEED_AUTO_TYPE SEED_AUTO;
	extern SEED_NONE_TYPE SEED_NONE;

	namespace RNGs {
		class vRNG {
		public:
		//constructors, destructors, seeding, serialization, & low level state manipulation:
			//vRNG(uint64_t seed_) {seed_64(seed_);}
			//vRNG(vRNG *rng_) {seed(rng_);}
			//vRNG(_dummy_SeedingTypeAuto *) {autoseed();}
			//vRNG(_dummy_SeedingTypeNone *) {}
			vRNG() = default;
			virtual ~vRNG();
			vRNG(const vRNG&) = delete;
			vRNG& operator=(const vRNG&) = delete;
			vRNG(vRNG&&) = delete;
			vRNG& operator=(vRNG&&) = delete;
			virtual void seed(uint64_t seed);
			virtual void seed_fast(uint64_t seed);
			virtual void seed(vRNG* rng);
			virtual void autoseed();
			long serialize( char* buffer, long buffer_size );//returns serialized size, or zero on failure
			char* serialize( size_t* size );//returns malloced block, or NULL on error, sets *size to size of block
			bool deserialize( const char* buffer, size_t size );//returns true on success, false on failure
			std::string print_state();//returns RNG state as a comma-delimited sequence of numbers
			virtual void walk_state(StateWalkingObject*) = 0;


		//raw random bits
			virtual uint8_t  raw8 () = 0;
			virtual uint16_t raw16() = 0;
			virtual uint32_t raw32() = 0;
			virtual uint64_t raw64() = 0;
			//virtual void raw_N(uint8_t *, size_t length) = 0;

		//uniform distributions
			uint32_t randi(uint32_t max);
			uint32_t randi(uint32_t min, uint32_t max) {return randi(max-min)+min;}
			uint32_t randi_fast(uint32_t max);
			uint32_t randi_fast(uint32_t min, uint32_t max) {return randi_fast(max-min)+min;}
			uint64_t randli(uint64_t max);
			uint64_t randli(uint64_t min, uint64_t max) {return randli(max-min)+min;}
			float randf();
			float randf(float max) {return randf() * max;}
			float randf(float min, float max) {return randf() * (max-min) + min;}
			double randlf();
			double randlf(double max) {return randlf() * max;}
			double randlf(double min, double max) {return randlf() * (max-min) + min;}

		//non-uniform distributions
			double gaussian();//mean 0.0, stddev 1.0
			double gaussian(double mean, double stddev) { return gaussian() * stddev + mean; }

		//metadata functions
			[[nodiscard]] virtual uint64_t get_flags() const;
			[[nodiscard]] virtual std::string get_name() const = 0;
			[[nodiscard]] virtual int get_native_output_size() const = 0;//generally 8, 16, 32, 64, or -1 (unknown)

		//exotic methods (not supported by many implementations - check flags to see if they support it):
		//exotic methods 1: random access
			virtual void seek_forward128 (uint64_t how_far_low64, uint64_t how_far_high64);
			virtual void seek_backward128(uint64_t how_far_low64, uint64_t how_far_high64);
			void seek_forward (uint64_t how_far) {seek_forward128 (how_far, 0);}
			void seek_backward(uint64_t how_far) {seek_backward128(how_far, 0);}

		//exotic methods 2: entropy pooling
			virtual void reset_entropy();//returns an entropy pool to its default state
			virtual void add_entropy8 (uint8_t );
			virtual void add_entropy16(uint16_t);
			virtual void add_entropy32(uint32_t);
			virtual void add_entropy64(uint64_t);
			//note that "add_entropy_N(&byte_buffer[0], 13)" will typically NOT produce the same state transition
			//  as "add_entropy_N(&byte_buffer[0], 7);add_entropy_N(&byte_buffer[0], 6);"
			virtual void add_entropy_N(const void*, size_t length);

			//add_entropy_automatically returns true if a good amount (>= 128 bits) of entropy was added
			//the milliseconds parameter is the maximum amount of time it is allowed to block while waiting for entropy
			virtual bool add_entropy_automatically(int milliseconds);

			virtual void flush_buffers();// some entropy pooling PRNGs have internal buffers that need to be flushed before inputs can effect outputs - this flushes both input and output buffers

		// C++2011 compatibility:
#if defined PRACTRAND_BOOST_COMPATIBILITY
			typedef uint64_t result_type;
			result_type operator()() {return raw64();}
			static constexpr bool has_fixed_value = true;
			static constexpr result_type min_value = 0;
			static constexpr result_type max_value = ~(result_type)0;
			result_type min() const {return min_value;}
			result_type max() const {return max_value;}
#endif
		};
		class vRNG8 : public vRNG {
		public:
			static constexpr int OUTPUT_BITS = 8;
			uint16_t raw16() override;
			uint32_t raw32() override;
			uint64_t raw64() override;
			[[nodiscard]] int get_native_output_size() const override;
		};
		class vRNG16 : public vRNG {
		public:
			static constexpr int OUTPUT_BITS = 16;
			uint8_t  raw8 () override;
			uint32_t raw32() override;
			uint64_t raw64() override;
			[[nodiscard]] int get_native_output_size() const override;
		};
		class vRNG32 : public vRNG {
		public:
			static constexpr int OUTPUT_BITS = 32;
			uint8_t  raw8 () override;
			uint16_t raw16() override;
			uint64_t raw64() override;
			[[nodiscard]] int get_native_output_size() const override;
		};
		class vRNG64 : public vRNG {
		public:
			static constexpr int OUTPUT_BITS = 64;
			uint8_t  raw8 () override;
			uint16_t raw16() override;
			uint32_t raw32() override;
			[[nodiscard]] int get_native_output_size() const override;
		};
		namespace OUTPUT_TYPES {
		//constexpr int SIMPLE_1 = 0;     //one of 8,16,32,64 as _raw()
		constexpr int NORMAL_1 = 1;       //one of 8,16,32,64 as raw ## X ()
		constexpr int NORMAL_ALL = 2;     //all of 8,16,32,64 as raw ## X ()
		//constexpr int TEMPLATED_ALL = 3;//all of 8,16,32,64 as raw ## X() AND as _raw<X>()
		}
//		enum DISTRIBUTIONS_TYPE {
//			DISTRIBUTIONS_TYPE__NONE = 0,
//			DISTRIBUTIONS_TYPE__NORMAL = 1
//		};
//		namespace SEEDING_TYPES { enum {
//			SEEDING_TYPE_INT = 1,//seed(uint64_t)
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
		namespace FLAG {
		constexpr int SUPPORTS_FASTFORWARD = 1<<0;//also includes rewind
		constexpr int SUPPORTS_ENTROPY_ACCUMULATION = 1<<1;//supports add_entropy*
		constexpr int CRYPTOGRAPHIC_SECURITY = 1<<2;
		constexpr int USES_SPECIFIED = 1<<3;//true if all the other USES_* flags are properly set
		constexpr int USES_MULTIPLICATION = 1<<4;
		constexpr int USES_COMPLEX_INSTRUCTIONS = 1<<5;//division, sqrt, exp, log, etc
		constexpr int USES_VARIABLE_SHIFTS = 1<<6;
		constexpr int USES_INDIRECTION = 1<<7;
		constexpr int USES_CYCLIC_BUFFER = 1<<8;
		constexpr int USES_FLOW_CONTROL = 1<<9;//very simple flow control is not counted
		constexpr int USES_BIT_SCANS = 1<<10;//bsf & bsr opcodes on x86
		constexpr int USES_OTHER_WORD_SIZES = 1 << 11;//uses mathematical primitives that do not match the size of its output
		constexpr int ENDIAN_SAFE = 1<<12;//single flag for output (raw*) and input (add_entropy*)
		constexpr int OUTPUT_IS_BUFFERED = 1<<13;
		constexpr int OUTPUT_IS_HASHED = 1<<14;
		constexpr int STATE_UNAVAILABLE = 1<<15;//don't trust any state-walking operations other than simple seeding (never true on recommended RNGs)
		constexpr int SEEDING_UNSUPPORTED = 1<<16;//PRNG does not support conventional seeding (example: an RNG that just returns data from standard input)
		constexpr int NEEDS_GENERIC_SEEDING = 1<<31;
		}
		using PolymorphicRNG = vRNG;
		using PolymorphicRNG8 = vRNG8;
		using PolymorphicRNG16 = vRNG16;
		using PolymorphicRNG32 = vRNG32;
		using PolymorphicRNG64 = vRNG64;
		namespace Polymorphic {
			using PractRand::RNGs::vRNG;
			using PractRand::RNGs::vRNG8;
			using PractRand::RNGs::vRNG16;
			using PractRand::RNGs::vRNG32;
			using PractRand::RNGs::vRNG64;
			using PolymorphicRNG = vRNG;
			using PolymorphicRNG8 = vRNG8;
			using PolymorphicRNG16 = vRNG16;
			using PolymorphicRNG32 = vRNG32;
			using PolymorphicRNG64 = vRNG64;
		}
	}//namespace RNGs
	namespace Tests { union TestBlock; }
}//namespace PractRand
