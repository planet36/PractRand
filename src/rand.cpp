#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

//for use in seeding & self-tests:
#include "PractRand/RNGs/all.h"

namespace PractRand {
	const char *version_str = "0.95" "-planet36";
	SEED_AUTO_TYPE SEED_AUTO;
	SEED_NONE_TYPE SEED_NONE;
	void print_err(const char* msg) { (void)std::fputs(msg, stderr); (void)std::fputc('\n', stderr); }
	void (*error_callback)(const char *) = nullptr;
	void issue_error ( const char *msg) {
		if (error_callback) error_callback(msg);
		else {
			if (msg) (void)std::fprintf(stderr, "%s\n", msg);
			std::exit(1);
		}
	}
	void hook_error_handler(void (*callback)(const char *msg)) {
		error_callback = callback;
	}
	class SerializingStateWalker final : public StateWalkingObject {
	public:
		char *buffer;
		std::size_t max_size;
		std::size_t size_used{};
		SerializingStateWalker( char *buffer_, std::size_t max_size_ )
			: buffer(buffer_), max_size(max_size_)
		{}
		void push(Uint8 value) {
			std::size_t index = size_used++;
			if (index < max_size) buffer[index] = value;
		}
		void handle(bool   &v) override {push(v ? 1 : 0);}
		void handle(Uint8  &v) override {push(v);}
		void handle(Uint16 &v) override {auto tmp=Uint8 (v); handle(tmp); tmp=Uint8 (v>> 8); handle(tmp);}
		void handle(Uint32 &v) override {auto tmp=Uint16(v); handle(tmp); tmp=Uint16(v>>16); handle(tmp);}
		void handle(Uint64 &v) override {auto tmp=Uint32(v); handle(tmp); tmp=Uint32(v>>32); handle(tmp);}
		void handle(float  &v) override {
			//uses excess bits to hopefully safely handle floats that might not be exactly IEEE
			bool sign = v < 0;
			v = std::abs(v);
			int exp;
			double n = std::frexp(v, &exp);
			auto tmp_exp = Uint16( (exp<<1) + (sign?1:0));
			handle(tmp_exp);
			auto tmp_sig = Uint32(std::ldexp(n - 0.5, 32));
			handle(tmp_sig);
		}
		void handle(double &v) override {
			//uses excess bits to hopefully safely handle floats that might not be exactly IEEE
			bool sign = v < 0;
			v = std::abs(v);
			int exp;
			double n = std::frexp(v, &exp);
			Uint32 tmp_exp = exp + (sign?0x80000000:0);
			handle(tmp_exp);
			auto tmp_sig = Uint64(std::ldexp(n - 0.5, 64));
			handle(tmp_sig);
		}
		Uint32 get_properties() const override {return FLAG_READ_ONLY;}
	};
	class DeserializingStateWalker : public StateWalkingObject {
	public:
		const char *buffer;
		std::size_t max_size;
		std::size_t size_used{};
		DeserializingStateWalker( const char *buffer_, std::size_t max_size_ )
			: buffer(buffer_), max_size(max_size_)
		{}
		Uint8 pop8() {
			std::size_t index = size_used++;
			if (index < max_size) return buffer[index];
			else return 0;
		}
		Uint16 pop16() {Uint16 tmp=pop8 (); tmp|=Uint16(pop8 ())<< 8; return tmp;}
		Uint32 pop32() {Uint32 tmp=pop16(); tmp|=Uint32(pop16())<<16; return tmp;}
		Uint64 pop64() {Uint64 tmp=pop32(); tmp|=Uint64(pop32())<<32; return tmp;}
		void handle(bool   &v) override {v = pop8() ? true : false;}
		void handle(Uint8  &v) override {v = pop8();}
		void handle(Uint16 &v) override {v = pop16();}
		void handle(Uint32 &v) override {v = pop32();}
		void handle(Uint64 &v) override {v = pop64();}
		void handle(float  &v) override {
			Uint16 tmp_exp = pop16();
			Uint32 tmp_sig = pop32();
			bool sign = (tmp_exp & 1) ? true : false;
			int exp = tmp_exp >> 1;
			if (exp >= 0x4000) exp -= 0x8000;
			v = (sign ? -1.0f : 1.0f) * float(std::ldexp(static_cast<double>(tmp_sig), exp-32));
		}
		void handle(double &v) override {
			Uint16 tmp_exp = pop16();
			Uint64 tmp_sig = pop64();
			bool sign = (tmp_exp & 1) ? true : false;
			int exp = tmp_exp >> 1;
			if (exp >= 0x4000) exp -= 0x8000;
			v = (sign ? -1.0 : 1.0) * std::ldexp(static_cast<double>(tmp_sig), exp-64);
		}
		Uint32 get_properties() const override {return 0;}
	};
	class PrintingStateWalker : public StateWalkingObject {
		std::ostringstream outbuf;
		bool first{true};
		void pre() {
			if (first) first = false;
			else outbuf << ",";
		}
	public:
		PrintingStateWalker() = default;
		void handle(bool   &v) override { pre(); outbuf << v; }
		void handle(Uint8  &v) override { pre(); outbuf << v; }
		void handle(Uint16 &v) override { pre(); outbuf << v; }
		void handle(Uint32 &v) override { pre(); outbuf << v; }
		void handle(Uint64 &v) override { pre(); outbuf << v; }
		void handle(float  &v) override { pre(); outbuf << v; }
		void handle(double &v) override { pre(); outbuf << v; }

		Uint32 get_properties() const override { return FLAG_CLUMSY; }
		std::string get_string() const { return outbuf.str(); }
		void reset() { first = true; outbuf.str(""); }
	};
	class GenericIntegerSeedingStateWalker : public StateWalkingObject {
	public:
		PractRand::RNGs::Raw::arbee seeder;
		GenericIntegerSeedingStateWalker(Uint64 seed) : seeder(seed) {}
		void handle(bool   &v) override {v = (seeder.raw8() & 1) ? true : false;}
		void handle(Uint8  &v) override {v = seeder.raw8 ();}
		void handle(Uint16 &v) override {v = seeder.raw16();}
		void handle(Uint32 &v) override {v = seeder.raw32();}
		void handle(Uint64 &v) override {v = seeder.raw64();}
		void handle(float  &) override {issue_error("RNGs with default integer seeding should not contain floating point values");}
		void handle(double &) override {issue_error("RNGs with default integer seeding should not contain floating point values");}
		Uint32 get_properties() const override { return FLAG_CLUMSY | FLAG_SEEDER; }
	};
	class GenericSeedingStateWalker : public StateWalkingObject {
	public:
		PractRand::RNGs::vRNG *seeder;
		GenericSeedingStateWalker(RNGs::vRNG *seeder_) : seeder(seeder_) {}
		void handle(bool   &v) override { v = (seeder->raw8() & 1) ? true : false; }
		void handle(Uint8  &v) override {v = seeder->raw8 ();}
		void handle(Uint16 &v) override {v = seeder->raw16();}
		void handle(Uint32 &v) override {v = seeder->raw32();}
		void handle(Uint64 &v) override {v = seeder->raw64();}
		void handle(float  &) override {issue_error("RNGs with default seeding should not contain floating point values");}
		void handle(double &) override {issue_error("RNGs with default seeding should not contain floating point values");}
		Uint32 get_properties() const override {return FLAG_CLUMSY | FLAG_SEEDER;}
	};
	namespace AutoSeeder {
		constexpr int POOL_SIZE = 5;
		static bool initialized = false;
		static bool enough_entropy_found;
		static Uint64 shared_entropy[POOL_SIZE] = {0};
		static void initialize() {// NOT thread-safe
			if (initialized) return;
			initialized = true;
			PractRand::RNGs::Polymorphic::sha2_based_pool entropy_pool;
			//PractRand::RNGs::Polymorphic::arbee entropy_pool;
			enough_entropy_found = entropy_pool.add_entropy_automatically(1);
			for (auto & i : shared_entropy) i = entropy_pool.raw64();
		}
		static void get_autoseed_fixed_entropy(Uint64 entropy[5], [[maybe_unused]] const void *target) {
			// NOT thread-safe the first time it's run
			if (!initialized) initialize();
#if defined PRACTRAND_THREAD_LOCAL_STORAGE
			static PRACTRAND_THREAD_LOCAL_STORAGE Uint64 intrathread_count = 0;
			static PRACTRAND_THREAD_LOCAL_STORAGE Uint64 thread_identifier = 0;
			if (!intrathread_count) thread_identifier = PractRand::Internals::issue_unique_identifier();
			++intrathread_count;
			entropy[0] = shared_entropy[0] ^ intrathread_count;
			entropy[1] = shared_entropy[1] ^ thread_identifier;
			entropy[2] = shared_entropy[2];
			entropy[3] = shared_entropy[3];
			entropy[4] = shared_entropy[4];
#else
			//I think the combination of an address and a non-zero time span guarantees uniqueness on a per-run basis
			//but it requires blocking on the time span
			std::clock_t start_clock = std::clock();
			std::time_t start_time = std::time(nullptr);
			std::clock_t end_clock;
			std::time_t end_time;
			while (true) {
				end_clock = std::clock();
				end_time = std::time(nullptr);
				if (start_clock != end_clock) break;
				if (start_time  != end_time) break;
			}
			Uint64 target_addr = (Uint64)target;
			//bool which_time = (end_clock == start_clock);
			//Uint64 timestamp = which_time ? end_time : end_clock;
			entropy[0] = shared_entropy[0] ^ target_addr;
			entropy[1] = shared_entropy[1] ^ end_time;
			entropy[2] = shared_entropy[2] ^ end_clock;
			entropy[3] = shared_entropy[3];
			entropy[4] = shared_entropy[4];
#endif
		}
		class AutoSeedingStateWalker : public StateWalkingObject {
		public:
			PractRand::RNGs::Polymorphic::arbee seeder;
			AutoSeedingStateWalker([[maybe_unused]] const void *target) {
				//get_autoseed_entropy(&seeder, target);
				Uint32 seed_and_iv[10] = {0};
				get_autoseed_fixed_entropy(reinterpret_cast<Uint64*>(&seed_and_iv[0]), &seeder);

				//I would prefer ChaCha, but the license is not 100% clear atm
				//PractRand::RNGs::Polymorphic::chacha bootstrap(seed_and_iv, false);
				PractRand::RNGs::Polymorphic::salsa bootstrap(seed_and_iv, false);

				std::memset(seed_and_iv, 0, sizeof(seed_and_iv));
				seeder.seed(bootstrap.raw64(), bootstrap.raw64(), bootstrap.raw64(), bootstrap.raw64());
			}
			void handle(bool   &v) override {v = (seeder.raw8() & 1) ? true : false;}
			void handle(Uint8  &v) override {v = seeder.raw8 ();}
			void handle(Uint16 &v) override {v = seeder.raw16();}
			void handle(Uint32 &v) override {v = seeder.raw32();}
			void handle(Uint64 &v) override {v = seeder.raw64();}
			void handle([[maybe_unused]] float  &v) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			void handle([[maybe_unused]] double &v) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			Uint32 get_properties() const override {return FLAG_CLUMSY | FLAG_SEEDER;}
		};
		class CryptoAutoSeedingStateWalker : public StateWalkingObject {
		public:
			//PractRand::RNGs::Polymorphic::sha2_based_pool seeder;
			PractRand::RNGs::Polymorphic::trivium seeder;
			CryptoAutoSeedingStateWalker(void *ptr1) : seeder(PractRand::SEED_NONE) {
				PractRand::RNGs::Polymorphic::sha2_based_pool entropy_pool;
				if (!entropy_pool.add_entropy_automatically())
					issue_error("PractRand: failed to obtain entropy for cryptographic quality autoseeding");
				Uint64 extra[5];
				get_autoseed_fixed_entropy(extra, ptr1);
				for (const auto i : extra) entropy_pool.add_entropy64(i);
				std::memset(extra, 0, sizeof(extra));
				if constexpr (false) {
					seeder.seed(&entropy_pool);//probably stronger, but... not 100% sure with Trivium
				}
				else {
					//what we're supposed to do:
					constexpr int B = 20;
					Uint8 s[B];
					for (auto & i : s) i = entropy_pool.raw8();
					seeder.seed(s, B);
					std::memset(s, 0, B);
					for (int i = 0; i < 4; i++) seeder.raw64();//strength of Trivium might be improved by skipping a few outputs after seeding
				}
			}
			void handle(bool   &v) override {v = (seeder.raw8() & 1) ? true : false;}
			void handle(Uint8  &v) override {v = seeder.raw8 ();}
			void handle(Uint16 &v) override {v = seeder.raw16();}
			void handle(Uint32 &v) override {v = seeder.raw32();}
			void handle(Uint64 &v) override {v = seeder.raw64();}
			void handle(float  &) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			void handle(double &) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			Uint32 get_properties() const override {return FLAG_CLUMSY | FLAG_SEEDER;}
		};
	}
	Uint32 randi_fast_implementation(Uint32 random_value, Uint32 max) {
		return Uint32((Uint64(max) * random_value) >> 32);
	}
	StateWalkingObject *int_to_rng_seed(Uint64 i) {
		return new GenericIntegerSeedingStateWalker(i);
	}
	StateWalkingObject *vrng_to_rng_seeder(RNGs::vRNG *rng) {
		return new GenericSeedingStateWalker(rng);
	}
	StateWalkingObject *get_autoseeder(const void *target) {
		return new AutoSeeder::AutoSeedingStateWalker(target);
	}
	namespace Internals {
		void test_random_access(PractRand::RNGs::vRNG *rng, PractRand::RNGs::vRNG *known_good, Uint64 period_low64, Uint64 period_high64) {
			Uint64 seed = known_good->raw64();
			Uint8 a1, a2, a3, b1, b2, b3;
			//basic check
			rng->seed(seed);
			//a1 = rng->raw8(); a2 = rng->raw8(); a3 = rng->raw8();
			(void)rng->raw8(); (void)rng->raw8(); (void)rng->raw8();
			a1 = rng->raw8(); a2 = rng->raw8(); a3 = rng->raw8();
			rng->seed(seed);
			rng->seek_forward(3);
			b1 = rng->raw8(); b2 = rng->raw8(); b3 = rng->raw8();
			if (a1 != b1 || a2 != b2 || a3 != b3) PractRand::issue_error("PractRand::test_random_access failed (1)");
			//check a longer range seek
			seed = known_good->raw64();
			rng->seed(seed);
			int how_far = known_good->randi(13179);
			for (int i = 0; i < how_far; i++) rng->raw8();
			a1 = rng->raw8(); a2 = rng->raw8(); a3 = rng->raw8();
			rng->seed(seed);
			rng->seek_forward(how_far);
			b1 = rng->raw8(); b2 = rng->raw8(); b3 = rng->raw8();
			if (a1 != b1 || a2 != b2 || a3 != b3) PractRand::issue_error("PractRand::test_random_access failed (2)");
			//check a more exotic pattern of seeks, with some longer still seeks
			for (int i = 0; i < 10; i++) {
				seed = known_good->raw64();
				rng->seed(seed);
				Sint64 how_far_ = known_good->raw64();
				while (how_far_ == std::numeric_limits<decltype(how_far_)>::min()) how_far_ = known_good->raw64();//we can't negate this value, so the code would fail
				if (how_far_ > 0) rng->seek_forward(how_far_);
				else if (how_far_ < 0) rng->seek_backward(-how_far_);
				a1 = rng->raw8(); a2 = rng->raw8(); a3 = rng->raw8();
				Sint64 delta = how_far_;
				rng->seed(seed);
				while (delta) {
					if (delta > 0) {
						Uint64 adjust = known_good->randli(delta + 1);
						rng->seek_forward(adjust);
						delta -= adjust;
					}
					else {
						Uint64 adjust = known_good->randli(1 - delta);
						rng->seek_backward(adjust);
						delta += adjust;
					}
				}
				b1 = rng->raw8(); b2 = rng->raw8(); b3 = rng->raw8();
				if (a1 != b1 || a2 != b2 || a3 != b3) PractRand::issue_error("PractRand::test_random_access failed (3)");
			}
			//check cycle length if one was reported
			//could add a check on prime factorization of cycle lengths, but that would be more trouble than it's worth right now
			if (period_low64 || period_high64) {
				rng->seed(seed);
				a1 = rng->raw8(); a2 = rng->raw8(); a3 = rng->raw8();
				rng->seed(seed);
				rng->seek_forward128(period_low64, period_high64);
				b1 = rng->raw8(); b2 = rng->raw8(); b3 = rng->raw8();
				if (a1 != b1 || a2 != b2 || a3 != b3) PractRand::issue_error("PractRand::test_random_access failed (4)");
				rng->seed(seed);
				rng->seek_backward128(period_low64, period_high64);
				b1 = rng->raw8(); b2 = rng->raw8(); b3 = rng->raw8();
				if (a1 != b1 || a2 != b2 || a3 != b3) PractRand::issue_error("PractRand::test_random_access failed (5)");
			}
		}
	}
	namespace RNGs {
		vRNG::~vRNG() = default;
		long vRNG::serialize( char *buffer, long buffer_size ) {//returns serialized size, or zero on failure
			SerializingStateWalker serializer(buffer, buffer_size);
			walk_state(&serializer);
			if (serializer.size_used <= static_cast<std::size_t>(buffer_size)) return serializer.size_used;
			return 0;
		}
		char *vRNG::serialize( std::size_t *size_ ) {//returns malloced block, or NULL on error, sets *size to size of block
			SerializingStateWalker byte_counter(nullptr, 0);
			walk_state(&byte_counter);
			std::size_t size = byte_counter.size_used;
			*size_ = size;
			if (!size) return nullptr;
			char *buffer = static_cast<char*>(std::malloc(size));
			SerializingStateWalker serializer(buffer, size);
			if (serializer.size_used != size) {
				std::free(buffer);
				return nullptr;
			}
			return buffer;
		}
		std::string vRNG::print_state() {
			PrintingStateWalker printer;
			walk_state(&printer);
			return printer.get_string();
		}
		bool vRNG::deserialize( const char *buffer, size_t size ) {//returns number of bytes used, or zero on error
			DeserializingStateWalker deserializer(buffer, size);
			walk_state(&deserializer);
			return deserializer.size_used == size;
		}
		void vRNG::seed(Uint64 i) {
			GenericIntegerSeedingStateWalker walker(i);
			walk_state(&walker);
			flush_buffers();
		}
		void vRNG::seed_fast(Uint64 i) {
			GenericIntegerSeedingStateWalker walker(i);
			walk_state(&walker);
			flush_buffers();
		}
		void vRNG::seed(vRNG *rng) {
			GenericSeedingStateWalker walker(rng);
			walk_state(&walker);
			flush_buffers();
		}
		void vRNG::autoseed() {
			AutoSeeder::AutoSeedingStateWalker walker(this);
			walk_state(&walker);
			flush_buffers();
		}
		Uint32 vRNG::randi(Uint32 max) {
			PRACTRAND_RANDI_IMPLEMENTATION(max)
		}
		Uint64 vRNG::randli(Uint64 max) {
			PRACTRAND_RANDLI_IMPLEMENTATION(max)
		}
		Uint32 vRNG::randi_fast(Uint32 max) {
			return randi_fast_implementation(raw32(), max);
		}
		float vRNG::randf() {PRACTRAND_RANDF_IMPLEMENTATION(*this)}
		double vRNG::randlf() {PRACTRAND_RANDLF_IMPLEMENTATION(*this)}
		double vRNG::gaussian() { return Internals::generate_gaussian_fast(raw64()); }
		Uint64 vRNG::get_flags() const {return 0;}
		void vRNG::seek_forward128 (Uint64, Uint64) {}
		void vRNG::seek_backward128(Uint64, Uint64) {}
		void vRNG::flush_buffers() {}


		int vRNG8::get_native_output_size() const {return 8;}
		Uint16 vRNG8::raw16() {
			Uint16 r = raw8();
			return r | (Uint16(raw8())<<8);
		}
		Uint32 vRNG8::raw32() {
			Uint32 r = raw8();
			r = r | (Uint32(raw8()) << 8);
			r = r | (Uint32(raw8()) << 16);
			return r | (Uint32(raw8()) << 24);
		}
		Uint64 vRNG8::raw64() {
			Uint64 r = raw8();
			r = r | (Uint64(raw8()) << 8);
			r = r | (Uint64(raw8()) << 16);
			r = r | (Uint64(raw8()) << 24);
			r = r | (Uint64(raw8()) << 32);
			r = r | (Uint64(raw8()) << 40);
			r = r | (Uint64(raw8()) << 48);
			return r | (Uint64(raw8()) << 56);
		}

		int vRNG16::get_native_output_size() const {return 16;}
		Uint8  vRNG16::raw8()  {return Uint8(raw16());}
		Uint32 vRNG16::raw32() {
			Uint32 r = raw16();
			return r | (Uint32(raw16()) << 16);
		}
		Uint64 vRNG16::raw64() {
			Uint64 r = raw16();
			r = r | (Uint64(raw16()) << 16);
			r = r | (Uint64(raw16()) << 32);
			return r | (Uint64(raw16()) << 48);
		}

		int vRNG32::get_native_output_size() const {return 32;}
		Uint8  vRNG32::raw8()  {return Uint8 (raw32());}
		Uint16 vRNG32::raw16() {return Uint16(raw32());}
		Uint64 vRNG32::raw64() {
			Uint32 r = raw32();
			return r | (Uint64(raw32()) << 32);
		}

		int vRNG64::get_native_output_size() const {return 64;}
		Uint8  vRNG64::raw8()  {return Uint8 (raw64());}
		Uint16 vRNG64::raw16() {return Uint16(raw64());}
		Uint32 vRNG64::raw32() {return Uint32(raw64());}

		void vRNG::reset_entropy()       {}
		void vRNG::add_entropy8 (Uint8)  {}
		void vRNG::add_entropy16(Uint16) {}
		void vRNG::add_entropy32(Uint32) {}
		void vRNG::add_entropy64(Uint64) {}
		void vRNG::add_entropy_N(const void *_data, std::size_t length) {
			auto *data = static_cast<const Uint8*>(_data);
			for (unsigned long i = 0; i < length; i++) add_entropy8(data[i]);
		}
		bool vRNG::add_entropy_automatically(int milliseconds) {
			return PractRand::Internals::add_entropy_automatically(this, milliseconds);
		}
	}
	void self_test_PractRand() {
		RNGs::Raw::mt19937::self_test();
		RNGs::Raw::hc256::self_test();
		RNGs::Raw::isaac32x256::self_test();
		RNGs::Raw::trivium::self_test();
		RNGs::Raw::chacha::self_test();
		RNGs::Raw::salsa::self_test();

		RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);

		{RNGs::Polymorphic::chacha rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 1ull << 36); }
		{RNGs::Polymorphic::salsa rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 1ull << 36); }
		{RNGs::Polymorphic::xsm32 rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 1); }
		{RNGs::Polymorphic::xsm64 rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 0); }
	}
	bool initialize_PractRand() {
		if (!AutoSeeder::initialized)
			AutoSeeder::initialize();
		union {
			Uint64 as64[1];
			Uint32 as32[2];
			Uint16 as16[4];
			Uint8   as8[8];
		};
		as64[0] = 0x0123456789ABCDEFull;
#if defined PRACTRAND_TARGET_IS_LITTLE_ENDIAN
		if (as8[7] != 0x01) {
#elif defined PRACTRAND_TARGET_IS_BIG_ENDIAN
		if (as8[0] != 0x01) {
#endif
			issue_error("PractRand - endianness configured incorrectly");
		}
		return AutoSeeder::enough_entropy_found;
	}
}
