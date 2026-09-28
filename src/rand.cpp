#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

//for use in seeding & self-tests:
#include "PractRand/RNGs/all.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <limits>
#include <print>
#include <sstream>
#include <string>
#include <vector>

namespace PractRand {
	const char* version_str = "0.95" "-planet36";
	SEED_AUTO_TYPE SEED_AUTO;
	SEED_NONE_TYPE SEED_NONE;
	void print_err(const char* msg) { std::println(stderr, "{}", msg ? msg : "PractRand: internal error"); }
	void (*error_callback)(const char*) = nullptr;
	void issue_error ( const char* msg) {
		if (error_callback) { error_callback(msg); }
		else {
			if (msg) std::println(stderr, "{}", msg);
			std::exit(1);
		}
	}
	void hook_error_handler(void (*callback)(const char* msg)) {
		error_callback = callback;
	}
	class SerializingStateWalker final : public StateWalkingObject {
	public:
		char* buffer;
		std::size_t max_size;
		std::size_t size_used{};
		SerializingStateWalker( char* buffer_, std::size_t max_size_ )
			: buffer(buffer_), max_size(max_size_)
		{}
		void push(uint8_t value) {
			std::size_t index = size_used++;
			if (index < max_size) buffer[index] = value;
		}
		void handle(bool& v) override {push(v ? 1 : 0);}
		void handle(uint8_t& v) override {push(v);}
		void handle(uint16_t& v) override {auto tmp=uint8_t (v); handle(tmp); tmp=uint8_t (v>> 8); handle(tmp);}
		void handle(uint32_t& v) override {auto tmp=uint16_t(v); handle(tmp); tmp=uint16_t(v>>16); handle(tmp);}
		void handle(uint64_t& v) override {auto tmp=uint32_t(v); handle(tmp); tmp=uint32_t(v>>32); handle(tmp);}
		void handle(float& v) override {
			//uses excess bits to hopefully safely handle floats that might not be exactly IEEE
			bool sign = v < 0;
			v = std::abs(v);
			int exp = 0;
			double n = std::frexp(v, &exp);
			auto tmp_exp = uint16_t( (exp<<1) + (sign?1:0));
			handle(tmp_exp);
			auto tmp_sig = uint32_t(std::ldexp(n - 0.5, 32));
			handle(tmp_sig);
		}
		void handle(double& v) override {
			//uses excess bits to hopefully safely handle floats that might not be exactly IEEE
			bool sign = v < 0;
			v = std::abs(v);
			int exp = 0;
			double n = std::frexp(v, &exp);
			uint32_t tmp_exp = exp + (sign?0x80000000:0);
			handle(tmp_exp);
			auto tmp_sig = uint64_t(std::ldexp(n - 0.5, 64));
			handle(tmp_sig);
		}
		[[nodiscard]] uint32_t get_properties() const override {return FLAG_READ_ONLY;}
	};
	class DeserializingStateWalker : public StateWalkingObject {
	public:
		const char* buffer;
		std::size_t max_size;
		std::size_t size_used{};
		DeserializingStateWalker( const char* buffer_, std::size_t max_size_ )
			: buffer(buffer_), max_size(max_size_)
		{}
		uint8_t pop8() {
			std::size_t index = size_used++;
			if (index < max_size) return buffer[index];
			else return 0;
		}
		uint16_t pop16() {uint16_t tmp=pop8 (); tmp|=uint16_t(pop8 ())<< 8; return tmp;}
		uint32_t pop32() {uint32_t tmp=pop16(); tmp|=uint32_t(pop16())<<16; return tmp;}
		uint64_t pop64() {uint64_t tmp=pop32(); tmp|=uint64_t(pop32())<<32; return tmp;}
		void handle(bool& v) override {v = pop8() != 0;}
		void handle(uint8_t& v) override {v = pop8();}
		void handle(uint16_t& v) override {v = pop16();}
		void handle(uint32_t& v) override {v = pop32();}
		void handle(uint64_t& v) override {v = pop64();}
		void handle(float& v) override {
			uint16_t tmp_exp = pop16();
			uint32_t tmp_sig = pop32();
			bool sign = (tmp_exp & 1) != 0;
			int exp = tmp_exp >> 1;
			if (exp >= 0x4000) exp -= 0x8000;
			v = (sign ? -1.0F : 1.0F) * float(std::ldexp(static_cast<double>(tmp_sig), exp-32));
		}
		void handle(double& v) override {
			uint16_t tmp_exp = pop16();
			uint64_t tmp_sig = pop64();
			bool sign = (tmp_exp & 1) != 0;
			int exp = tmp_exp >> 1;
			if (exp >= 0x4000) exp -= 0x8000;
			v = (sign ? -1.0 : 1.0) * std::ldexp(static_cast<double>(tmp_sig), exp-64);
		}
		[[nodiscard]] uint32_t get_properties() const override {return 0;}
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
		void handle(bool& v) override { pre(); outbuf << v; }
		void handle(uint8_t& v) override { pre(); outbuf << static_cast<unsigned int>(v); }
		void handle(uint16_t& v) override { pre(); outbuf << v; }
		void handle(uint32_t& v) override { pre(); outbuf << v; }
		void handle(uint64_t& v) override { pre(); outbuf << v; }
		void handle(float& v) override { pre(); outbuf << v; }
		void handle(double& v) override { pre(); outbuf << v; }

		uint32_t get_properties() const override { return FLAG_CLUMSY; }
		std::string get_string() const { return outbuf.str(); }
		void reset() { first = true; outbuf.str(""); }
	};
	class GenericIntegerSeedingStateWalker : public StateWalkingObject {
	public:
		PractRand::RNGs::Raw::arbee seeder;
		explicit GenericIntegerSeedingStateWalker(uint64_t seed) : seeder(seed) {}
		void handle(bool& v) override {v = (seeder.raw8() & 1) != 0;}
		void handle(uint8_t& v) override {v = seeder.raw8 ();}
		void handle(uint16_t& v) override {v = seeder.raw16();}
		void handle(uint32_t& v) override {v = seeder.raw32();}
		void handle(uint64_t& v) override {v = seeder.raw64();}
		void handle(float&) override {issue_error("RNGs with default integer seeding should not contain floating point values");}
		void handle(double&) override {issue_error("RNGs with default integer seeding should not contain floating point values");}
		[[nodiscard]] uint32_t get_properties() const override { return FLAG_CLUMSY | FLAG_SEEDER; }
	};
	class GenericSeedingStateWalker : public StateWalkingObject {
	public:
		PractRand::RNGs::vRNG* seeder;
		explicit GenericSeedingStateWalker(RNGs::vRNG* seeder_) : seeder(seeder_) {}
		void handle(bool& v) override { v = (seeder->raw8() & 1) != 0; }
		void handle(uint8_t& v) override {v = seeder->raw8 ();}
		void handle(uint16_t& v) override {v = seeder->raw16();}
		void handle(uint32_t& v) override {v = seeder->raw32();}
		void handle(uint64_t& v) override {v = seeder->raw64();}
		void handle(float&) override {issue_error("RNGs with default seeding should not contain floating point values");}
		void handle(double&) override {issue_error("RNGs with default seeding should not contain floating point values");}
		[[nodiscard]] uint32_t get_properties() const override {return FLAG_CLUMSY | FLAG_SEEDER;}
	};
	namespace AutoSeeder {
		constexpr int POOL_SIZE = 5;
		static bool initialized = false;
		static bool enough_entropy_found;
		static uint64_t shared_entropy[POOL_SIZE] = {0};
		static void initialize() {// NOT thread-safe
			if (initialized) return;
			initialized = true;
			PractRand::RNGs::Polymorphic::sha2_based_pool entropy_pool;
			//PractRand::RNGs::Polymorphic::arbee entropy_pool;
			enough_entropy_found = entropy_pool.add_entropy_automatically(1);
			for (auto& i : shared_entropy) i = entropy_pool.raw64();
		}
		static void get_autoseed_fixed_entropy(uint64_t entropy[5], [[maybe_unused]] const void* target) {
			// NOT thread-safe the first time it's run
			if (!initialized) initialize();
			static thread_local uint64_t intrathread_count = 0;
			static thread_local uint64_t thread_identifier = 0;
			if (!intrathread_count) thread_identifier = PractRand::Internals::issue_unique_identifier();
			++intrathread_count;
			entropy[0] = shared_entropy[0] ^ intrathread_count;
			entropy[1] = shared_entropy[1] ^ thread_identifier;
			entropy[2] = shared_entropy[2];
			entropy[3] = shared_entropy[3];
			entropy[4] = shared_entropy[4];
#if 0
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
			uint64_t target_addr = (uint64_t)target;
			//bool which_time = (end_clock == start_clock);
			//uint64_t timestamp = which_time ? end_time : end_clock;
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
			explicit AutoSeedingStateWalker([[maybe_unused]] const void* target) {
				//get_autoseed_entropy(&seeder, target);
				uint32_t seed_and_iv[10] = {0};
				get_autoseed_fixed_entropy(reinterpret_cast<uint64_t*>(&seed_and_iv[0]), &seeder);

				//I would prefer ChaCha, but the license is not 100% clear atm
				//PractRand::RNGs::Polymorphic::chacha bootstrap(seed_and_iv, false);
				PractRand::RNGs::Polymorphic::salsa bootstrap(seed_and_iv, false);

				std::memset(seed_and_iv, 0, sizeof(seed_and_iv));
				seeder.seed(bootstrap.raw64(), bootstrap.raw64(), bootstrap.raw64(), bootstrap.raw64());
			}
			void handle(bool& v) override {v = (seeder.raw8() & 1) != 0;}
			void handle(uint8_t& v) override {v = seeder.raw8 ();}
			void handle(uint16_t& v) override {v = seeder.raw16();}
			void handle(uint32_t& v) override {v = seeder.raw32();}
			void handle(uint64_t& v) override {v = seeder.raw64();}
			void handle([[maybe_unused]] float& v) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			void handle([[maybe_unused]] double& v) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			[[nodiscard]] uint32_t get_properties() const override {return FLAG_CLUMSY | FLAG_SEEDER;}
		};
		class CryptoAutoSeedingStateWalker : public StateWalkingObject {
		public:
			//PractRand::RNGs::Polymorphic::sha2_based_pool seeder;
			PractRand::RNGs::Polymorphic::trivium seeder;
			explicit CryptoAutoSeedingStateWalker(void* ptr1) : seeder(PractRand::SEED_NONE) {
				PractRand::RNGs::Polymorphic::sha2_based_pool entropy_pool;
				if (!entropy_pool.add_entropy_automatically(0))
					issue_error("PractRand: failed to obtain entropy for cryptographic quality autoseeding");
				uint64_t extra[5];
				get_autoseed_fixed_entropy(extra, ptr1);
				for (const auto i : extra) entropy_pool.add_entropy64(i);
				std::memset(extra, 0, sizeof(extra));
				if constexpr (false) {
					seeder.seed(&entropy_pool);//probably stronger, but... not 100% sure with Trivium
				}
				else {
					//what we're supposed to do:
					constexpr int B = 20;
					uint8_t s[B];
					for (auto& i : s) i = entropy_pool.raw8();
					seeder.seed(s, B);
					std::memset(s, 0, B);
					for (int i = 0; i < 4; i++) seeder.raw64();//strength of Trivium might be improved by skipping a few outputs after seeding
				}
			}
			void handle(bool& v) override {v = (seeder.raw8() & 1) != 0;}
			void handle(uint8_t& v) override {v = seeder.raw8 ();}
			void handle(uint16_t& v) override {v = seeder.raw16();}
			void handle(uint32_t& v) override {v = seeder.raw32();}
			void handle(uint64_t& v) override {v = seeder.raw64();}
			void handle(float&) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			void handle(double&) override {issue_error("RNGs with auto-seeding should not contain floating point values");}
			[[nodiscard]] uint32_t get_properties() const override {return FLAG_CLUMSY | FLAG_SEEDER;}
		};
	}
	uint32_t randi_fast_implementation(uint32_t random_value, uint32_t max) {
		return uint32_t((uint64_t(max) * random_value) >> 32);
	}
	StateWalkingObject* int_to_rng_seed(uint64_t i) {
		return new GenericIntegerSeedingStateWalker(i);
	}
	StateWalkingObject* vrng_to_rng_seeder(RNGs::vRNG* rng) {
		return new GenericSeedingStateWalker(rng);
	}
	StateWalkingObject* get_autoseeder(const void* target) {
		return new AutoSeeder::AutoSeedingStateWalker(target);
	}
	namespace Internals {
		void test_random_access(PractRand::RNGs::vRNG* rng, PractRand::RNGs::vRNG* known_good, uint64_t period_low64, uint64_t period_high64) {
			uint64_t seed = known_good->raw64();
			uint8_t a1 = 0, a2 = 0, a3 = 0, b1 = 0, b2 = 0, b3 = 0;
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
				int64_t how_far_ = known_good->raw64();
				while (how_far_ == std::numeric_limits<decltype(how_far_)>::min()) how_far_ = known_good->raw64();//we can't negate this value, so the code would fail
				if (how_far_ > 0) rng->seek_forward(how_far_);
				else if (how_far_ < 0) rng->seek_backward(-how_far_);
				a1 = rng->raw8(); a2 = rng->raw8(); a3 = rng->raw8();
				int64_t delta = how_far_;
				rng->seed(seed);
				while (delta) {
					if (delta > 0) {
						uint64_t adjust = known_good->randli(delta + 1);
						rng->seek_forward(adjust);
						delta -= adjust;
					}
					else {
						uint64_t adjust = known_good->randli(1 - delta);
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
		long vRNG::serialize( char* buffer, long buffer_size ) {//returns serialized size, or zero on failure
			SerializingStateWalker serializer(buffer, buffer_size);
			walk_state(&serializer);
			if (serializer.size_used <= static_cast<std::size_t>(buffer_size)) return serializer.size_used;
			return 0;
		}
		char* vRNG::serialize( std::size_t* size_ ) {//returns malloced block, or NULL on error, sets *size to size of block
			SerializingStateWalker byte_counter(nullptr, 0);
			walk_state(&byte_counter);
			std::size_t size = byte_counter.size_used;
			*size_ = size;
			if (!size) return nullptr;
			char* buffer = static_cast<char*>(std::malloc(size));
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
		bool vRNG::deserialize( const char* buffer, size_t size ) {//returns number of bytes used, or zero on error
			DeserializingStateWalker deserializer(buffer, size);
			walk_state(&deserializer);
			return deserializer.size_used == size;
		}
		void vRNG::seed(uint64_t seed) {
			GenericIntegerSeedingStateWalker walker(seed);
			walk_state(&walker);
			flush_buffers();
		}
		void vRNG::seed_fast(uint64_t seed) {
			GenericIntegerSeedingStateWalker walker(seed);
			walk_state(&walker);
			flush_buffers();
		}
		void vRNG::seed(vRNG* rng) {
			GenericSeedingStateWalker walker(rng);
			walk_state(&walker);
			flush_buffers();
		}
		void vRNG::autoseed() {
			AutoSeeder::AutoSeedingStateWalker walker(this);
			walk_state(&walker);
			flush_buffers();
		}
		uint32_t vRNG::randi(uint32_t max) {
			max -= 1;
			uint32_t mask = max;
			mask |= mask >> 1; mask |= mask >>  2; mask |= mask >> 4;
			mask |= mask >> 8; mask |= mask >> 16;
			uint32_t tmp = 0;
			do {
				tmp = raw32() & mask;
			} while (tmp > max);
			return tmp;
		}
		uint64_t vRNG::randli(uint64_t max) {
			max -= 1;
			uint64_t mask = max;
			mask |= mask >> 1; mask |= mask >>  2; mask |= mask >>  4;
			mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;
			uint64_t tmp = 0;
			do {
				tmp = raw64() & mask;
			} while (tmp > max);
			return tmp;
		}
		uint32_t vRNG::randi_fast(uint32_t max) {
			return randi_fast_implementation(raw32(), max);
		}
		float vRNG::randf() {return static_cast<float>(raw32() & ((static_cast<uint32_t>(1) << 24)-1)) * static_cast<float>(1.0/16777216.0);}
		double vRNG::randlf() {return static_cast<double>(raw64() & ((static_cast<uint64_t>(1) << 53)-1)) * (1.0/9007199254740992.0);}
		double vRNG::gaussian() { return Internals::generate_gaussian_fast(raw64()); }
		uint64_t vRNG::get_flags() const {return 0;}
		void vRNG::seek_forward128 (uint64_t, uint64_t) {}
		void vRNG::seek_backward128(uint64_t, uint64_t) {}
		void vRNG::flush_buffers() {}


		int vRNG8::get_native_output_size() const {return 8;}
		uint16_t vRNG8::raw16() {
			uint16_t r = raw8();
			return r | (uint16_t(raw8())<<8);
		}
		uint32_t vRNG8::raw32() {
			uint32_t r = raw8();
			r = r | (uint32_t(raw8()) << 8);
			r = r | (uint32_t(raw8()) << 16);
			return r | (uint32_t(raw8()) << 24);
		}
		uint64_t vRNG8::raw64() {
			uint64_t r = raw8();
			r = r | (uint64_t(raw8()) << 8);
			r = r | (uint64_t(raw8()) << 16);
			r = r | (uint64_t(raw8()) << 24);
			r = r | (uint64_t(raw8()) << 32);
			r = r | (uint64_t(raw8()) << 40);
			r = r | (uint64_t(raw8()) << 48);
			return r | (uint64_t(raw8()) << 56);
		}

		int vRNG16::get_native_output_size() const {return 16;}
		uint8_t  vRNG16::raw8()  {return uint8_t(raw16());}
		uint32_t vRNG16::raw32() {
			uint32_t r = raw16();
			return r | (uint32_t(raw16()) << 16);
		}
		uint64_t vRNG16::raw64() {
			uint64_t r = raw16();
			r = r | (uint64_t(raw16()) << 16);
			r = r | (uint64_t(raw16()) << 32);
			return r | (uint64_t(raw16()) << 48);
		}

		int vRNG32::get_native_output_size() const {return 32;}
		uint8_t  vRNG32::raw8()  {return uint8_t (raw32());}
		uint16_t vRNG32::raw16() {return uint16_t(raw32());}
		uint64_t vRNG32::raw64() {
			uint32_t r = raw32();
			return r | (uint64_t(raw32()) << 32);
		}

		int vRNG64::get_native_output_size() const {return 64;}
		uint8_t  vRNG64::raw8()  {return uint8_t (raw64());}
		uint16_t vRNG64::raw16() {return uint16_t(raw64());}
		uint32_t vRNG64::raw32() {return uint32_t(raw64());}

		void vRNG::reset_entropy()       {}
		void vRNG::add_entropy8 (uint8_t)  {}
		void vRNG::add_entropy16(uint16_t) {}
		void vRNG::add_entropy32(uint32_t) {}
		void vRNG::add_entropy64(uint64_t) {}
		void vRNG::add_entropy_N(const void* _data, std::size_t length) {
			const auto* data = static_cast<const uint8_t*>(_data);
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

		{RNGs::Polymorphic::chacha rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 1ULL << 36); }
		{RNGs::Polymorphic::salsa rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 1ULL << 36); }
		{RNGs::Polymorphic::xsm32 rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 1); }
		{RNGs::Polymorphic::xsm64 rng(PractRand::SEED_NONE); PractRand::Internals::test_random_access(&rng, &known_good, 0, 0); }
	}
	bool initialize_PractRand() {
		if (!AutoSeeder::initialized)
			AutoSeeder::initialize();
#if 0
		union {
			uint64_t as64[1]{};
			uint32_t as32[2];
			uint16_t as16[4];
			uint8_t   as8[8];
		};
		as64[0] = 0x0123456789ABCDEFULL;
		if (as8[7] != 0x01) {
			issue_error("PractRand - endianness configured incorrectly");
		}
#endif
		return AutoSeeder::enough_entropy_found;
	}
}
