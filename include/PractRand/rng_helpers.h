#pragma once

#include <cstdint>

namespace PractRand {
	class StateWalkingObject {
	public:
		virtual ~StateWalkingObject() = default;
		virtual void handle(bool&) = 0;
		/*virtual void handle(unsigned char &) = 0;
		virtual void handle(unsigned short &) = 0;
		virtual void handle(unsigned int &) = 0;
		virtual void handle(unsigned long &) = 0;
		virtual void handle(unsigned long long &) = 0;*/
		virtual void handle(uint8_t&) = 0;
		virtual void handle(uint16_t&) = 0;
		virtual void handle(uint32_t&) = 0;
		virtual void handle(uint64_t&) = 0;

		virtual void handle(float&) = 0;
		virtual void handle(double&) = 0;

		//purposes:
		// 1. seeding
		// 2. serialization (serialize, deserialize, measure state size)
		// 3. avalanche style testing (measure state size / props, tweak bits in state, etc)
		// 4. ???
		static constexpr int FLAG_READ_ONLY = 1; // does not make changes
		static constexpr int FLAG_WRITE_ONLY = 2;// result does not depend upon prior state
		static constexpr int FLAG_CLUMSY = 4;    // may violate invariants (if also FLAG_READ_ONLY then only wants to see state visible to clumsy writers)
		static constexpr int FLAG_SEEDER = 8;    // some RNGs may have extra invariants enforced only on seeded states (minimum distance away from other seeded states on cycle)
		[[nodiscard]] virtual uint32_t get_properties() const = 0;
		[[nodiscard]] bool is_read_only() const { return (get_properties() & FLAG_READ_ONLY) != 0; }
		[[nodiscard]] bool is_write_only() const { return (get_properties() & FLAG_WRITE_ONLY) != 0; }
		[[nodiscard]] bool is_clumsy() const { return (get_properties() & FLAG_CLUMSY) != 0; }
		[[nodiscard]] bool is_seeder() const {return (get_properties() & FLAG_SEEDER) != 0;}

		/*void handle(signed char      &v) {handle((unsigned char)v);}
		void handle(signed short     &v) {handle((unsigned short)v);}
		void handle(signed int       &v) {handle((unsigned int)v);}
		void handle(signed long      &v) {handle((unsigned long)v);}
		void handle(signed long long &v) {handle((unsigned long long)v);}*/
		void handle(int8_t& v) {handle(reinterpret_cast<uint8_t&>(v));}
		void handle(int16_t&v) {handle(reinterpret_cast<uint16_t&>(v));}
		void handle(int32_t&v) {handle(reinterpret_cast<uint32_t&>(v));}
		void handle(int64_t&v) {handle(reinterpret_cast<uint64_t&>(v));}

		StateWalkingObject& operator<<(uint8_t& v) {handle(v);return *this;}
		StateWalkingObject& operator<<(uint16_t&v) {handle(v);return *this;}
		StateWalkingObject& operator<<(uint32_t&v) {handle(v);return *this;}
		StateWalkingObject& operator<<(uint64_t&v) {handle(v);return *this;}
		StateWalkingObject& operator<<(int8_t& v) {handle(v);return *this;}
		StateWalkingObject& operator<<(int16_t&v) {handle(v);return *this;}
		StateWalkingObject& operator<<(int32_t&v) {handle(v);return *this;}
		StateWalkingObject& operator<<(int64_t&v) {handle(v);return *this;}
		StateWalkingObject& operator<<(float& v){handle(v);return *this;}
		StateWalkingObject& operator<<(double& v){handle(v);return *this;}
	};
	namespace RNGs {
		class vRNG;
	}
	uint32_t randi_fast_implementation(uint32_t random_value, uint32_t max);
	StateWalkingObject* int_to_rng_seeder(uint64_t);//must be deleted after use
	StateWalkingObject* vrng_to_rng_seeder(RNGs::vRNG*);//must be deleted after use
	//An autoseeding constructor passes its own unseeded object, which GCC's LTO then flags with
	//-Wmaybe-uninitialized.  The pointer is never read, and access(none) tells GCC so.
#if __has_cpp_attribute(gnu::access)
	[[gnu::access(none, 1)]]
#endif
	StateWalkingObject* get_autoseeder(const void*);//must be deleted after use
}
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PRACTRAND_POLYMORPHIC_RNG_BASICS_H(RNG) public:\
		static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;\
		static constexpr int OUTPUT_BITS = Raw:: RNG ::OUTPUT_BITS;\
		static constexpr int FLAGS = Raw:: RNG ::FLAGS;\
		Raw:: RNG implementation{};\
		explicit RNG (uint64_t s) {seed(s);}\
		explicit RNG (vRNG* seeder) {seed(seeder);}\
		explicit RNG (SEED_AUTO_TYPE ) {autoseed();}\
		explicit RNG (SEED_NONE_TYPE ) {}\
		uint8_t  raw8 () override;\
		uint16_t raw16() override;\
		uint32_t raw32() override;\
		uint64_t raw64() override;\
		using vRNG::seed;\
		uint64_t get_flags() const override;\
		std::string get_name() const override;\
		void walk_state(StateWalkingObject* walker) override;

#if defined PRACTRAND_NO_LIGHT_WEIGHT_RNGS
#define PRACTRAND_LIGHT_WEIGHT_RNG(RNG)
#define PRACTRAND_LIGHT_WEIGHT_ENTROPY_POOL(RNG)
#else // ! PRACTRAND_NO_LIGHT_WEIGHT_RNGS
#include "rng_adaptors.h"
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PRACTRAND_LIGHT_WEIGHT_RNG(RNG) 	\
	namespace LightWeight {\
		typedef PractRand::RNGs::Adaptors::RAW_TO_LIGHT_WEIGHT_RNG<PractRand::RNGs::Raw:: RNG > RNG;\
	}
#endif//PRACTRAND_PROVIDE_LIGHT_WEIGHT_RNGS
