#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			class arbee {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;
				static constexpr int OUTPUT_BITS = 64;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::ENDIAN_SAFE | FLAG::SUPPORTS_ENTROPY_ACCUMULATION;
			protected:
				uint64_t a{}, b{}, c{}, d{}, i{};
				void mix();
			public:
				arbee() {reset_entropy();}
				explicit arbee(uint64_t s) {seed(s);}
				arbee(uint64_t s1, uint64_t s2, uint64_t s3, uint64_t s4) {seed(s1,s2,s3,s4);}
				explicit arbee(SEED_NONE_TYPE ) {}
				explicit arbee(SEED_AUTO_TYPE ) {StateWalkingObject* walker = get_autoseeder(this); this->walk_state(walker); delete walker;}
				uint8_t  raw8 () {return uint8_t (raw64());}
				uint16_t raw16() {return uint16_t(raw64());}
				uint32_t raw32() {return uint32_t(raw64());}
				uint64_t raw64();
				void seed(uint64_t s);
				void seed(uint64_t seed1, uint64_t seed2, uint64_t seed3, uint64_t seed4);//custom seeding
				void seed(vRNG* rng);
				void walk_state(StateWalkingObject* walker);
				void reset_entropy();
				void add_entropy8 (uint8_t  value);
				void add_entropy16(uint16_t value);
				void add_entropy32(uint32_t value);
				void add_entropy64(uint64_t value);
				void add_entropy_N(const void*, size_t length);
				void flush_buffers() {mix();}
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class arbee final : public vRNG64 {
			public:
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::ENDIAN_SAFE | FLAG::SUPPORTS_ENTROPY_ACCUMULATION;
				Raw::arbee implementation;
				[[nodiscard]] uint64_t get_flags() const override;
				[[nodiscard]] std::string get_name() const override;
				explicit arbee(uint64_t s) : implementation(s) {}
				arbee(uint64_t s1, uint64_t s2, uint64_t s3, uint64_t s4) : implementation(s1,s2,s3,s4) {}
				explicit arbee(vRNG* seeder) {seed(seeder);}
				explicit arbee(SEED_AUTO_TYPE ) {autoseed();}
				explicit arbee(SEED_NONE_TYPE ) {}
				arbee() = default;
				uint8_t  raw8 () override;
				uint16_t raw16() override;
				uint32_t raw32() override;
				uint64_t raw64() override;
				void seed(uint64_t s) override;
				void seed(uint64_t s1, uint64_t s2, uint64_t s3, uint64_t s4);
				void seed(vRNG* rng) override;
				void walk_state(StateWalkingObject* walker) override;
				void reset_entropy() override;
				void add_entropy8 (uint8_t  value) override;
				void add_entropy16(uint16_t value) override;
				void add_entropy32(uint32_t value) override;
				void add_entropy64(uint64_t value) override;
				void add_entropy_N(const void*, size_t length) override;
				void flush_buffers() override;
			};
		}
		namespace LightWeight {
			using Raw::arbee;
		};
}
