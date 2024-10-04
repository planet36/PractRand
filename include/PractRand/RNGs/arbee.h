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
				Uint64 a, b, c, d, i;
				void mix();
			public:
				arbee() {reset_entropy();}
				arbee(Uint64 s) {seed(s);}
				arbee(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4) {seed(s1,s2,s3,s4);}
				arbee(SEED_NONE_TYPE ) {}
				arbee(SEED_AUTO_TYPE ) {StateWalkingObject *walker = get_autoseeder(this); this->walk_state(walker); delete walker;}
				Uint8  raw8 () {return Uint8 (raw64());}
				Uint16 raw16() {return Uint16(raw64());}
				Uint32 raw32() {return Uint32(raw64());}
				Uint64 raw64();
				void seed(Uint64 s);
				void seed(Uint64 seed1, Uint64 seed2, Uint64 seed3, Uint64 seed4);//custom seeding
				void seed(vRNG *rng);
				void walk_state(StateWalkingObject *walker);
				void reset_entropy();
				void add_entropy8 (Uint8  value);
				void add_entropy16(Uint16 value);
				void add_entropy32(Uint32 value);
				void add_entropy64(Uint64 value);
				void add_entropy_N(const void *, size_t length);
				void flush_buffers() {mix();}
				//static void self_test();
			};
		}

		namespace Polymorphic {
			class arbee final : public vRNG64 {
			public:
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::ENDIAN_SAFE | FLAG::SUPPORTS_ENTROPY_ACCUMULATION;
				Raw::arbee implementation;
				Uint64 get_flags() const override;
				std::string get_name() const override;
				arbee(Uint64 s) : implementation(s) {}
				arbee(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4) : implementation(s1,s2,s3,s4) {}
				arbee(vRNG *seeder) {seed(seeder);}
				arbee(SEED_AUTO_TYPE ) {autoseed();}
				arbee(SEED_NONE_TYPE ) {}
				arbee() = default;
				Uint8  raw8 () override;
				Uint16 raw16() override;
				Uint32 raw32() override;
				Uint64 raw64() override;
				void seed(Uint64 s) override;
				void seed(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4);
				void seed(vRNG *rng) override;
				void walk_state(StateWalkingObject *walker) override;
				void reset_entropy() override;
				void add_entropy8 (Uint8  value) override;
				void add_entropy16(Uint16 value) override;
				void add_entropy32(Uint32 value) override;
				void add_entropy64(Uint64 value) override;
				void add_entropy_N(const void *, size_t length) override;
				void flush_buffers() override;
			};
		}
		namespace LightWeight {
			using Raw::arbee;
		};
}
