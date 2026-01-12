namespace PractRand {
	namespace RNGs {
		namespace Raw {
			namespace NotRecommended {
				class mt19937 {
				public:
					enum {
						OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
						OUTPUT_BITS = 32,
						FLAGS = FLAG::OUTPUT_IS_BUFFERED | FLAG::OUTPUT_IS_HASHED | FLAG::ENDIAN_SAFE
					};
				protected:
					enum { ARRAY_SIZE = 624, OFFSET = 397 };
					Uint32 state[ARRAY_SIZE];
					Uint32 used;
					void _advance_state();
				public:
					void flush_buffers() { used = ARRAY_SIZE; }
					Uint32 raw32();
					void walk_state(StateWalkingObject *walker);

					//seeds < 2**32 use the standard MT19937 seeding algorithm
					//seeds >= 2**32 use a nonstandard MT19937 seeding algorithm
					void seed(Uint64 seed_low, Uint64 seed_high=0);
					void seed(Uint32 s[], int seed_length);//alternate seeding algorithm added to MT19937 in 2002

					Uint32 untempered_raw32() {
						if (used >= ARRAY_SIZE) {
							_advance_state();
							return state[used++];
						}
						else return state[used++];
					}
					static void self_test();
				};
			}
		}
		
		namespace Polymorphic {
			namespace NotRecommended {
				class mt19937 : public vRNG32 {
					Raw::NotRecommended::mt19937 implementation;
				public:
					enum { OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL, OUTPUT_BITS = Raw::NotRecommended::mt19937::OUTPUT_BITS, FLAGS = Raw::NotRecommended::mt19937::FLAGS };
					mt19937(Uint64 seed_low, Uint64 seed_high=0) { seed(seed_low, seed_high); }
					mt19937(vRNG *seeder) { seed(seeder); }
					mt19937(SEED_AUTO_TYPE) { autoseed(); }
					mt19937(SEED_NONE_TYPE) {}
					mt19937() {}
					mt19937(Uint32 s[], int seed_length) { seed(s, seed_length); }
					Uint32 raw32();
					using vRNG::seed;
					Uint64 get_flags() const;
					std::string get_name() const;
					void walk_state(StateWalkingObject *walker);
					void seed(Uint64 seed_low, Uint64 seed_high=0);
					void seed(Uint32 s[], int seed_length);//alternate seeding algorithm added to MT in 2002
					void flush_buffers();
				};
			}
		}
		namespace LightWeight {
			namespace NotRecommended {
				typedef PractRand::RNGs::Adaptors::RAW_TO_LIGHT_WEIGHT_RNG<PractRand::RNGs::Raw::NotRecommended::mt19937 > mt19937;
			}
		}
	}
}
