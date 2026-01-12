namespace PractRand {
	namespace RNGs {
		namespace Raw {
			class efiix8x48 {
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 8,
					FLAGS = FLAG::USES_SPECIFIED | FLAG::USES_INDIRECTION | FLAG::USES_CYCLIC_BUFFER | FLAG::ENDIAN_SAFE
				};
			protected:
				typedef Uint8 Word;
				enum {
					ITERATION_SIZE_L2 = 5,
					ITERATION_SIZE = 1 << ITERATION_SIZE_L2,
					INDIRECTION_SIZE_L2 = 4,
					INDIRECTION_SIZE = 1 << INDIRECTION_SIZE_L2
				};
				Word indirection_table[INDIRECTION_SIZE], iteration_table[ITERATION_SIZE];
				Word i, a, b, c;
			public:
				~efiix8x48();
				Uint8 raw8();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				//void seed(vRNG *source_rng);
				void walk_state(StateWalkingObject *walker);

				void reset_entropy();
				void add_entropy8(Uint8  value);
				//void add_entropy16(Uint16 value);
				//void add_entropy32(Uint32 value);
				//void add_entropy64(Uint64 value);
				void add_entropy_N(const void *, size_t length);
				void flush_buffers();
				//				static void self_test();
			};
		}

		namespace Polymorphic {
			class efiix8x48 : public vRNG8 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(efiix8x48)
				//void seed(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4);
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				//void seed(vRNG *source_rng);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(efiix8x48)
	}
}
