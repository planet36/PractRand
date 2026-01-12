namespace PractRand {
	namespace RNGs {
		namespace Raw {
			class efiix8min {
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 8,
					FLAGS = FLAG::USES_SPECIFIED | FLAG::USES_INDIRECTION | FLAG::USES_CYCLIC_BUFFER | FLAG::ENDIAN_SAFE
				};
			protected:
				typedef Uint8 Word;
				enum {
					ITERATION_SIZE_L2 = 0,
					ITERATION_SIZE = 1 << ITERATION_SIZE_L2,
					INDIRECTION_SIZE_L2 = 0,
					INDIRECTION_SIZE = 1 << INDIRECTION_SIZE_L2
				};
				//Word indirection_table[INDIRECTION_SIZE], iteration_table[ITERATION_SIZE];
				Word indirection_table, iteration_table;//with array size reduced to 1, no need to declare them as arrays anymore?
				Word i, a, b, c;
			public:
				Uint8 raw8();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void walk_state(StateWalkingObject *walker);
				//				static void self_test();
			};
		}

		namespace Polymorphic {
			class efiix8min : public vRNG8 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(efiix8min)
				void seed(Uint64 s1, Uint64 s2 = 0);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(efiix8min)
	}
}
