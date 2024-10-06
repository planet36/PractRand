#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand {
	namespace RNGs {
		namespace Raw {
			class efiix8x384 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 8;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::USES_INDIRECTION | FLAG::USES_CYCLIC_BUFFER | FLAG::ENDIAN_SAFE;
			protected:
				typedef Uint8 Word;
				static constexpr int ITERATION_SIZE_L2 = 7;
				static constexpr int ITERATION_SIZE = 1 << ITERATION_SIZE_L2;
				static constexpr int INDIRECTION_SIZE_L2 = 8;
				static constexpr int INDIRECTION_SIZE = 1 << INDIRECTION_SIZE_L2;
				Word indirection_table[INDIRECTION_SIZE], iteration_table[ITERATION_SIZE];
				Word i, a, b, c;
			public:
				~efiix8x384();
				Uint8 raw8();
				void seed(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4);
				void seed(Uint64 s) {seed(s,s,s,s);}
				void seed(vRNG *source_rng);
				void walk_state(StateWalkingObject *walker);
//				static void self_test();
			};
		}

		namespace Polymorphic {
			class efiix8x384 : public vRNG8 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(efiix8x384)
				void seed(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4);
				void seed(Uint64 s);
				void seed(vRNG *source_rng);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(efiix8x384)
	}
}
