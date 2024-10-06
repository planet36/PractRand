#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand {
	namespace RNGs {
		namespace Raw {
			class efiix32x48 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 32;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::USES_INDIRECTION | FLAG::USES_CYCLIC_BUFFER | FLAG::ENDIAN_SAFE;
			protected:
				typedef Uint32 Word;
				static constexpr int ITERATION_SIZE = 32;
				static constexpr int INDIRECTION_SIZE = 16;
				Word iteration_table[ITERATION_SIZE], indirection_table[INDIRECTION_SIZE];
				Word i, a, b, c;
			public:
				~efiix32x48();
				Uint32 raw32();
				//void seed(const Word *seeds, int num_seeds, int seeding_quality=4);
				void seed(Uint64 s) { seed(s, s, s, s); }
				void seed(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4);
				void seed(vRNG *source_rng);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class efiix32x48 : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(efiix32x48)
				//void seed(const Uint32 *seeds, int num_seeds, int seeding_quality=4) {implementation.seed(seeds, num_seeds, seeding_quality);}
				void seed(Uint64 s);
				void seed(Uint64 s1, Uint64 s2, Uint64 s3, Uint64 s4);
				void seed(vRNG *source_rng);
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(efiix32x48)
	}
}
