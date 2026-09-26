#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand::RNGs {
		namespace Raw {
			class efiix8x48 {
			public:
				static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1;
				static constexpr int OUTPUT_BITS = 8;
				static constexpr int FLAGS = FLAG::USES_SPECIFIED | FLAG::USES_INDIRECTION | FLAG::USES_CYCLIC_BUFFER | FLAG::ENDIAN_SAFE;
			protected:
				using Word = uint8_t;
				static constexpr int ITERATION_SIZE_L2 = 5;
				static constexpr int ITERATION_SIZE = 1 << ITERATION_SIZE_L2;
				static constexpr int INDIRECTION_SIZE_L2 = 4;
				static constexpr int INDIRECTION_SIZE = 1 << INDIRECTION_SIZE_L2;
				Word indirection_table[INDIRECTION_SIZE], iteration_table[ITERATION_SIZE];
				Word i, a, b, c;
			public:
				~efiix8x48();
				uint8_t raw8();
				void seed(uint64_t s1, uint64_t s2, uint64_t s3, uint64_t s4);
				void seed(uint64_t s) { seed(s, s, s, s); }
				void seed(vRNG* source_rng);
				void walk_state(StateWalkingObject* walker);
				//				static void self_test();
			};
		}

		namespace Polymorphic {
			class efiix8x48 final : public vRNG8 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(efiix8x48)
				void seed(uint64_t s1, uint64_t s2, uint64_t s3, uint64_t s4);
				void seed(uint64_t s) override;
				void seed(vRNG* seeder_rng) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(efiix8x48)
}
