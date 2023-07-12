
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/jsf.cpp
			class jsf32 {
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 32,
					FLAGS = FLAG::USES_SPECIFIED | FLAG::ENDIAN_SAFE
				};
			protected:
				Uint32 a, b, c, d;
			public:
				Uint32 raw32();
				void seed(Uint64 s);
				void seed_fast(Uint64 s);
				void seed(vRNG *seeder_rng);
				void seed(Uint32 seed1, Uint32 seed2, Uint32 seed3, Uint32 seed4);//custom seeding
				void walk_state(StateWalkingObject *walker);
				//static void self_test();
			};
		}
		
		namespace Polymorphic {
			class jsf32 final : public vRNG32 {
				PRACTRAND_POLYMORPHIC_RNG_BASICS_H(jsf32)
				void seed(Uint64 s) override;
				void seed_fast(Uint64 s) override;
			};
		}
		PRACTRAND_LIGHT_WEIGHT_RNG(jsf32)
	}
}
