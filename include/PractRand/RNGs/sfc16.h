namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/sfc.cpp
			class sfc16 {
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 16,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED
				};
			protected:
				Uint16 a, b, c, counter;
			public:
				Uint16 raw16();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void seed_fast(Uint64 s);
				void walk_state(StateWalkingObject *walker);
			};
		}
		
		namespace Polymorphic {
			class sfc16 : public vRNG16 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(sfc16)
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void seed_fast(Uint64 s);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(sfc16)
	}
}
