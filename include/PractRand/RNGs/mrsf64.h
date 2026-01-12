namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/mtsf.cpp
			class mrsf64 {// Mrs. Fast 64: Multiply Rotate Small Fast 64-bit PRNG
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 64,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION
				};
			protected:
				Uint64 a, b;
			public:
				Uint64 raw64();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void seed_fast(Uint64 s);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class mrsf64 : public vRNG64 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(mrsf64)
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void seed_fast(Uint64 s);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(mrsf64)
	}
}
