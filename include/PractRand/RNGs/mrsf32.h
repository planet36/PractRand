namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/mtsf.cpp
			class mrsf32 {// Mrs. Fast 32: Multiply Rotate Small Fast 32-bit PRNG
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 32,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION
				};
			protected:
				Uint32 a, b;
			public:
				Uint32 raw32();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void seed_fast(Uint64 s);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class mrsf32 : public vRNG32 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(mrsf32)
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void seed_fast(Uint64 s);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(mrsf32)

	}
}
