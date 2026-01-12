namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/mtc.cpp
			class mrc64 {// Mr. Count: or Multiply Rotate Count 64-bit PRNG
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 64,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION
				};
			protected:
				Uint64 a, b, counter;
			public:
				Uint64 raw64();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class mrc64 : public vRNG64 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(mrc64)
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(mrc64)
	}
}
