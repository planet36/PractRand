namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/mtc.cpp
			class mrc32 {// Mr. Count: or Multiply Rotate Count 32-bit PRNG
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 32,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION
				};
			protected:
				Uint32 a, b, counter;
			public:
				Uint32 raw32();
				void seed(Uint64 seed_value);
				void seed(Uint64 seed_low, Uint64 seed_high);//only the bottom 64 bits are actually used
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class mrc32 : public vRNG32 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(mrc32)
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(mrc32)
	}
}
