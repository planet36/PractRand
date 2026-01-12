namespace PractRand {
	namespace RNGs {
		namespace Raw {
			//implemented in RNGs/mrc.cpp
			class mrc16 {// Mr. Count: or Multiply Rotate Count 16-bit PRNG
			public:
				enum {
					OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_1,
					OUTPUT_BITS = 16,
					FLAGS = FLAG::ENDIAN_SAFE | FLAG::USES_SPECIFIED | FLAG::USES_MULTIPLICATION
				};
			protected:
				Uint16 a, b, counter;
			public:
				Uint16 raw16();
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
				void walk_state(StateWalkingObject *walker);
			};
		}

		namespace Polymorphic {
			class mrc16 : public vRNG16 {
				PRACTRAND__POLYMORPHIC_RNG_BASICS_H(mrc16)
				void seed(Uint64 seed_low, Uint64 seed_high = 0);
			};
		}
		PRACTRAND__LIGHT_WEIGHT_RNG(mrc16)
	}
}
