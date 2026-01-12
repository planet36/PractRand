
/*
RNGs in the other directory (and the NotRecommended namespace) are not intended for real world use,
only for research; as such they may get pretty sloppy in some areas

This particular subset of them is RNGs that don't fit in to the other categories.  
Most often this means use of multiplication in non-LCG-like ways, but it could mean a variety of other things too.  
*/

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				//
				class garthy16 : public vRNG16 {
					Uint16 value, scale, counter;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class garthy32 : public vRNG32 {
					Uint32 value, scale, counter;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//both sides of the multiply are pseudo-random values in this RNG
				class binarymult16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class binarymult32 : public vRNG32 {
					Uint32 a, b, c, d;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class mmr16 : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mmr32 : public vRNG32 {
					Uint32 a, b, c;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//uses multiplication, rightshifts, xors, that kind of stuff
				class rxmult16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//these are similar to my mwlac algorithm, but lower quality
				class multish2x64 : public vRNG64 {
					Uint64 a, b;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class multish3x32 : public vRNG32 {
					Uint32 a, b, c;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class multish4x16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class mwrca16 : public vRNG16 {
					Uint16 a, b;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrca32 : public vRNG32 {
					Uint32 a, b;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrca64 : public vRNG64 {
					Uint64 a, b, counter;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrcc16 : public vRNG16 {
					Uint16 a, b, counter;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrcc32 : public vRNG32 {
					Uint32 a, b, counter;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrcca16 : public vRNG16 {
					Uint16 a, b, counter;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrcca32 : public vRNG32 {
					Uint32 a, b, counter;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwrcca64 : public vRNG64 {
					Uint64 a, b, counter;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//the 16 bit variant of the old version of my mwlac algorithm
				class old_mwlac16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varA : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varB : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varC : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varD : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varE : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varF : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


				class mwc64x : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class cxm64_varqual : public vRNG64 {
					Uint64 low, high;
					int num_mult;
				public:
					cxm64_varqual(int num_mult_) : num_mult(num_mult_) {}
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//see http://www.drdobbs.com/tools/229625477
				class mo_Cmfr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Cmr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Cmr32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class murmlac32 : public vRNG32 {
					Uint32 state1, state2;
					int rounds;
				public:
					Uint32 raw32();
					murmlac32(int rounds_) : rounds(rounds_) {}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//multiplication (by a counter), rotate
				class mulcr64 : public vRNG64 {
					Uint64 a, b, count;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcr32 : public vRNG32 {
					Uint32 a, b, count;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcr16 : public vRNG16 {
					Uint32 a, b, count;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class mulcrax8 : public vRNG8 {
					Uint8 a, b, count_low, count_high;
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcrax16 : public vRNG16 {
					Uint16 a, b, count_low, count_high;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcrax32 : public vRNG32 {
					Uint32 a, b, count_low, count_high;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcrx16 : public vRNG16 {
					Uint16 a, count_low, count_high;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcrx32 : public vRNG32 {
					Uint32 a, count_low, count_high;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcrx64 : public vRNG64 {
					Uint64 a, count_low, count_high;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class varrotA : public vRNG16 {
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class varrotB : public vRNG16 {
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class varrotC : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class varrotD : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class varrotE : public vRNG16 {
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class varrotF : public vRNG16 {
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lxwm16 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lxwm32 : public vRNG32 {
					Uint32 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mrsf16 : public vRNG16 {
					Uint16 a, b;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					void seed(Uint64 s);
				};
				class hierarchyA : public vRNG32 {
					Uint32 a, b, c;
					Uint64 top_a, top_b, top_c;
					void set_frame();
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
			}
		}
	}
}