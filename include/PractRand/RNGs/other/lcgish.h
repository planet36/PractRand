
/*
RNGs in the other directory (and the NotRecommended namespace) are not intended for real world use,
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that:
1. use multiplication
2. don't use much indirection, flow control, variable shifts, etc
3. have only a few words of state
4. are likely to have easily detectable bias
*/

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				//similar to the classic LCGs, but with a longer period
				class lcg16of32_extended : public vRNG16 {
					Uint32 state, add;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg32_extended : public vRNG32 {
					Uint32 state, add;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//simple classic LCGs
				class lcg32of64_varqual : public vRNG32 {
					Uint64 state;
					int outshift;
				public:
					lcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg16of64_varqual : public vRNG16 {
					Uint64 state;
					int outshift;
				public:
					lcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg8of64_varqual : public vRNG8 {
					Uint64 state;
					int outshift;
				public:
					lcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg32of128_varqual : public vRNG32 {
					Uint64 low, high;
					int outshift;
				public:
					lcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg16of128_varqual : public vRNG16 {
					Uint64 low, high;
					int outshift;
				public:
					lcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg8of128_varqual : public vRNG8 {
					Uint64 low, high;
					int outshift;
				public:
					lcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//some helpers:
				namespace NP2LCG_Helpers {
					template<int shift, int post_shift_multiplier, int multiplier, typename StateType, typename IntermediateType> class NearP2LCG {
					public:
						StateType state;
						NearP2LCG() : state(1) {}
						void advance() {
							if (sizeof(StateType) == sizeof(IntermediateType)) {
								state *= multiplier;
								state = (state >> shift) * post_shift_multiplier + (state & ((StateType(1) << shift) - 1));
								state = (state >> shift) * post_shift_multiplier + (state & ((StateType(1) << shift) - 1));
							}
							else {
								IntermediateType tmp = state;
								tmp *= multiplier;
								tmp = (tmp >> shift) * post_shift_multiplier + (tmp & ((IntermediateType(1) << shift) - 1));
								tmp = (tmp >> shift) * post_shift_multiplier + (tmp & ((IntermediateType(1) << shift) - 1));
								state = tmp;
							}
						}
						void sanitize() {
							while (state >= (IntermediateType(1) << shift)) {
								IntermediateType tmp = state;
								tmp = (tmp >> shift) * post_shift_multiplier + (tmp & ((IntermediateType(1) << shift) - 1));
								state = tmp;
							}
							if (state < post_shift_multiplier) state += (IntermediateType(1) << shift) - post_shift_multiplier;
						}
						Uint64 seed(Uint64 value) {
							const IntermediateType one_less_than_modulus = (IntermediateType(1) << shift) - post_shift_multiplier - 1;
							state = 1 + (value % one_less_than_modulus);
							sanitize();
							return value / one_less_than_modulus;
						}
					};
					typedef NearP2LCG<10, 3, 177, Uint32, Uint32> np2lcg_10;
					typedef NearP2LCG<11, 9, 913, Uint32, Uint32> np2lcg_11;
					typedef NearP2LCG<13, 1, 2185, Uint32, Uint32> np2lcg_13;
					typedef NearP2LCG<17, 1, 8485, Uint32, Uint32> np2lcg_17;
					typedef NearP2LCG<19, 1, 5601, Uint32, Uint32> np2lcg_19;
					typedef NearP2LCG<22, 3, 916, Uint32, Uint32> np2lcg_22;
					typedef NearP2LCG<24, 3, 17880, Uint32, Uint64> np2lcg_24;
					typedef NearP2LCG<26, 5, 14791, Uint32, Uint64> np2lcg_26;
					typedef NearP2LCG<29, 3, 78944, Uint32, Uint64> np2lcg_29;
					typedef NearP2LCG<31, 1, 16807, Uint32, Uint64> np2lcg_31;
					typedef NearP2LCG<33, 9, 10349, Uint64, Uint64> np2lcg_33;
					typedef NearP2LCG<36, 5, 2027339, Uint64, Uint64> np2lcg_36;
					typedef NearP2LCG<39, 7, 42017, Uint64, Uint64> np2lcg_39;
					typedef NearP2LCG<44, 17, 89164, Uint64, Uint64> np2lcg_44;
				}// 10, 11, 13, 17, 19, 22, 24, 26, 29, 31, 33, 36, 39, 44, 
				/*
					45: 55
					46: 21
					47: 115
					48: 59
					49: 81
					50: 27
					51: 129
					52: 47
					53: 111
					54: 33
					55: 55 (cute)
					56: 5
					57: 13
				*/
				class np2clcg2_32 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_13 lcg_A;
					NP2LCG_Helpers::np2lcg_19 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_36 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_17 lcg_A;
					NP2LCG_Helpers::np2lcg_19 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_39 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_17 lcg_A;
					NP2LCG_Helpers::np2lcg_22 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_45 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_19 lcg_A;
					NP2LCG_Helpers::np2lcg_26 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_53 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_22 lcg_A;
					NP2LCG_Helpers::np2lcg_31 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_57 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_26 lcg_A;
					NP2LCG_Helpers::np2lcg_31 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_60 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_29 lcg_A;
					NP2LCG_Helpers::np2lcg_31 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_64 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_31 lcg_A;
					NP2LCG_Helpers::np2lcg_33 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_67 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_31 lcg_A;
					NP2LCG_Helpers::np2lcg_36 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg2_70 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_31 lcg_A;
					NP2LCG_Helpers::np2lcg_39 lcg_B;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class np2clcg3_40 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_10 lcg_A;
					NP2LCG_Helpers::np2lcg_13 lcg_B;
					NP2LCG_Helpers::np2lcg_17 lcg_C;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg3_49 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_13 lcg_A;
					NP2LCG_Helpers::np2lcg_17 lcg_B;
					NP2LCG_Helpers::np2lcg_19 lcg_C;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg3_58 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_17 lcg_A;
					NP2LCG_Helpers::np2lcg_19 lcg_B;
					NP2LCG_Helpers::np2lcg_22 lcg_C;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2clcg3_67 : public vRNG16 {
					NP2LCG_Helpers::np2lcg_19 lcg_A;
					NP2LCG_Helpers::np2lcg_22 lcg_B;
					NP2LCG_Helpers::np2lcg_26 lcg_C;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class np2lcg32_varqual : public vRNG32 {//LCG is limited to 32-bit constants
					Uint32 x;
					Uint32 m, a, c;
				public:
					np2lcg32_varqual(Uint32 _m, Uint32 _a, Uint32 _c);
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2lcg16_varqual : public vRNG16 {//LCG is limited to 32-bit constants ; output is lowest 16 bits of LCG state
					Uint32 x;
					Uint32 m, a, c;// mm;
				public:
					np2lcg16_varqual(Uint32 _m, Uint32 _a, Uint32 _c);
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class np2lcg8_varqual : public vRNG8 {//as above, but output is lowest 8 bits
					Uint32 x;
					Uint32 m, a, c;// mm;
				public:
					np2lcg8_varqual(Uint32 _m, Uint32 _a, Uint32 _c);
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//two LCGs combined
				class clcg8of95_varqual : public vRNG8 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg8of95_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg16of95_varqual : public vRNG16 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg16of95_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg32of95_varqual : public vRNG32 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg32of95_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg8of108_varqual : public vRNG8 {
					Uint64 p2mlcg;
					NP2LCG_Helpers::np2lcg_44 np2lcg;
					int outshift;
				public:
					enum { NP2LCG_BITS = 44 };
					clcg8of108_varqual(int p2mlcg_discard_bits) : outshift(p2mlcg_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg16of108_varqual : public vRNG16 {
					Uint64 p2mlcg;
					NP2LCG_Helpers::np2lcg_44 np2lcg;
					int outshift;
				public:
					enum { NP2LCG_BITS = 44 };
					clcg16of108_varqual(int p2mlcg_discard_bits) : outshift(p2mlcg_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg32of108_varqual : public vRNG32 {
					Uint64 p2mlcg;
					NP2LCG_Helpers::np2lcg_44 np2lcg;
					int outshift;
				public:
					enum { NP2LCG_BITS = 44 };
					clcg32of108_varqual(int p2mlcg_discard_bits) : outshift(p2mlcg_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//LCGs modified by suppressing the carries
				//note: similar to a power of 2 modulus LCG, but uses XOR instead of ADD
				// for P2M LCGs, the multiplier modulo 8 is usually supposed to be 5.  But for this variant, instead the multiplier modulo 4 should be 3.  
				// for P2M LCGs, the additive constant should be an odd number.  That remains true for this variant.  
				class xlcg32of64_varqual : public vRNG32 {
					Uint64 state;
					int outshift;
				public:
					xlcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg16of64_varqual : public vRNG16 {
					Uint64 state;
					int outshift;
				public:
					xlcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg8of64_varqual : public vRNG8 {
					Uint64 state;
					int outshift;
				public:
					xlcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg32of128_varqual : public vRNG32 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg16of128_varqual : public vRNG16 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg8of128_varqual : public vRNG8 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//modified LCG combined with regular LCG
				class cxlcg8of96_varqual : public vRNG8 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class cxlcg16of96_varqual : public vRNG16 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class cxlcg32of96_varqual : public vRNG32 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//large-state LCGs with very poor constants
				class bigbadlcg64X : public vRNG64 {
					enum { MAX_N = 16 };
					Uint64 state[MAX_N];
					int n;
				public:
					int discard_bits;
					int shift_i;
					int shift_b;
					Uint64 raw64();
					bigbadlcg64X(int discard_bits_, int shift_);
					//~bigbadlcgX();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class bigbadlcg32X : public vRNG32 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg32X(int discard_bits_, int shift_);
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class bigbadlcg16X : public vRNG16 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg16X(int discard_bits_, int shift_);
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class bigbadlcg8X : public vRNG8 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg8X(int discard_bits_, int shift_);
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class pcg32 : public vRNG32 {
					Uint64 state, inc;
				public:
					pcg32() : state(0x853c49e6748fea9bULL), inc(0xda3e39cb94b95bdbULL) {}
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class pcg32_norot : public vRNG32 {
					Uint64 state, inc;
				public:
					pcg32_norot() : state(0x853c49e6748fea9bULL), inc(0xda3e39cb94b95bdbULL) {}
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class cmrg32of192 : public vRNG32 {//I originally encountered this under the name lecuyer3by2b
					//presumably by L'Ecuyer, I adjusted it slightly to output a full 32 bits (instead of ~31.9 bits)
					//it is a Combined Multiple Recursive Generator (the moduli are 2**32-209 and 2**32-22853)
					Uint32 n1m0, n1m1, n1m2, n2m0, n2m1, n2m2;//why is one dimension zero-based and the other not?  no idea, it was that way in the code I based this off of
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class xsh_lcg_bad : public vRNG32 {//name was xorwowPlus, I changed it because I wasn't sure it actually qualified as an xorwow
					Uint64 lcg, x0, x1, x2, x3;
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);// I also changed the seeding function, because the original permitted the bad all-zeroes state
					void walk_state(StateWalkingObject *);
				};
			}
		}
	}
}
