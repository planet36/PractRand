#pragma once

#include "PractRand/rng_helpers.h"

/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that:
1. use multiplication
2. don't use much indirection, flow control, variable shifts, etc
3. have only a few words of state
4. are likely to have easily detectable bias
*/

namespace PractRand::RNGs::Polymorphic::NotRecommended {
				//similar to the classic LCGs, but with a longer period
				class lcg16of32_extended : public vRNG16 {
					Uint32 state, add;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class lcg32_extended : public vRNG32 {
					Uint32 state, add;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//simple classic LCGs
				class lcg32of64_varqual : public vRNG32 {
					Uint64 state;
					int outshift;
				public:
					lcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class lcg16of64_varqual : public vRNG16 {
					Uint64 state;
					int outshift;
				public:
					lcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class lcg8of64_varqual : public vRNG8 {
					Uint64 state;
					int outshift;
				public:
					lcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class lcg32of128_varqual : public vRNG32 {
					Uint64 low, high;
					int outshift;
				public:
					lcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class lcg16of128_varqual : public vRNG16 {
					Uint64 low, high;
					int outshift;
				public:
					lcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class lcg8of128_varqual : public vRNG8 {
					Uint64 low, high;
					int outshift;
				public:
					lcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//two LCGs combined
				class clcg8of96_varqual : public vRNG8 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class clcg16of96_varqual : public vRNG16 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class clcg32of96_varqual : public vRNG32 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//LCGs modified by suppressing the carries
				class xlcg32of64_varqual : public vRNG32 {
					Uint64 state;
					int outshift;
				public:
					xlcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class xlcg16of64_varqual : public vRNG16 {
					Uint64 state;
					int outshift;
				public:
					xlcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class xlcg8of64_varqual : public vRNG8 {
					Uint64 state;
					int outshift;
				public:
					xlcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class xlcg32of128_varqual : public vRNG32 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class xlcg16of128_varqual : public vRNG16 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class xlcg8of128_varqual : public vRNG8 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//modified LCG combined with regular LCG
				class cxlcg8of96_varqual : public vRNG8 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class cxlcg16of96_varqual : public vRNG16 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class cxlcg32of96_varqual : public vRNG32 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				class pcg32 final : public vRNG32 {
					Uint64 state{0x853c49e6748fea9bULL}, inc{0xda3e39cb94b95bdbULL};
				public:
					pcg32() = default;
					Uint32 raw32() override;
					std::string get_name() const override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
				};
				class pcg32_norot final : public vRNG32 {
					Uint64 state{0x853c49e6748fea9bULL}, inc{0xda3e39cb94b95bdbULL};
				public:
					pcg32_norot() = default;
					Uint32 raw32() override;
					std::string get_name() const override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
				};
				class cmrg32of192 : public vRNG32 {//I originally encountered this under the name lecuyer3by2b
					//presumably by L'Ecuyer, I adjusted it slightly to output a full 32 bits (instead of ~31.9 bits)
					//it is a Combined Multiple Recursive Generator (the moduli are 2**32-209 and 2**32-22853)
					Uint32 n1m0, n1m1, n1m2, n2m0, n2m1, n2m2;//why is one dimension zero-based and the other not?  no idea, it was that way in the code I based this off of
					Uint32 raw32() override;
					std::string get_name() const override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
				};
				class xsh_lcg_bad final : public vRNG32 {//name was xorwowPlus, I changed it because I wasn't sure it actually qualified as an xorwow
					Uint64 lcg, x0, x1, x2, x3;
					Uint32 raw32() override;
					std::string get_name() const override;
					void seed(Uint64 s) override;// I also changed the seeding function, because the original permitted the bad all-zeroes state
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
				};



				//
				class garthy16 : public vRNG16 {
					Uint16 value, scale, counter;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class garthy32 : public vRNG32 {
					Uint32 value, scale, counter;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//both sides of the multiply are pseudo-random values in this RNG
				class binarymult16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class binarymult32 : public vRNG32 {
					Uint32 a, b, c, d;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				class mmr16 : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mmr32 : public vRNG32 {
					Uint32 a, b, c;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				//uses multiplication, rightshifts, xors, that kind of stuff
				class rxmult16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				//these are similar to my mwlac algorithm, but lower quality
				class multish2x64 : public vRNG64 {
					Uint64 a, b;
					Uint64 raw64() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class multish3x32 : public vRNG32 {
					Uint32 a, b, c;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class multish4x16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				class mwrca16 : public vRNG16 {
					Uint16 a, b;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwrca32 : public vRNG32 {
					Uint32 a, b;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwrcc16 : public vRNG16 {
					Uint16 a, b, counter;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwrcc32 : public vRNG32 {
					Uint32 a, b, counter;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwrcca16 : public vRNG16 {
					Uint16 a, b, counter;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwrcca32 : public vRNG32 {
					Uint32 a, b, counter;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//the 16 bit variant of the old version of my mwlac algorithm
				class old_mwlac16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwlac_varA : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwlac_varB : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwlac_varC : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwlac_varD : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mwlac_varE : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};


				class mwc64x : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class cxm64_varqual : public vRNG64 {
					Uint64 low, high;
					int num_mult;
				public:
					cxm64_varqual(int num_mult_) : num_mult(num_mult_) {}
					Uint64 raw64() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				class mo_Cmfr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mo_Cmr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mo_Cmr32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class murmlac32 : public vRNG32 {
					Uint32 state1, state2;
					int rounds;
				public:
					Uint32 raw32() override;
					murmlac32(int rounds_) : rounds(rounds_) {}
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				//multiplication (by a counter), rotate
				class mulcr64 : public vRNG64 {
					Uint64 a, b, count;
					Uint64 raw64() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mulcr32 : public vRNG32 {
					Uint32 a, b, count;
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mulcr16 : public vRNG16 {
					Uint32 a, b, count;
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
}
