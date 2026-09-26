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
					uint32_t state{}, add{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class lcg32_extended : public vRNG32 {
					uint32_t state{}, add{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//simple classic LCGs
				class lcg32of64_varqual : public vRNG32 {
					uint64_t state{};
					int outshift;
				public:
					explicit lcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class lcg16of64_varqual : public vRNG16 {
					uint64_t state{};
					int outshift;
				public:
					explicit lcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class lcg8of64_varqual : public vRNG8 {
					uint64_t state{};
					int outshift;
				public:
					explicit lcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class lcg32of128_varqual : public vRNG32 {
					uint64_t low{}, high{};
					int outshift;
				public:
					explicit lcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class lcg16of128_varqual : public vRNG16 {
					uint64_t low{}, high{};
					int outshift;
				public:
					explicit lcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class lcg8of128_varqual : public vRNG8 {
					uint64_t low{}, high{};
					int outshift;
				public:
					explicit lcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//two LCGs combined
				class clcg8of96_varqual : public vRNG8 {
					uint64_t lcg1{};
					uint32_t lcg2{};
					int outshift;
				public:
					explicit clcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class clcg16of96_varqual : public vRNG16 {
					uint64_t lcg1{};
					uint32_t lcg2{};
					int outshift;
				public:
					explicit clcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class clcg32of96_varqual : public vRNG32 {
					uint64_t lcg1{};
					uint32_t lcg2{};
					int outshift;
				public:
					explicit clcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//LCGs modified by suppressing the carries
				class xlcg32of64_varqual : public vRNG32 {
					uint64_t state{};
					int outshift;
				public:
					explicit xlcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xlcg16of64_varqual : public vRNG16 {
					uint64_t state{};
					int outshift;
				public:
					explicit xlcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xlcg8of64_varqual : public vRNG8 {
					uint64_t state{};
					int outshift;
				public:
					explicit xlcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xlcg32of128_varqual : public vRNG32 {
					uint64_t low{}, high{};
					int outshift;
				public:
					explicit xlcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xlcg16of128_varqual : public vRNG16 {
					uint64_t low{}, high{};
					int outshift;
				public:
					explicit xlcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xlcg8of128_varqual : public vRNG8 {
					uint64_t low{}, high{};
					int outshift;
				public:
					explicit xlcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//modified LCG combined with regular LCG
				class cxlcg8of96_varqual : public vRNG8 {
					uint64_t lcg1{};
					uint32_t lcg2{};
					int outshift;
				public:
					explicit cxlcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class cxlcg16of96_varqual : public vRNG16 {
					uint64_t lcg1{};
					uint32_t lcg2{};
					int outshift;
				public:
					explicit cxlcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class cxlcg32of96_varqual : public vRNG32 {
					uint64_t lcg1{};
					uint32_t lcg2{};
					int outshift;
				public:
					explicit cxlcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class pcg32 final : public vRNG32 {
					uint64_t state{0x853c49e6748fea9bULL}, inc{0xda3e39cb94b95bdbULL};
				public:
					pcg32() = default;
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
				};
				class pcg32_norot final : public vRNG32 {
					uint64_t state{0x853c49e6748fea9bULL}, inc{0xda3e39cb94b95bdbULL};
				public:
					pcg32_norot() = default;
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
				};
				class cmrg32of192 : public vRNG32 {//I originally encountered this under the name lecuyer3by2b
					//presumably by L'Ecuyer, I adjusted it slightly to output a full 32 bits (instead of ~31.9 bits)
					//it is a Combined Multiple Recursive Generator (the moduli are 2**32-209 and 2**32-22853)
					uint32_t n1m0{}, n1m1{}, n1m2{}, n2m0{}, n2m1{}, n2m2{};//why is one dimension zero-based and the other not?  no idea, it was that way in the code I based this off of
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
				};
				class xsh_lcg_bad final : public vRNG32 {//name was xorwowPlus, I changed it because I wasn't sure it actually qualified as an xorwow
					uint64_t lcg{}, x0{}, x1{}, x2{}, x3{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void seed(uint64_t s) override;// I also changed the seeding function, because the original permitted the bad all-zeroes state
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
				};



				//
				class garthy16 : public vRNG16 {
					uint16_t value{}, scale{}, counter{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class garthy32 : public vRNG32 {
					uint32_t value{}, scale{}, counter{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//both sides of the multiply are pseudo-random values in this RNG
				class binarymult16 : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class binarymult32 : public vRNG32 {
					uint32_t a{}, b{}, c{}, d{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class mmr16 : public vRNG16 {
					uint16_t a{}, b{}, c{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mmr32 : public vRNG32 {
					uint32_t a{}, b{}, c{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//uses multiplication, rightshifts, xors, that kind of stuff
				class rxmult16 : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//these are similar to my mwlac algorithm, but lower quality
				class multish2x64 : public vRNG64 {
					uint64_t a{}, b{};
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class multish3x32 : public vRNG32 {
					uint32_t a{}, b{}, c{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class multish4x16 : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class mwrca16 : public vRNG16 {
					uint16_t a{}, b{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwrca32 : public vRNG32 {
					uint32_t a{}, b{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwrcc16 : public vRNG16 {
					uint16_t a{}, b{}, counter{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwrcc32 : public vRNG32 {
					uint32_t a{}, b{}, counter{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwrcca16 : public vRNG16 {
					uint16_t a{}, b{}, counter{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwrcca32 : public vRNG32 {
					uint32_t a{}, b{}, counter{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//the 16 bit variant of the old version of my mwlac algorithm
				class old_mwlac16 : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwlac_varA : public vRNG16 {
					uint16_t a{}, b{}, c{};//, d;
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwlac_varB : public vRNG16 {
					uint16_t a{}, b{}, c{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwlac_varC : public vRNG16 {
					uint16_t a{}, b{}, c{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwlac_varD : public vRNG16 {
					uint16_t a{}, b{}, c{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mwlac_varE : public vRNG16 {
					uint16_t a{}, b{}, c{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};


				class mwc64x : public vRNG32 {
					uint64_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class cxm64_varqual : public vRNG64 {
					uint64_t low{}, high{};
					int num_mult;
				public:
					explicit cxm64_varqual(int num_mult_) : num_mult(num_mult_) {}
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class mo_Cmfr32 : public vRNG32 {
					uint32_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mo_Cmr32 : public vRNG32 {
					uint32_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mo_Cmr32of64 : public vRNG32 {
					uint64_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class murmlac32 : public vRNG32 {
					uint32_t state1{}, state2{};
					int rounds;
				public:
					uint32_t raw32() override;
					explicit murmlac32(int rounds_) : rounds(rounds_) {}
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//multiplication (by a counter), rotate
				class mulcr64 : public vRNG64 {
					uint64_t a{}, b{}, count{};
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mulcr32 : public vRNG32 {
					uint32_t a{}, b{}, count{};
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mulcr16 : public vRNG16 {
					uint32_t a{}, b{}, count{};
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
}
