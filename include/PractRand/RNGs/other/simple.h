#pragma once

#include "PractRand/rng_helpers.h"

/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that do not make any significant use of:
	multiplication/division, arrays, flow control, complex math functions
*/

namespace PractRand::RNGs::Polymorphic::NotRecommended {
				class xsalta16x3 : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xsaltb16x3 : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xsaltc16x3 : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xsaltd16x3 : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xsalte16x3 : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//xorshift RNGs, a subset of LFSRs proposed by Marsaglia in 2003
				class xorshift32 : public vRNG32 {
					//constants those Marsaglia described as "one of my favorites" on page 4 of his 2003 paper
					uint32_t a{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorshift64 : public vRNG64 {
					//constants those Marsaglia used in his sample code on page 4 of his 2003 paper
					uint64_t a{};
				public:
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorshift64of128 : public vRNG64 {
					//the constants are still in need of tuning
					uint64_t high{}, low{};
					void xls(int bits);
					void xrs(int bits);
				public:
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorshift32of128 : public vRNG32 {
					xorshift64of128 impl;
				public:
					uint32_t raw32() override {return uint32_t(impl.raw64());}
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorshift16of32 : public vRNG16 {
					xorshift32 impl;
				public:
					uint16_t raw16() override {return uint16_t(impl.raw32());}
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorshift32of64 : public vRNG32 {
					xorshift64 impl;
				public:
					uint32_t raw32() override {return uint32_t(impl.raw64());}
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class xorshift32x4 : public vRNG32 {
					//recommended at the top of Marsaglia's 2003 xorshift paper
					uint32_t x{},y{},z{},w{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorwow32plus32 : public vRNG32 {
					uint32_t x{}, i{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorwow32of96 : public vRNG32 {
					xorshift64 impl;
					uint32_t a{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorwow32x6 : public vRNG32 {
					//recommended at the top of Marsaglia's 2003 xorshift paper
					uint32_t x{},y{},z{},w{},v{},d{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xorshift128plus : public vRNG64 {
					uint64_t state0{}, state1{};
				public:
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xoroshiro128plus : public vRNG64 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					uint64_t state0{}, state1{};
				public:
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class xoroshiro128plus_2p64 : public vRNG64 {
					// as xoroshiro128plus, but it skips 2**64-1 outputs between each pair of outputs (testing its recommended parallel sequences)
					uint64_t state0{}, state1{};
				public:
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class tinyMT : public vRNG32 {
					uint32_t state[4]{};
					uint32_t state_param1{}, state_param2{}, out_param{};
				public:
					void next_state();
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};


				//by Ilya O. Levin, see http://www.literatecode.com/2004/10/18/sapparot/
				class sapparot : public vRNG32 {
					uint32_t a{}, b{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//variants of sapparot created for testing purposes
				class sap16of48 : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class sap32of96 : public vRNG32 {
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//the low quality variant of the FLEA RNG by Robert Jenkins
				//(he named this variant flea0)
				class flea32x1 : public vRNG32 {
					static constexpr int SIZE = 1;
					uint32_t a[SIZE]{}, b{}, c{}, d{}, i{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//versions of my SFC algorithm, archived here
				class sfc_v1_16 : public vRNG16 {
					uint16_t a{}, b{}, counter{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class sfc_v1_32: public vRNG32 {
					uint32_t a{}, b{}, counter{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class sfc_v2_16 : public vRNG16 {
					uint16_t a{}, b{}, counter{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class sfc_v2_32 : public vRNG32 {
					uint32_t a{}, b{}, counter{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class sfc_v3_16 : public vRNG16 {
					uint16_t a{}, b{}, counter{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class sfc_v3_32 : public vRNG32 {
					uint32_t a{}, b{}, counter{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class jsf16 : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class tyche final : public vRNG32 {
					uint32_t a{}, b{}, c{}, d{};
				public:
					uint32_t raw32() override;
					void seed(uint64_t s) override;
					void seed(uint64_t s, uint32_t idx);
					using vRNG::seed;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class tyche16 : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//a few simple RNGs just for testing purposes
				class simpleA : public vRNG32 {
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class simpleB : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class simpleC : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class simpleD : public vRNG32 {
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class simpleE : public vRNG32 {
					//seems like a good combination of speed & quality
					//but falls flat when used on 16 bit words (irreversible, statespace issues)
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class simpleF : public vRNG16 {
					uint16_t a{}, b{}, c{}, d{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class simpleG : public vRNG32 {
					uint32_t a{}, b{}, c{}, d{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class trivium_weakenedA : public vRNG32 {
					uint64_t a{}, b{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class trivium_weakenedB : public vRNG16 {
					uint32_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//see http://www.drdobbs.com/tools/229625477
				class mo_Lesr32 : public vRNG32 {
					uint32_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mo_ResrRers32 : public vRNG32 {
					uint32_t a{};
					uint32_t b{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mo_Rers32of64 : public vRNG32 {
					uint64_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mo_Resr32of64 : public vRNG32 {
					uint64_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mo_Resdra32of64 : public vRNG32 {
					uint64_t state{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class murmlacish : public vRNG32 {
					uint32_t state1{}, state2{}, state3{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class gjishA : public vRNG16 {
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;//broadly similar to gjrand, but 16 bit and no counter
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class gjishB : public vRNG16 {
					uint16_t a{}, b{}, c{}, counter{};
				public:
					uint16_t raw16() override;//broadly similar to gjrand, but 16 bit
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class gjishC : public vRNG32 {
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;//broadly similar to gjrand, but 32 bit, no counter, and no 3rd shift
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class gjishD : public vRNG32 {
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;//broadly similar to gjrand, but 32 bit, no counter
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ara16 : public vRNG16 {//add, bit rotate, add
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ara32 : public vRNG32 {//add, bit rotate, add
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class arx16 : public vRNG16 {//add, bit rotate, xor
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class arx32 : public vRNG32 {//add, bit rotate, xor
					uint32_t a{}, b{}, c{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class hara16 : public vRNG16 {//heterogeneous add, bit rotate, xor
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class harx16 : public vRNG16 {//heterogeneous add, bit rotate, xor
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class learx16 : public vRNG16 {//LEAs, bit rotates, xor
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class hlearx16 : public vRNG16 {//heterogeneous LEAs, bit rotates, xor
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class alearx16 : public vRNG16 {//add, LEA, bit rotate, xor
					uint16_t a{}, b{}, c{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class arac16 : public vRNG16 {//add, bit rotate, add (with counter)
					uint16_t a{}, b{}, c{}, counter{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class arxc16 : public vRNG16 {//add, bit rotate, xor (with counter)
					uint16_t a{}, b{}, c{}, counter{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
}
