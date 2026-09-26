#pragma once

#include "PractRand/RNGs/mt19937.h"
#include "PractRand/rng_helpers.h"

/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that:
1. use an array with repetitive access patterns - generally a Fibonacci-style cyclic buffer
2. don't use much flow control, variable shifts, etc
3. are likely to have easily detectable bias
*/

namespace PractRand::RNGs::Polymorphic::NotRecommended {
				//large-state LCGs with very poor constants
				class bigbadlcg64X : public vRNG64 {
					static constexpr int MAX_N = 16;
					uint64_t state[MAX_N]{};
					int n;
				public:
					int discard_bits;
					int shift_i;
					int shift_b;
					uint64_t raw64() override;
					bigbadlcg64X(int discard_bits_, int shift_);
					//~bigbadlcgX();
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class bigbadlcg32X : public vRNG32 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg32X(int discard_bits_, int shift_);
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class bigbadlcg16X : public vRNG16 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg16X(int discard_bits_, int shift_);
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class bigbadlcg8X : public vRNG8 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg8X(int discard_bits_, int shift_);
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				//Mitchell-Moore: LFib32(uint32_t, 55, 24, ADD)
				class mm32 : public vRNG32 {
					uint32_t cbuf[55]{};
					uint8_t index1{}, index2{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//Mitchell-Moore modified: LFib16(uint32_t, 55, 24, ADD) >> 16
				class mm16of32 : public vRNG16 {
					uint32_t cbuf[55]{};
					uint8_t index1{}, index2{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//Mitchell-Moore modified: LFib32(uint32_t, 55, 24, ADC)
				class mm32_awc : public vRNG32 {
					uint32_t cbuf[55]{};
					uint8_t index1{}, index2{}, carry{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//Mitchell-Moore modified: LFib16(uint32_t, 55, 24, ADC)
				class mm16of32_awc : public vRNG16 {
					uint32_t cbuf[55]{};
					uint8_t index1{}, index2{}, carry{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};

				class lfsr_medium : public vRNG8 {
					static constexpr int SIZE = 55;
					static constexpr int LAG = 25; //0 < LAG < SIZE-2
					uint8_t cbuf[55]{};
					uint8_t table1[256]{}, table2[256]{};
					uint8_t used{};
				public:
					lfsr_medium();
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};


				//proposed by Marsaglia
				class mwc4691 : public vRNG32 {
					uint32_t cbuf[4691]{};
					unsigned int index{}, carry{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				//proposed by Marsaglia
				//class cwsb4288;

				class cbuf_accum : public vRNG32 {
					static constexpr int L = 32;
					uint32_t cbuf[L]{}, accum{};
					uint8_t index{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class cbuf_accum_big : public vRNG32 {
					static constexpr int L = 128;
					uint32_t cbuf[L]{}, accum{};
					uint32_t index{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class cbuf_2accum_small : public vRNG32 {
					static constexpr int L = 3;
					uint32_t cbuf[L]{}, accum1{}, accum2{};
					uint8_t index{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class cbuf_2accum : public vRNG32 {
					static constexpr int L = 12;
					uint32_t cbuf[L]{}, accum1{}, accum2{};
					uint8_t index{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class dual_cbuf_small : public vRNG32 {
					static constexpr int L1 = 3;
					static constexpr int L2 = 5;
					uint32_t cbuf1[L1]{}, cbuf2[L2]{};
					uint8_t index1{}, index2{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class dual_cbuf : public vRNG32 {
					static constexpr int L1 = 13;
					static constexpr int L2 = 19;
					uint32_t cbuf1[L1]{}, cbuf2[L2]{};
					uint8_t index1{}, index2{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class dual_cbufa_small : public vRNG32 {
					static constexpr int L1 = 4;
					static constexpr int L2 = 5;
					uint32_t cbuf1[L1]{}, cbuf2[L2]{}, accum{};
					uint8_t index1{}, index2{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class dual_cbuf_accum : public vRNG32 {
					static constexpr int L1 = 13;
					static constexpr int L2 = 19;
					uint32_t cbuf1[L1]{}, cbuf2[L2]{}, accum{};
					uint8_t index1{}, index2{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot32small : public vRNG32 {
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 3;
					static constexpr int ROT1 = 9;
					static constexpr int ROT2 = 13;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > 0
					uint8_t position{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot32 : public vRNG32 {
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 9;
					static constexpr int ROT1 = 9;
					static constexpr int ROT2 = 13;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > 0
					uint8_t position{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot32big : public vRNG32 {
					static constexpr int LAG1 = 57;
					static constexpr int LAG2 = 13;
					static constexpr int ROT1 = 9;
					static constexpr int ROT2 = 13;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > 0
					uint8_t position{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot3tap32small : public vRNG32 {
					//7,3:29, 9,4:33, 11,5:34, 13,6:34, 15,7:35, 17,9:38
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 3;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > LAG3, LAG3 = 1
					uint8_t position{};
					static uint32_t func(uint32_t a, uint32_t b, uint32_t c);
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot3tap32 : public vRNG32 {
					//7,3:29, 9,4:33, 11,5:34, 13,6:34, 15,7:35, 17,9:38
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 9;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > LAG3, LAG3 = 1
					uint8_t position{};
					static uint32_t func(uint32_t a, uint32_t b, uint32_t c);
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot3tap32big : public vRNG32 {
					//7,3:29, 9,4:33, 11,5:34, 13,6:34, 15,7:35, 17,9:38
					static constexpr int LAG1 = 57;
					static constexpr int LAG2 = 13;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > LAG3, LAG3 = 1
					uint8_t position{};
					static uint32_t func(uint32_t a, uint32_t b, uint32_t c);
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot32hetsmall : public vRNG32 {
					//7,3:32, 9,4:36, 11,5:37, 13,6:38-, 15,6:38, 17,9:40?
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 4;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > LAG3, LAG3 = 1
					uint8_t position{};
					static uint32_t func(uint32_t a, uint32_t b, uint32_t c);
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot32het : public vRNG32 {
					//7,3:32, 9,4:36, 11,5:37, 13,6:38-, 15,6:38, 17,9:40?
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 9;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > LAG3, LAG3 = 1
					uint8_t position{};
					static uint32_t func(uint32_t a, uint32_t b, uint32_t c);
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class ranrot32hetbig : public vRNG32 {
					//7,3:32, 9,4:36, 11,5:37, 13,6:38-, 15,6:38, 17,9:40?
					static constexpr int LAG1 = 57;
					static constexpr int LAG2 = 13;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > LAG3, LAG3 = 1
					uint8_t position{};
					static uint32_t func(uint32_t a, uint32_t b, uint32_t c);
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class fibmul16of32 : public vRNG16 {// 31 @ 17/9
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 5;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > 0
					uint8_t position{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class fibmul32of64 : public vRNG32 {// 35 @ 3/2, 39 @ 7/5
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 5;
					uint16_t buffer[LAG1]{}; // LAG1 > LAG2 > 0
					uint8_t position{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class fibmulmix16 : public vRNG16 {
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 3;
					uint32_t buffer[LAG1]{}; // LAG1 > LAG2 > 0
					uint8_t position{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
				class mt19937_unhashed : public vRNG32 {//
					PractRand::RNGs::Raw::mt19937 implementation{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
				};
}
