#pragma once

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
					Uint64 state[MAX_N];
					int n;
				public:
					int discard_bits;
					int shift_i;
					int shift_b;
					Uint64 raw64() override;
					bigbadlcg64X(int discard_bits_, int shift_);
					//~bigbadlcgX();
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class bigbadlcg32X : public vRNG32 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg32X(int discard_bits_, int shift_);
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class bigbadlcg16X : public vRNG16 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg16X(int discard_bits_, int shift_);
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class bigbadlcg8X : public vRNG8 {
				public:
					bigbadlcg64X base_lcg;
					bigbadlcg8X(int discard_bits_, int shift_);
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				//Mitchell-Moore: LFib32(Uint32, 55, 24, ADD)
				class mm32 : public vRNG32 {
					Uint32 cbuf[55];
					Uint8 index1, index2;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//Mitchell-Moore modified: LFib16(Uint32, 55, 24, ADD) >> 16
				class mm16of32 : public vRNG16 {
					Uint32 cbuf[55];
					Uint8 index1, index2;
				public:
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//Mitchell-Moore modified: LFib32(Uint32, 55, 24, ADC)
				class mm32_awc : public vRNG32 {
					Uint32 cbuf[55];
					Uint8 index1, index2, carry;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//Mitchell-Moore modified: LFib16(Uint32, 55, 24, ADC)
				class mm16of32_awc : public vRNG16 {
					Uint32 cbuf[55];
					Uint8 index1, index2, carry;
				public:
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};

				class lfsr_medium : public vRNG8 {
					static constexpr int SIZE = 55;
					static constexpr int LAG = 25; //0 < LAG < SIZE-2
					Uint8 cbuf[55];
					Uint8 table1[256], table2[256];
					Uint8 used;
				public:
					lfsr_medium();
					Uint8 raw8() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};


				//proposed by Marsaglia
				class mwc4691 : public vRNG32 {
					Uint32 cbuf[4691];
					unsigned int index, carry;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				//proposed by Marsaglia
				//class cwsb4288;

				class cbuf_accum : public vRNG32 {
					static constexpr int L = 32;
					Uint32 cbuf[L], accum;
					Uint8 index;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class cbuf_accum_big : public vRNG32 {
					static constexpr int L = 128;
					Uint32 cbuf[L], accum;
					Uint32 index;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class cbuf_2accum_small : public vRNG32 {
					static constexpr int L = 3;
					Uint32 cbuf[L], accum1, accum2;
					Uint8 index;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class cbuf_2accum : public vRNG32 {
					static constexpr int L = 12;
					Uint32 cbuf[L], accum1, accum2;
					Uint8 index;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class dual_cbuf_small : public vRNG32 {
					static constexpr int L1 = 3;
					static constexpr int L2 = 5;
					Uint32 cbuf1[L1], cbuf2[L2];
					Uint8 index1, index2;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class dual_cbuf : public vRNG32 {
					static constexpr int L1 = 13;
					static constexpr int L2 = 19;
					Uint32 cbuf1[L1], cbuf2[L2];
					Uint8 index1, index2;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class dual_cbufa_small : public vRNG32 {
					static constexpr int L1 = 4;
					static constexpr int L2 = 5;
					Uint32 cbuf1[L1], cbuf2[L2], accum;
					Uint8 index1, index2;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class dual_cbuf_accum : public vRNG32 {
					static constexpr int L1 = 13;
					static constexpr int L2 = 19;
					Uint32 cbuf1[L1], cbuf2[L2], accum;
					Uint8 index1, index2;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot32small : public vRNG32 {
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 3;
					static constexpr int ROT1 = 9;
					static constexpr int ROT2 = 13;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > 0
					Uint8 position;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot32 : public vRNG32 {
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 9;
					static constexpr int ROT1 = 9;
					static constexpr int ROT2 = 13;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > 0
					Uint8 position;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot32big : public vRNG32 {
					static constexpr int LAG1 = 57;
					static constexpr int LAG2 = 13;
					static constexpr int ROT1 = 9;
					static constexpr int ROT2 = 13;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > 0
					Uint8 position;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot3tap32small : public vRNG32 {
					//7,3:29, 9,4:33, 11,5:34, 13,6:34, 15,7:35, 17,9:38
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 3;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > LAG3, LAG3 = 1
					Uint8 position;
					static Uint32 func(Uint32 a, Uint32 b, Uint32 c);
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot3tap32 : public vRNG32 {
					//7,3:29, 9,4:33, 11,5:34, 13,6:34, 15,7:35, 17,9:38
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 9;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > LAG3, LAG3 = 1
					Uint8 position;
					static Uint32 func(Uint32 a, Uint32 b, Uint32 c);
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot3tap32big : public vRNG32 {
					//7,3:29, 9,4:33, 11,5:34, 13,6:34, 15,7:35, 17,9:38
					static constexpr int LAG1 = 57;
					static constexpr int LAG2 = 13;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > LAG3, LAG3 = 1
					Uint8 position;
					static Uint32 func(Uint32 a, Uint32 b, Uint32 c);
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot32hetsmall : public vRNG32 {
					//7,3:32, 9,4:36, 11,5:37, 13,6:38-, 15,6:38, 17,9:40?
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 4;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > LAG3, LAG3 = 1
					Uint8 position;
					static Uint32 func(Uint32 a, Uint32 b, Uint32 c);
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot32het : public vRNG32 {
					//7,3:32, 9,4:36, 11,5:37, 13,6:38-, 15,6:38, 17,9:40?
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 9;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > LAG3, LAG3 = 1
					Uint8 position;
					static Uint32 func(Uint32 a, Uint32 b, Uint32 c);
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class ranrot32hetbig : public vRNG32 {
					//7,3:32, 9,4:36, 11,5:37, 13,6:38-, 15,6:38, 17,9:40?
					static constexpr int LAG1 = 57;
					static constexpr int LAG2 = 13;
					static constexpr int LAG3 = 1;
					static constexpr int ROT1 = 3;
					static constexpr int ROT2 = 17;
					static constexpr int ROT3 = 9;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > LAG3, LAG3 = 1
					Uint8 position;
					static Uint32 func(Uint32 a, Uint32 b, Uint32 c);
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class fibmul16of32 : public vRNG16 {// 31 @ 17/9
					static constexpr int LAG1 = 17;
					static constexpr int LAG2 = 5;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > 0
					Uint8 position;
				public:
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class fibmul32of64 : public vRNG32 {// 35 @ 3/2, 39 @ 7/5
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 5;
					Uint16 buffer[LAG1]; // LAG1 > LAG2 > 0
					Uint8 position;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class fibmulmix16 : public vRNG16 {
					static constexpr int LAG1 = 7;
					static constexpr int LAG2 = 3;
					Uint32 buffer[LAG1]; // LAG1 > LAG2 > 0
					Uint8 position;
				public:
					Uint16 raw16() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
				class mt19937_unhashed : public vRNG32 {//
					PractRand::RNGs::Raw::mt19937 implementation;
				public:
					Uint32 raw32() override;
					std::string get_name() const override;
					void walk_state(StateWalkingObject *) override;
				};
}
