#pragma once

#include "PractRand/rng_helpers.h"

/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that:
1. use an array with dynamic access patterns
2. don't use complex flow control or advanced math functions (sin/cos, log/exp, etc)
3. are likely to have detectable bias
*/

namespace PractRand::RNGs::Polymorphic::NotRecommended {
				//the classic cryptographic RNG
				class rc4 : public vRNG8 {
				protected:
					uint8_t arr[256]{};
					uint8_t a{}, b{};
				public:
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					//void seed (uint64_t s);
					void walk_state(StateWalkingObject*) override;
				};

				//weaker variant of the classic cryptographic RNG
				//(same state transitions, different output)
				class rc4_weakenedA : public rc4 {
				public:
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class rc4_weakenedB : public rc4 {
				public:
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class rc4_weakenedC : public rc4 {
				public:
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class rc4_weakenedD : public rc4 {
				public:
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};


				//based upon the IBAA algorithm by Robert Jenkins
				//but reduced strength via smaller integers and smaller table sizes
				class ibaa8 : public vRNG8 {
					int table_size_L2;
					uint8_t* table;
					uint8_t a{}, b{}, left{};
				public:
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
					explicit ibaa8(int table_size_L2_);
					~ibaa8() override;
				};
				class ibaa16 : public vRNG16 {
					int table_size_L2;
					uint16_t* table;
					uint16_t a{}, b{}, left{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
					explicit ibaa16(int table_size_L2_);
					~ibaa16() override;
				};
				class ibaa32 : public vRNG32 {
					int table_size_L2;
					uint32_t* table;
					uint32_t a{}, b{}, left{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
					explicit ibaa32(int table_size_L2_);
					~ibaa32() override;
				};
				//based upon the ISAAC algorithm by Robert Jenkins
				//but adapted slightly to permit smaller minimum table sizes
				class isaac32_varqual : public vRNG32 {
					int table_size_L2;
					uint32_t* table;
					uint32_t a{}, b{}, c{}, left{};
				public:
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
					explicit isaac32_varqual(int table_size_L2_);
					~isaac32_varqual() override;
				};
				//as isaac32_small, but adapted to use 16 bit integers instead of 32
				class isaac16_varqual : public vRNG16 {
					int table_size_L2;
					uint16_t* table;
					uint16_t a{}, b{}, c{}, left{};
				public:
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
					void walk_state(StateWalkingObject*) override;
					explicit isaac16_varqual(int table_size_L2_);
					~isaac16_varqual() override;
				};
				class efiix8_varqual : public vRNG8 {
					static constexpr int SHIFT_AMOUNT = 2;
					//average cycle length, in bytes: 2**(4 * W - 2), where W=(4 + iteration_table_size + indirection_table_size)
					int iteration_table_size_m1;
					int indirection_table_size_m1;
					uint8_t* iteration_table;
					uint8_t* indirection_table;
					uint8_t a{}, b{}, c{}, i{};
				public:
					uint8_t raw8() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
					efiix8_varqual(int iteration_table_size_L2, int indirection_table_size_L2);
					~efiix8_varqual() override;
				};
				//efiix algorithm, shrunk down to operate on 4 bit integers (reports itself as an 8 bit PRNG)
				class efiix4_varqual : public vRNG8 {
					//average cycle length, in bytes: 2**(4 * W - 2), where W=(4 + iteration_table_size + indirection_table_size)
					int iteration_table_size_m1;
					int indirection_table_size_m1;
					uint8_t* iteration_table;
					uint8_t* indirection_table;
					uint8_t a{}, b{}, c{}, i{};
					uint8_t raw4();
					static uint8_t rotate4(uint8_t value, int bits);
				public:
					uint8_t raw8() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
					efiix4_varqual(int iteration_table_size_L2, int indirection_table_size_L2);
					~efiix4_varqual() override;
				};

				//generic indirection-based PRNGs, just for testing test suites
				class genindA : public vRNG16 {
					// 4:25, 5:27, 6:32, 7:38, 8:41, 9:46
					uint16_t* table;
					uint16_t a{}, i{};
					int table_size_mask;
					int shift;
				public:
					explicit genindA(int size_L2);
					~genindA() override;
					uint16_t raw16() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
				};
				class genindB : public vRNG16 {
					//efiix-style, single pool, irreversible
					// 0*:34, 1:30, 2:34, 3:39, 4:42
					int shift;
					int table_size_mask;
					uint16_t* table;
					uint16_t a{}, b{}, i{};
				public:
					explicit genindB(int size_L2);
					~genindB() override;
					uint16_t raw16() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
				};
				class genindC : public vRNG16 {
					//split pool, bits fed into the accumulator is limited to tsL2 per output
					// 1*:23, 2:31, 3:36, 4:42
					uint16_t* table;
					int16_t left;
					uint16_t a{};
					int table_size_L2;
				public:
					uint16_t raw16() override {
						if (left >= 0) return table[left--];
						else return refill();
					}
					explicit genindC(int size_L2);
					~genindC() override;
					uint16_t refill();
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
				};
				class genindD : public vRNG16 {
					//fairly solid, RC4-style
					// 2:25, 3:27, 4:30, 5:32, 6:35, 7:37, 8:39, 9:41, 10:43
					uint16_t* table;
					uint16_t a{}, i{};
					uint16_t mask;
					uint16_t table_size_L2;
				public:
					uint16_t raw16() override;
					explicit genindD(int size_L2);
					~genindD() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
				};
				class genindE : public vRNG16 {
					//split pool
					// 1:33, 2:39, 3:46
					uint16_t* table1;
					uint16_t* table2;
					uint16_t a{}, i{};
					uint16_t mask;
					uint16_t table_size_L2;
				public:
					uint16_t raw16() override;
					explicit genindE(int size_L2);
					~genindE() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
				};
				class genindF : public vRNG16 {
					//split pool
					// 1:28, 2:34, 3:37, 4:40, 5:45
					uint16_t* table1;
					uint16_t* table2;
					uint16_t a{}, i{};
					uint16_t mask;
					uint16_t table_size_L2;
				public:
					uint16_t raw16() override;
					explicit genindF(int size_L2);
					~genindF() override;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
				};
}
