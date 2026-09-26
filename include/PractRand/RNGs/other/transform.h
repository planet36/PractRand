#pragma once

#include "PractRand/rng_helpers.h"

#include <deque>
#include <utility>
#include <vector>
//RNGs in the "other" directory are not intended for real world use
//only for research; as such they may get pretty sloppy in some areas
//and are usually not optimized
namespace PractRand::RNGs::Polymorphic::NotRecommended {
				class Transform64 : public vRNG64 {
				public:
					vRNG* base_rng;
					void seed(uint64_t seed) override;
					using vRNG::seed;
					[[nodiscard]] uint64_t get_flags() const override;
					void walk_state(StateWalkingObject* walker) override;
					explicit Transform64(vRNG* rng) : base_rng(rng) {}
					~Transform64() override;
				};
				class Transform32 : public vRNG32 {
				public:
					vRNG* base_rng;
					void seed(uint64_t seed) override;
					using vRNG::seed;
					[[nodiscard]] uint64_t get_flags() const override;
					void walk_state(StateWalkingObject* walker) override;
					explicit Transform32(vRNG* rng) : base_rng(rng) {}
					~Transform32() override;
				};
				class Transform16 : public vRNG16 {
				public:
					vRNG* base_rng;
					void seed(uint64_t seed) override;
					using vRNG::seed;
					[[nodiscard]] uint64_t get_flags() const override;
					void walk_state(StateWalkingObject* walker) override;
					explicit Transform16(vRNG* rng) : base_rng(rng) {}
					~Transform16() override;
				};
				class Transform8 : public vRNG8 {
				public:
					vRNG* base_rng;
					void seed(uint64_t seed) override;
					using vRNG::seed;
					[[nodiscard]] uint64_t get_flags() const override;
					void walk_state(StateWalkingObject* walker) override;
					explicit Transform8(vRNG* rng) : base_rng(rng) {}
					~Transform8() override;
				};
				class MultiplexTransformRNG : public vRNG {
				public:
					PractRand::Tests::TestBlock* buffer;
					int index{999999};
					virtual void refill();
					std::vector<vRNG*> source_rngs;

					uint8_t raw8() override;
					uint16_t raw16() override;
					uint32_t raw32() override;
					uint64_t raw64() override;
					void seed(uint64_t seedval) override;
					void seed(vRNG* seeder) override;
					[[nodiscard]] uint64_t get_flags() const override;
					void walk_state(StateWalkingObject* walker) override;
					explicit MultiplexTransformRNG(const std::vector<vRNG*>& sources);
					~MultiplexTransformRNG() override;
					[[nodiscard]] int get_native_output_size() const override;
				};

				class GeneralizedTableTransform : public vRNG8 {//written for self-shrinking-generators, but also useful for other transforms
				public:
					struct Entry {
						uint8_t data;
						uint8_t count;
					};
					//the transform table
					const Entry* table;

					//for buffering fractional bytes of output
					uint32_t buf_data{};
					uint32_t buf_count{};

					//for buffering full bytes of output
					std::deque<uint8_t> finished_bytes;

					std::string name;
					vRNG* base_rng;
					void seed(uint64_t seed) override;
					using vRNG::seed;
					[[nodiscard]] uint64_t get_flags() const override;
					[[nodiscard]] std::string get_name() const override;
					GeneralizedTableTransform(vRNG* rng, const Entry* table_, std::string name_) : table(table_), name(std::move(name_)), base_rng(rng) {}
					~GeneralizedTableTransform() override;
					void walk_state(StateWalkingObject*) override;
					uint8_t raw8() override;
				};
				vRNG* apply_SelfShrinkTransform(vRNG* base_rng);
				//vRNG *apply_SimpleShrinkTransform(vRNG *base_rng);

				class ReinterpretAsUnknown : public Transform8 {
					uint8_t* buffer;//don't feel like requiring a header for TestBlock
					int index{8192 / OUTPUT_BITS};
					void refill();
				public:
					explicit ReinterpretAsUnknown( vRNG* rng );
					~ReinterpretAsUnknown() override;
					uint8_t raw8() override;
					//to do: fix endianness issues
					[[nodiscard]] std::string get_name() const override;
					[[nodiscard]] int get_native_output_size() const override {return -1;}
				};
				class ReinterpretAs8 : public Transform8 {
					uint8_t* buffer;//don't feel like requiring a header for TestBlock
					int index{8192 / OUTPUT_BITS};
					void refill();
				public:
					explicit ReinterpretAs8( vRNG* rng );
					~ReinterpretAs8() override;
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class ReinterpretAs16 : public Transform16 {
					uint16_t* buffer;
					int index{8192 / OUTPUT_BITS};
					void refill();
				public:
					explicit ReinterpretAs16( vRNG* rng );
					~ReinterpretAs16() override;
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class ReinterpretAs32 : public Transform32 {
					uint32_t* buffer;
					int index{8192 / OUTPUT_BITS};
					void refill();
				public:
					explicit ReinterpretAs32( vRNG* rng );
					~ReinterpretAs32() override;
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class ReinterpretAs64 : public Transform64 {
					uint64_t* buffer;
					int index{8192 / OUTPUT_BITS};
					void refill();
				public:
					explicit ReinterpretAs64( vRNG* rng );
					~ReinterpretAs64() override;
					uint64_t raw64() override;
					[[nodiscard]] std::string get_name() const override;
				};

				class Xor : public MultiplexTransformRNG {
					void refill() override;
				public:
					explicit Xor(const std::vector<vRNG*>& sources) : MultiplexTransformRNG(sources) {}
					[[nodiscard]] std::string get_name() const override;
				};
				/*class Interleave8 : public MultiplexTransformRNG {
					virtual void refill() override;
				public:
					Interleave8(const std::vector<vRNG*> &sources) : MultiplexTransformRNG(sources) {}
					std::string get_name() const override;
				};
				class Interleave16 : public MultiplexTransformRNG {
					virtual void refill() override;
				public:
					Interleave16(const std::vector<vRNG*> &sources) : MultiplexTransformRNG(sources) {}
					std::string get_name() const override;
				};
				class Interleave32 : public MultiplexTransformRNG {
					virtual void refill() override;
				public:
					Interleave32(const std::vector<vRNG*> &sources) : MultiplexTransformRNG(sources) {}
					std::string get_name() const override;
				};
				class Interleave64 : public MultiplexTransformRNG {
					virtual void refill() override;
				public:
					Interleave64(const std::vector<vRNG*> &sources) : MultiplexTransformRNG(sources) {}
					std::string get_name() const override;
				};*/

				class Discard16to8 : public Transform8 {
					using InWord = uint16_t;
					using OutWord = uint8_t;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord* buffer;
					int index{8192 / INPUT_BITS};
					void refill();
				public:
					explicit Discard16to8(vRNG* base_rng_);
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class Discard32to8 : public Transform8 {
					using InWord = uint32_t;
					using OutWord = uint8_t;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord* buffer;
					int index{8192 / 32};
					void refill();
				public:
					explicit Discard32to8(vRNG* base_rng_);
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class Discard64to8 : public Transform8 {
					using InWord = uint64_t;
					using OutWord = uint8_t;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord* buffer;
					int index{8192 / INPUT_BITS};
					void refill();
				public:
					explicit Discard64to8(vRNG* base_rng_);
					uint8_t raw8() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class Discard32to16 : public Transform16 {
					using InWord = uint32_t;
					using OutWord = uint16_t;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord* buffer;
					int index{8192 / INPUT_BITS};
					void refill();
				public:
					explicit Discard32to16(vRNG* base_rng_);
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class Discard64to16 : public Transform16 {
					using InWord = uint64_t;
					using OutWord = uint16_t;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord* buffer;
					int index{8192 / INPUT_BITS};
					void refill();
				public:
					explicit Discard64to16(vRNG* base_rng_);
					uint16_t raw16() override;
					[[nodiscard]] std::string get_name() const override;
				};
				class Discard64to32 : public Transform32 {
					using InWord = uint64_t;
					using OutWord = uint32_t;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord* buffer;
					int index{8192 / INPUT_BITS};
					void refill();
				public:
					explicit Discard64to32(vRNG* base_rng_);
					uint32_t raw32() override;
					[[nodiscard]] std::string get_name() const override;
				};

				class BaysDurhamShuffle64 final : public Transform64 {
					uint64_t table[256]{};
					uint8_t prev{};
					uint8_t index_mask;
					uint8_t index_shift;
				public:
					uint64_t raw64() override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
					BaysDurhamShuffle64(vRNG64* rng, int table_size_L2, int shift=0)
						: Transform64(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				class BaysDurhamShuffle32 final : public Transform32 {
					uint32_t table[256]{};
					uint8_t prev{};
					uint8_t index_mask;
					uint8_t index_shift;
				public:
					uint32_t raw32() override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
					BaysDurhamShuffle32(vRNG32* rng, int table_size_L2, int shift=0)
						: Transform32(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				class BaysDurhamShuffle16 final : public Transform16 {
					uint16_t table[256]{};
					uint8_t prev{};
					uint8_t index_mask;
					uint8_t index_shift;
				public:
					uint16_t raw16() override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
					BaysDurhamShuffle16(vRNG16* rng, int table_size_L2, int shift=0)
						: Transform16(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				class BaysDurhamShuffle8 final : public Transform8 {
					uint8_t table[256]{};
					uint8_t prev{};
					uint8_t index_mask;
					uint8_t index_shift;
				public:
					uint8_t raw8() override;
					void seed(uint64_t s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject*) override;
					[[nodiscard]] std::string get_name() const override;
					BaysDurhamShuffle8(vRNG8* rng, int table_size_L2, int shift=0)
						: Transform8(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				vRNG* apply_BaysDurhamShuffle(vRNG* base_rng, int table_size_L2=8, int shift=-1);
}
