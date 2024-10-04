#pragma once

#include "PractRand/rng_helpers.h"

#include <deque>
#include <vector>
//RNGs in the "other" directory are not intended for real world use
//only for research; as such they may get pretty sloppy in some areas
//and are usually not optimized
namespace PractRand::RNGs::Polymorphic::NotRecommended {
				class Transform64 : public vRNG64 {
				public:
					vRNG *base_rng;
					void seed(Uint64 seed) override;
					using vRNG::seed;
					Uint64 get_flags() const override;
					void walk_state(StateWalkingObject *walker) override;
					Transform64(vRNG *rng) : base_rng(rng) {}
					~Transform64() override;
				};
				class Transform32 : public vRNG32 {
				public:
					vRNG *base_rng;
					void seed(Uint64 seed) override;
					using vRNG::seed;
					Uint64 get_flags() const override;
					void walk_state(StateWalkingObject *walker) override;
					Transform32(vRNG *rng) : base_rng(rng) {}
					~Transform32() override;
				};
				class Transform16 : public vRNG16 {
				public:
					vRNG *base_rng;
					void seed(Uint64 seed) override;
					using vRNG::seed;
					Uint64 get_flags() const override;
					void walk_state(StateWalkingObject *walker) override;
					Transform16(vRNG *rng) : base_rng(rng) {}
					~Transform16() override;
				};
				class Transform8 : public vRNG8 {
				public:
					vRNG *base_rng;
					void seed(Uint64 seed) override;
					using vRNG::seed;
					Uint64 get_flags() const override;
					void walk_state(StateWalkingObject *walker) override;
					Transform8(vRNG *rng) : base_rng(rng) {}
					~Transform8() override;
				};
				class MultiplexTransformRNG : public vRNG {
				public:
					PractRand::Tests::TestBlock *buffer;
					int index;
					virtual void refill();
					std::vector<vRNG *> source_rngs;
				public:
					Uint8 raw8() override;
					Uint16 raw16() override;
					Uint32 raw32() override;
					Uint64 raw64() override;
					void seed(Uint64 seedval) override;
					void seed(vRNG *seeder) override;
					Uint64 get_flags() const override;
					void walk_state(StateWalkingObject *walker) override;
					MultiplexTransformRNG(const std::vector<vRNG*> &sources);
					~MultiplexTransformRNG() override;
					int get_native_output_size() const override;
				};

				class GeneralizedTableTransform : public vRNG8 {//written for self-shrinking-generators, but also useful for other transforms
				public:
					struct Entry {
						Uint8 data;
						Uint8 count;
					};
					//the transform table
					const Entry *table;

					//for buffering fractional bytes of output
					Uint32 buf_data{};
					Uint32 buf_count{};

					//for buffering full bytes of output
					std::deque<Uint8> finished_bytes;

					std::string name;
					vRNG *base_rng;
					void seed(Uint64 seed) override;
					using vRNG::seed;
					Uint64 get_flags() const override;
					std::string get_name() const override;
					GeneralizedTableTransform(vRNG *rng, const Entry *table_, std::string name_) : table(table_), name(name_), base_rng(rng) {}
					~GeneralizedTableTransform() override;
					void walk_state(StateWalkingObject *) override;
					Uint8 raw8() override;
				};
				vRNG *apply_SelfShrinkTransform(vRNG *base_rng);
				//vRNG *apply_SimpleShrinkTransform(vRNG *base_rng);

				class ReinterpretAsUnknown : public Transform8 {
					Uint8 *buffer;//don't feel like requiring a header for TestBlock
					int index;
					void refill();
				public:
					ReinterpretAsUnknown( vRNG *rng );
					~ReinterpretAsUnknown() override;
					Uint8 raw8() override;
					//to do: fix endianness issues
					std::string get_name() const override;
					int get_native_output_size() const override {return -1;}
				};
				class ReinterpretAs8 : public Transform8 {
					Uint8 *buffer;//don't feel like requiring a header for TestBlock
					int index;
					void refill();
				public:
					ReinterpretAs8( vRNG *rng );
					~ReinterpretAs8() override;
					Uint8 raw8() override;
					std::string get_name() const override;
				};
				class ReinterpretAs16 : public Transform16 {
					Uint16 *buffer;
					int index;
					void refill();
				public:
					ReinterpretAs16( vRNG *rng );
					~ReinterpretAs16() override;
					Uint16 raw16() override;
					std::string get_name() const override;
				};
				class ReinterpretAs32 : public Transform32 {
					Uint32 *buffer;
					int index;
					void refill();
				public:
					ReinterpretAs32( vRNG *rng );
					~ReinterpretAs32() override;
					Uint32 raw32() override;
					std::string get_name() const override;
				};
				class ReinterpretAs64 : public Transform64 {
					Uint64 *buffer;
					int index;
					void refill();
				public:
					ReinterpretAs64( vRNG *rng );
					~ReinterpretAs64() override;
					Uint64 raw64() override;
					std::string get_name() const override;
				};

				class Xor : public MultiplexTransformRNG {
					void refill() override;
				public:
					Xor(const std::vector<vRNG*> &sources) : MultiplexTransformRNG(sources) {}
					std::string get_name() const override;
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
					typedef Uint16 InWord;
					typedef Uint8 OutWord;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord *buffer;
					int index;
					void refill();
				public:
					Discard16to8(vRNG *rng);
					Uint8 raw8() override;
					std::string get_name() const override;
				};
				class Discard32to8 : public Transform8 {
					typedef Uint32 InWord;
					typedef Uint8 OutWord;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord *buffer;
					int index;
					void refill();
				public:
					Discard32to8(vRNG *rng);
					Uint8 raw8() override;
					std::string get_name() const override;
				};
				class Discard64to8 : public Transform8 {
					typedef Uint64 InWord;
					typedef Uint8 OutWord;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord *buffer;
					int index;
					void refill();
				public:
					Discard64to8(vRNG *rng);
					Uint8 raw8() override;
					std::string get_name() const override;
				};
				class Discard32to16 : public Transform16 {
					typedef Uint32 InWord;
					typedef Uint16 OutWord;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord *buffer;
					int index;
					void refill();
				public:
					Discard32to16(vRNG *rng);
					Uint16 raw16() override;
					std::string get_name() const override;
				};
				class Discard64to16 : public Transform16 {
					typedef Uint64 InWord;
					typedef Uint16 OutWord;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord *buffer;
					int index;
					void refill();
				public:
					Discard64to16(vRNG *rng);
					Uint16 raw16() override;
					std::string get_name() const override;
				};
				class Discard64to32 : public Transform32 {
					typedef Uint64 InWord;
					typedef Uint32 OutWord;
					static constexpr int INPUT_BITS = sizeof(InWord)* 8;
					InWord *buffer;
					int index;
					void refill();
				public:
					Discard64to32(vRNG *rng);
					Uint32 raw32() override;
					std::string get_name() const override;
				};

				class BaysDurhamShuffle64 final : public Transform64 {
					Uint64 table[256];
					Uint8 prev;
					Uint8 index_mask;
					Uint8 index_shift;
				public:
					Uint64 raw64() override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
					std::string get_name() const override;
					BaysDurhamShuffle64(vRNG64 *rng, int table_size_L2, int shift=0)
						: Transform64(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				class BaysDurhamShuffle32 final : public Transform32 {
					Uint32 table[256];
					Uint8 prev;
					Uint8 index_mask;
					Uint8 index_shift;
				public:
					Uint32 raw32() override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
					std::string get_name() const override;
					BaysDurhamShuffle32(vRNG32 *rng, int table_size_L2, int shift=0)
						: Transform32(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				class BaysDurhamShuffle16 final : public Transform16 {
					Uint16 table[256];
					Uint8 prev;
					Uint8 index_mask;
					Uint8 index_shift;
				public:
					Uint16 raw16() override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
					std::string get_name() const override;
					BaysDurhamShuffle16(vRNG16 *rng, int table_size_L2, int shift=0)
						: Transform16(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				class BaysDurhamShuffle8 final : public Transform8 {
					Uint8 table[256];
					Uint8 prev;
					Uint8 index_mask;
					Uint8 index_shift;
				public:
					Uint8 raw8() override;
					void seed(Uint64 s) override;
					using vRNG::seed;
					void walk_state(StateWalkingObject *) override;
					std::string get_name() const override;
					BaysDurhamShuffle8(vRNG8 *rng, int table_size_L2, int shift=0)
						: Transform8(rng), index_mask((1<<table_size_L2)-1), index_shift(shift) {}
				};
				vRNG *apply_BaysDurhamShuffle(vRNG *base_rng, int table_size_L2=8, int shift=-1);
}
