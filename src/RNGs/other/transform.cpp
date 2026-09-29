#include "PractRand/RNGs/other/transform.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"
#include "PractRand/tests.h"

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace PractRand::RNGs::Polymorphic::NotRecommended {
				void Transform64::seed(uint64_t s) {base_rng->seed(s);}
				uint64_t Transform64::get_flags() const {return base_rng->get_flags() | FLAG::USES_INDIRECTION;}
				void Transform64::walk_state(StateWalkingObject* walker) {base_rng->walk_state(walker);}
				Transform64::~Transform64() {delete base_rng;}
				void Transform32::seed(uint64_t s) {base_rng->seed(s);}
				uint64_t Transform32::get_flags() const {return base_rng->get_flags() | FLAG::USES_INDIRECTION;}
				void Transform32::walk_state(StateWalkingObject* walker) {base_rng->walk_state(walker);}
				Transform32::~Transform32() {delete base_rng;}
				void Transform16::seed(uint64_t s) {base_rng->seed(s);}
				uint64_t Transform16::get_flags() const {return base_rng->get_flags() | FLAG::USES_INDIRECTION;}
				void Transform16::walk_state(StateWalkingObject* walker) {base_rng->walk_state(walker);}
				Transform16::~Transform16() {delete base_rng;}
				void Transform8::seed(uint64_t s) {base_rng->seed(s);}
				uint64_t Transform8::get_flags() const {return base_rng->get_flags() | FLAG::USES_INDIRECTION;}
				void Transform8::walk_state(StateWalkingObject* walker) {base_rng->walk_state(walker);}
				Transform8::~Transform8() {delete base_rng;}
				void MultiplexTransformRNG::refill() { index = 0; }
				uint8_t MultiplexTransformRNG::raw8() {
					if (index >= Tests::TestBlock::SIZE) {
						refill();
						index = 1;
						return buffer->as8[0];
					}
					return buffer->as8[index++];
				}
				uint16_t MultiplexTransformRNG::raw16() {
					index += 3; index &= ~1;//round up to force alignment, and also increment position
					if (index > Tests::TestBlock::SIZE) {
						refill();
						index = 2;
						return buffer->as16[0];
					}
					uint16_t rv = *reinterpret_cast<uint16_t*>(&buffer->as8[index - 2]);//read 16 aligned bits
					return rv;
				}
				uint32_t MultiplexTransformRNG::raw32() {
					index += 7; index &= ~3;
					if (index > Tests::TestBlock::SIZE) {
						refill();
						index = 4;
						return buffer->as32[0];
					}
					uint32_t rv = *reinterpret_cast<uint32_t*>(&buffer->as8[index - 4]);
					return rv;
				}
				uint64_t MultiplexTransformRNG::raw64() {
					index += 15; index &= ~7;
					if (index > Tests::TestBlock::SIZE) {
						refill();
						index = 8;
						return buffer->as64[0];
					}
					uint64_t rv = *reinterpret_cast<uint64_t*>(&buffer->as8[index - 8]);
					return rv;
				}
				void MultiplexTransformRNG::seed(uint64_t seedval) {
					index = 999999;
					for (auto* vrng : source_rngs) {
							vrng->seed(seedval);
					}
				}
				void MultiplexTransformRNG::seed(vRNG* seeder) {
					index = 999999;
					for (auto* vrng : source_rngs) {
							vrng->seed(seeder);
					}
					/*static bool first = true;
					if (first) {
						std::printf("\n{\n");
						for (std::vector<vRNG*>::iterator it = source_rngs.begin(); it != source_rngs.end(); it++) std::printf("\t{%s:%s}\n", (*it)->get_name().c_str(), (*it)->print_state().c_str());
						std::printf("}\n");
						first = false;
					}*/
				}
				int MultiplexTransformRNG::get_native_output_size() const {
					int lowest = 9999, highest = -1;
					for (auto* source_rng : source_rngs) {
						int ls = source_rng->get_native_output_size();
						if (ls < lowest) lowest = ls;
						if (ls > highest) highest = ls;
					}
					if (lowest == highest) return lowest;
					return -1;
				}
				uint64_t MultiplexTransformRNG::get_flags() const {
					auto anded_bits = uint64_t(-1);
					auto ored_bits = uint64_t(0);
					for (auto* source_rng : source_rngs) {
						uint64_t lf = source_rng->get_flags();
						anded_bits &= lf;
						ored_bits |= lf;
					}
					using namespace PractRand::RNGs::FLAG;
					return (OUTPUT_IS_BUFFERED | STATE_UNAVAILABLE) |
						(anded_bits & (SUPPORTS_FASTFORWARD | SUPPORTS_ENTROPY_ACCUMULATION | USES_SPECIFIED | ENDIAN_SAFE)) |
						(ored_bits & (CRYPTOGRAPHIC_SECURITY | USES_MULTIPLICATION | USES_COMPLEX_INSTRUCTIONS | USES_VARIABLE_SHIFTS | USES_INDIRECTION | USES_CYCLIC_BUFFER | USES_FLOW_CONTROL | USES_BIT_SCANS | USES_OTHER_WORD_SIZES | OUTPUT_IS_HASHED));
				}
				void MultiplexTransformRNG::walk_state(StateWalkingObject* walker) {
					for (auto* source_rng : source_rngs) {
						source_rng->walk_state(walker);
					}
					if (!walker->is_read_only()) index = 999999;
				}
				MultiplexTransformRNG::MultiplexTransformRNG(const std::vector<vRNG*>& sources) : buffer(new Tests::TestBlock), source_rngs(sources) {
				}
				MultiplexTransformRNG::~MultiplexTransformRNG() { delete buffer; buffer = nullptr; }


				ReinterpretAsUnknown::ReinterpretAsUnknown( vRNG* rng ) : Transform8(rng) {
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as8[0];
				}
				ReinterpretAsUnknown::~ReinterpretAsUnknown() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					delete block;
				}
				void ReinterpretAsUnknown::refill() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					block->fill(base_rng);
					index = 0;
				}
				uint8_t ReinterpretAsUnknown::raw8() {
					if (index >= 8192 / OUTPUT_BITS) refill();
					return buffer[index++];
				}
				std::string ReinterpretAsUnknown::get_name() const {return std::string("AsUnknown(") + base_rng->get_name() + ")";}

				ReinterpretAs8::ReinterpretAs8( vRNG* rng ) : Transform8(rng) {
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as8[0];
				}
				ReinterpretAs8::~ReinterpretAs8() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					delete block;
				}
				void ReinterpretAs8::refill() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					block->fill(base_rng);
					index = 0;
				}
				uint8_t ReinterpretAs8::raw8() {
					if (index >= 8192 / OUTPUT_BITS) refill();
					return buffer[index++];
				}
				std::string ReinterpretAs8::get_name() const {return std::string("As8(") + base_rng->get_name() + ")";}

				ReinterpretAs16::ReinterpretAs16( vRNG* rng ) : Transform16(rng) {
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as16[0];
				}
				ReinterpretAs16::~ReinterpretAs16() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					delete block;
				}
				void ReinterpretAs16::refill() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					block->fill(base_rng);
					index = 0;
				}
				uint16_t ReinterpretAs16::raw16() {
					if (index >= 8192 / OUTPUT_BITS) refill();
					return buffer[index++];
				}
				std::string ReinterpretAs16::get_name() const {return std::string("As16(") + base_rng->get_name() + ")";}

				ReinterpretAs32::ReinterpretAs32( vRNG* rng ) : Transform32(rng) {
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as32[0];
				}
				ReinterpretAs32::~ReinterpretAs32() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					delete block;
				}
				void ReinterpretAs32::refill() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					block->fill(base_rng);
					index = 0;
				}
				uint32_t ReinterpretAs32::raw32() {
					if (index >= 8192 / OUTPUT_BITS) refill();
					return buffer[index++];
				}
				std::string ReinterpretAs32::get_name() const {return std::string("As32(") + base_rng->get_name() + ")";}

				ReinterpretAs64::ReinterpretAs64( vRNG* rng ) : Transform64(rng) {
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as64[0];
				}
				ReinterpretAs64::~ReinterpretAs64() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					delete block;
				}
				void ReinterpretAs64::refill() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer);
					block->fill(base_rng);
					index = 0;
				}
				uint64_t ReinterpretAs64::raw64() {
					if (index >= 8192 / OUTPUT_BITS) refill();
					return buffer[index++];
				}
				std::string ReinterpretAs64::get_name() const {return std::string("As64(") + base_rng->get_name() + ")";}

				void Xor::refill() {
					MultiplexTransformRNG::refill();
					PractRand::Tests::TestBlock tmp{};
					buffer->fill(source_rngs[0]);
					for (unsigned int sri = 1; sri < source_rngs.size(); sri++) {
						tmp.fill(source_rngs[sri]);
						for (int i = 0; i < PractRand::Tests::TestBlock::SIZE / 8; i++) buffer->as64[i] ^= tmp.as64[i];
					}
				}
				std::string Xor::get_name() const {
					std::string rv = "xor(";
					for (unsigned int sri = 0; sri < source_rngs.size(); sri++) {
						if (sri) rv += ",";
						rv += source_rngs[sri]->get_name();
					}
					rv += ")";
					return rv;
				}

				Discard16to8::Discard16to8(vRNG* base_rng_) : Transform8(base_rng_) {
					//if (base_rng_->get_native_output_size() != INPUT_BITS) std::cerr << "* warning: Discard16to8 using incorrect input size?\n";
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as16[0];
				}
				void Discard16to8::refill() { auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer); block->fill(base_rng); index = 0; }
				std::string Discard16to8::get_name() const { return std::string("Discard16to8(") + base_rng->get_name() + ")"; }
				uint8_t Discard16to8::raw8() {
					if (index >= 8192 / INPUT_BITS) refill();
					return OutWord(buffer[index++]);
				}
				Discard32to8::Discard32to8(vRNG* base_rng_) : Transform8(base_rng_) {
					//if (base_rng_->get_native_output_size() != INPUT_BITS) std::cerr << "* warning: Discard16to8 using incorrect input size?\n";
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as32[0];
				}
				std::string Discard32to8::get_name() const { return std::string("Discard32to8(") + base_rng->get_name() + ")"; }
				void Discard32to8::refill() { auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer); block->fill(base_rng); index = 0; }
				uint8_t Discard32to8::raw8() {
					if (index >= 8192 / INPUT_BITS) refill();
					return OutWord(buffer[index++]);
				}
				Discard64to8::Discard64to8(vRNG* base_rng_) : Transform8(base_rng_) {
					//if (base_rng_->get_native_output_size() != INPUT_BITS) std::cerr << "* warning: Discard16to8 using incorrect input size?\n";
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as64[0];
				}
				std::string Discard64to8::get_name() const { return std::string("Discard64to8(") + base_rng->get_name() + ")"; }
				void Discard64to8::refill() {
					auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer); block->fill(base_rng); index = 0;
				}
				uint8_t Discard64to8::raw8() {
					if (index >= 8192 / INPUT_BITS) refill();
					return OutWord(buffer[index++]);
				}
				Discard32to16::Discard32to16(vRNG* base_rng_) : Transform16(base_rng_) {
					//if (base_rng_->get_native_output_size() != INPUT_BITS) std::cerr << "* warning: Discard16to8 using incorrect input size?\n";
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as32[0];
				}
				std::string Discard32to16::get_name() const { return std::string("Discard32to16(") + base_rng->get_name() + ")"; }
				void Discard32to16::refill() { auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer); block->fill(base_rng); index = 0; }
				uint16_t Discard32to16::raw16() {
					if (index >= 8192 / INPUT_BITS) refill();
					return OutWord(buffer[index++]);
				}
				Discard64to16::Discard64to16(vRNG* base_rng_) : Transform16(base_rng_) {
					//if (base_rng_->get_native_output_size() != INPUT_BITS) std::cerr << "* warning: Discard16to8 using incorrect input size?\n";
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as64[0];
				}
				std::string Discard64to16::get_name() const { return std::string("Discard64to16(") + base_rng->get_name() + ")"; }
				void Discard64to16::refill() { auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer); block->fill(base_rng); index = 0; }
				uint16_t Discard64to16::raw16() {
					if (index >= 8192 / INPUT_BITS) refill();
					return OutWord(buffer[index++]);
				}
				Discard64to32::Discard64to32(vRNG* base_rng_) : Transform32(base_rng_) {
					//if (base_rng_->get_native_output_size() != INPUT_BITS) std::cerr << "* warning: Discard16to8 using incorrect input size?\n";
					auto* block = new PractRand::Tests::TestBlock;
					buffer = &block->as64[0];
				}
				std::string Discard64to32::get_name() const { return std::string("Discard64to32(") + base_rng->get_name() + ")"; }
				void Discard64to32::refill() { auto* block = reinterpret_cast<PractRand::Tests::TestBlock*>(buffer); block->fill(base_rng); index = 0; }
				uint32_t Discard64to32::raw32() {
					if (index >= 8192 / INPUT_BITS) refill();
					return OutWord(buffer[index++]);
				}

				void GeneralizedTableTransform::seed(uint64_t s) {base_rng->seed(s);}
				uint64_t GeneralizedTableTransform::get_flags() const {
					return base_rng->get_flags() | FLAG::USES_FLOW_CONTROL | FLAG::STATE_UNAVAILABLE;//not exactly, but close enough
				}
				GeneralizedTableTransform::~GeneralizedTableTransform() {delete base_rng;}
				void GeneralizedTableTransform::walk_state(StateWalkingObject* walker) {
					base_rng->walk_state(walker);
					buf_data = 0;
					buf_count = 0;
					finished_bytes.clear();
				}
				uint8_t GeneralizedTableTransform::raw8() {
					while (true) {
						if (!finished_bytes.empty()) {
							uint8_t rv = finished_bytes.front();
							finished_bytes.pop_front();
							return rv;
						}
						uint64_t in = base_rng->raw64();
						for (int i = 0; i < 8; i++) {
							const Entry& e = table[in & 255];
							in >>= 8;
							buf_data |= uint32_t(e.data) << buf_count;
							buf_count += e.count;
							if (buf_count >= 8) {
								finished_bytes.push_back(buf_data & 255);
								buf_count -= 8;
								buf_data >>= 8;
							}
						}
					}
				}
				static constexpr GeneralizedTableTransform::Entry self_shrinking_table11[256] = {
					{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },
					{ .data=0, .count=2 },{ .data=0, .count=2 },{ .data=0, .count=3 },{ .data=1, .count=3 },{ .data=1, .count=2 },{ .data=1, .count=2 },{ .data=2, .count=3 },{ .data=3, .count=3 },
					{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=2, .count=2 },{ .data=2, .count=2 },{ .data=4, .count=3 },{ .data=5, .count=3 },{ .data=3, .count=2 },{ .data=3, .count=2 },{ .data=6, .count=3 },{ .data=7, .count=3 },
					{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },{ .data=0, .count=0 },{ .data=0, .count=0 },{ .data=0, .count=1 },{ .data=1, .count=1 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },
					{ .data=0, .count=2 },{ .data=0, .count=2 },{ .data=0, .count=3 },{ .data=1, .count=3 },{ .data=1, .count=2 },{ .data=1, .count=2 },{ .data=2, .count=3 },{ .data=3, .count=3 },
					{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=2, .count=2 },{ .data=2, .count=2 },{ .data=4, .count=3 },{ .data=5, .count=3 },{ .data=3, .count=2 },{ .data=3, .count=2 },{ .data=6, .count=3 },{ .data=7, .count=3 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },
					{ .data=0, .count=2 },{ .data=0, .count=2 },{ .data=0, .count=3 },{ .data=1, .count=3 },{ .data=1, .count=2 },{ .data=1, .count=2 },{ .data=2, .count=3 },{ .data=3, .count=3 },
					{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },{ .data=0, .count=1 },{ .data=0, .count=1 },{ .data=0, .count=2 },{ .data=1, .count=2 },
					{ .data=0, .count=2 },{ .data=0, .count=2 },{ .data=0, .count=3 },{ .data=1, .count=3 },{ .data=1, .count=2 },{ .data=1, .count=2 },{ .data=2, .count=3 },{ .data=3, .count=3 },
					{ .data=0, .count=2 },{ .data=0, .count=2 },{ .data=0, .count=3 },{ .data=1, .count=3 },{ .data=0, .count=2 },{ .data=0, .count=2 },{ .data=0, .count=3 },{ .data=1, .count=3 },
					{ .data=0, .count=3 },{ .data=0, .count=3 },{ .data=0, .count=4 },{ .data=1, .count=4 },{ .data=1, .count=3 },{ .data=1, .count=3 },{ .data=2, .count=4 },{ .data=3, .count=4 },
					{ .data=1, .count=2 },{ .data=1, .count=2 },{ .data=2, .count=3 },{ .data=3, .count=3 },{ .data=1, .count=2 },{ .data=1, .count=2 },{ .data=2, .count=3 },{ .data=3, .count=3 },
					{ .data=2, .count=3 },{ .data=2, .count=3 },{ .data=4, .count=4 },{ .data=5, .count=4 },{ .data=3, .count=3 },{ .data=3, .count=3 },{ .data=6, .count=4 },{ .data=7, .count=4 },
					{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=2, .count=2 },{ .data=2, .count=2 },{ .data=4, .count=3 },{ .data=5, .count=3 },{ .data=3, .count=2 },{ .data=3, .count=2 },{ .data=6, .count=3 },{ .data=7, .count=3 },
					{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },{ .data=1, .count=1 },{ .data=1, .count=1 },{ .data=2, .count=2 },{ .data=3, .count=2 },
					{ .data=2, .count=2 },{ .data=2, .count=2 },{ .data=4, .count=3 },{ .data=5, .count=3 },{ .data=3, .count=2 },{ .data=3, .count=2 },{ .data=6, .count=3 },{ .data=7, .count=3 },
					{ .data=2, .count=2 },{ .data=2, .count=2 },{ .data=4, .count=3 },{ .data=5, .count=3 },{ .data=2, .count=2 },{ .data=2, .count=2 },{ .data=4, .count=3 },{ .data=5, .count=3 },
					{ .data=4, .count=3 },{ .data=4, .count=3 },{ .data=8, .count=4 },{ .data=9, .count=4 },{ .data=5, .count=3 },{ .data=5, .count=3 },{.data=10, .count=4 },{.data=11, .count=4 },
					{ .data=3, .count=2 },{ .data=3, .count=2 },{ .data=6, .count=3 },{ .data=7, .count=3 },{ .data=3, .count=2 },{ .data=3, .count=2 },{ .data=6, .count=3 },{ .data=7, .count=3 },
					{ .data=6, .count=3 },{ .data=6, .count=3 },{.data=12, .count=4 },{.data=13, .count=4 },{ .data=7, .count=3 },{ .data=7, .count=3 },{.data=14, .count=4 },{.data=15, .count=4 }
				};
				std::string GeneralizedTableTransform::get_name() const {return name;}
				vRNG* apply_SelfShrinkTransform(vRNG* base_rng) {
					return new GeneralizedTableTransform (base_rng, self_shrinking_table11, std::string("[SShrink]") + base_rng->get_name() );
				}
				//vRNG *apply_SimpleShrinkTransform(vRNG *base_rng) {
				//	return new GeneralizedTableTransform (base_rng, NULL, std::string("[Shrink3of4]") + base_rng->get_name() );
				//}



				uint64_t BaysDurhamShuffle64::raw64() {
					uint64_t& storage = table[prev];
					uint64_t rv = storage;
					storage = base_rng->raw64();
					prev = uint8_t(storage >> index_shift) & index_mask;
					return rv;
				}
				void BaysDurhamShuffle64::seed(uint64_t s) {
					base_rng->seed(s);
					for (int i = 0; i <= int{index_mask}; i++)
						table[i] = base_rng->raw64();
					for (int i = 0; i <= int{index_mask}; i++) {raw64();raw64();}
				}
				void BaysDurhamShuffle64::walk_state(StateWalkingObject* walker) {
					base_rng->walk_state(walker);
					if (!(walker->get_properties() & StateWalkingObject::FLAG_CLUMSY)) {
						walker->handle(index_mask);
						walker->handle(index_shift);
					}
					for (int i = 0; i <= int{index_mask}; i++) walker->handle(table[i]);
					walker->handle(prev);
					prev &= index_mask;
				}
				std::string BaysDurhamShuffle64::get_name() const {
					std::ostringstream tmp;
					tmp << "[BDS" << (1+int(index_mask)) << "]" << base_rng->get_name();
					return tmp.str();
				}

				uint32_t BaysDurhamShuffle32::raw32() {
					uint32_t& storage = table[prev];
					uint32_t rv = storage;
					storage = base_rng->raw32();
					prev = uint8_t(storage >> index_shift) & index_mask;
					return rv;
				}
				void BaysDurhamShuffle32::seed(uint64_t s) {
					base_rng->seed(s);
					for (int i = 0; i <= int{index_mask}; i++)
						table[i] = base_rng->raw32();
					for (int i = 0; i <= int{index_mask}; i++)  {raw32();raw32();}
				}
				void BaysDurhamShuffle32::walk_state(StateWalkingObject* walker) {
					base_rng->walk_state(walker);
					if (!(walker->get_properties() & StateWalkingObject::FLAG_CLUMSY)) {
						walker->handle(index_mask);
						walker->handle(index_shift);
					}
					for (int i = 0; i <= int{index_mask}; i++) walker->handle(table[i]);
					walker->handle(prev);
					prev &= index_mask;
				}
				std::string BaysDurhamShuffle32::get_name() const {
					std::ostringstream tmp;
					tmp << "[BDS" << (1+int(index_mask)) << "]" << base_rng->get_name();
					return tmp.str();
				}

				uint16_t BaysDurhamShuffle16::raw16() {
					uint16_t& storage = table[prev];
					uint16_t rv = storage;
					storage = base_rng->raw32();
					prev = uint8_t(storage >> index_shift) & index_mask;
					return rv;
				}
				void BaysDurhamShuffle16::seed(uint64_t s) {
					base_rng->seed(s);
					for (int i = 0; i <= int{index_mask}; i++)
						table[i] = base_rng->raw32();
					for (int i = 0; i <= int{index_mask}; i++)  {raw16();raw16();}
				}
				void BaysDurhamShuffle16::walk_state(StateWalkingObject* walker) {
					base_rng->walk_state(walker);
					if (!(walker->get_properties() & StateWalkingObject::FLAG_CLUMSY)) {
						walker->handle(index_mask);
						walker->handle(index_shift);
					}
					for (int i = 0; i <= int{index_mask}; i++) walker->handle(table[i]);
					walker->handle(prev);
					prev &= index_mask;
				}
				std::string BaysDurhamShuffle16::get_name() const {
					std::ostringstream tmp;
					tmp << "[BDS" << (1+int(index_mask)) << "]" << base_rng->get_name();
					return tmp.str();
				}

				uint8_t BaysDurhamShuffle8::raw8() {
					uint8_t& storage = table[prev];
					uint8_t rv = storage;
					storage = base_rng->raw32();
					prev = uint8_t(storage >> index_shift) & index_mask;
					return rv;
				}
				void BaysDurhamShuffle8::seed(uint64_t s) {
					base_rng->seed(s);
					for (int i = 0; i <= int{index_mask}; i++)
						table[i] = base_rng->raw32();
					for (int i = 0; i <= int{index_mask}; i++)  {raw8();raw8();}
				}
				void BaysDurhamShuffle8::walk_state(StateWalkingObject* walker) {
					base_rng->walk_state(walker);
					if (!(walker->get_properties() & StateWalkingObject::FLAG_CLUMSY)) {
						walker->handle(index_mask);
						walker->handle(index_shift);
					}
					for (int i = 0; i <= int{index_mask}; i++) walker->handle(table[i]);
					walker->handle(prev);
					prev &= index_mask;
				}
				std::string BaysDurhamShuffle8::get_name() const {
					std::ostringstream tmp;
					tmp << "[BDS" << (1+int(index_mask)) << "]" << base_rng->get_name();
					return tmp.str();
				}
				vRNG* apply_BaysDurhamShuffle(vRNG* base_rng, int table_size_L2, int shift) {
					auto* tmp8 = dynamic_cast<vRNG8*>(base_rng);
					if (tmp8)  return new BaysDurhamShuffle8 (tmp8, table_size_L2, shift >= 0 ? shift : 8-table_size_L2);
					auto* tmp16 = dynamic_cast<vRNG16*>(base_rng);
					if (tmp16) return new BaysDurhamShuffle16(tmp16, table_size_L2, shift >= 0 ? shift : 16-table_size_L2);
					auto* tmp32 = dynamic_cast<vRNG32*>(base_rng);
					if (tmp32) return new BaysDurhamShuffle32(tmp32, table_size_L2, shift >= 0 ? shift : 32-table_size_L2);
					auto* tmp64 = dynamic_cast<vRNG64*>(base_rng);
					if (tmp64) return new BaysDurhamShuffle64(tmp64, table_size_L2, shift >= 0 ? shift : 64-table_size_L2);
					issue_error();
					return nullptr;//just to quiet the warnings
				}
}
