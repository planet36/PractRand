#include <string>
#include <sstream>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>

#include "PractRand/RNGs/other/lcgish.h"
//#include "PractRand/test_helpers.h"

namespace PractRand {
	using namespace Internals;
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				Uint32 lcg32of64_varqual::raw32() {
					state = state * 1103515245 + 12345;
					return Uint32(state >> outshift);
				}
				std::string lcg32of64_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void lcg32of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint16 lcg16of64_varqual::raw16() {
					state = state * 1103515245 + 12345;
					return Uint16(state >> outshift);
				}
				std::string lcg16of64_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void lcg16of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint8 lcg8of64_varqual::raw8() {
					state = state * 1103515245 + 12345;
					return Uint8(state >> outshift);
				}
				std::string lcg8of64_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void lcg8of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}

				Uint32 lcg32of128_varqual::raw32() {
					//large multiplication is much harder to write in C (than asm)
					//made some compromises here... restricting the multiplier to 32 bits
					//which would hurt quality some, being so small compared to the state
					//but I think the correct comparison is to the output window, not the state
					//which is only 16 bits, so it should all be good
					const Uint32 multiplier = 1103515245;
					const Uint64 adder = 1234567;
					Uint64 a = Uint32(low) * Uint64(multiplier);
					Uint64 b = (low >> 32) * Uint64(multiplier);
					Uint64 old = low;
					low = a + (b << 32);
					b += a >> 32;
					high *= multiplier;
					high += b >> 32;
					high += old;//adds 2**64 to the multiplier
					low += adder;
					if (a < adder) high++;
					if (outshift >= 64) return Uint32(high >> (outshift-64));
					if (outshift > 32) return Uint32( (low >> outshift) | (high << (64-outshift)) );
					return Uint32(low >> outshift);
				}
				std::string lcg32of128_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void lcg32of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
				}
				Uint16 lcg16of128_varqual::raw16() {
					const Uint32 multiplier = 1103515245;
					const Uint64 adder = 1234567;
					Uint64 a = Uint32(low) * Uint64(multiplier);
					Uint64 b = (low >> 32) * Uint64(multiplier);
					Uint64 old = low;
					low = a + (b << 32);
					b += a >> 32;
					high *= multiplier;
					high += b >> 32;
					high += old;//adds 2**64 to the multiplier
					low += adder;
					if (a < adder) high++;
					if (outshift >= 64) return Uint16(high >> (outshift-64));
					if (outshift > 48) return Uint16( (low >> outshift) | (high << (64-outshift)) );
					return Uint16(low >> outshift);
				}
				std::string lcg16of128_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void lcg16of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
				}
				Uint8 lcg8of128_varqual::raw8() {
					const Uint32 multiplier = 1103515245;
					const Uint64 adder = 1234567;
					Uint64 a = Uint32(low) * Uint64(multiplier);
					Uint64 b = (low >> 32) * Uint64(multiplier);
					Uint64 old = low;
					low = a + (b << 32);
					b += a >> 32;
					high *= multiplier;
					high += b >> 32;
					high += old;//adds 2**64 to the multiplier
					low += adder;
					if (a < adder) high++;
					if (outshift >= 64) return Uint8(high >> (outshift-64));
					if (outshift > 56) return Uint8( (low >> outshift) | (high << (64-outshift)) );
					return Uint8(low >> outshift);
				}
				std::string lcg8of128_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void lcg8of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
				}

				// non-power-of-2 LCGs

				//helpers are provided for sizes: 10, 11, 13, 17, 19, 22, 24, 26, 29, 31, 33, 36, 39, 44 (in bits)
				// basic constants to make more include:
				// 2^X-K is prime, and works well with multiplier M which is also small enough that overflowing 64 bit values isn't an issue, but reasonably large with the range of numbers meeting those criteria
				// X=			45			46			47			48			49			50			51			52			53			54			55			56			57
				// K=			55			21			115			59			81			27			129			47			111			33			55			5*			13
				//in binary		110111		10101		1110011		111011		1010001		11011		10000001	101111		1101111		100001		110111		101			1101
				// M=			?

				Uint16 np2clcg2_32::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_32::get_name() const { return "np2clcg2_32"; }
				void np2clcg2_32::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_36::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_36::get_name() const { return "np2clcg2_36"; }
				void np2clcg2_36::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_39::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_39::get_name() const { return "np2clcg2_39"; }
				void np2clcg2_39::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_45::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_45::get_name() const { return "np2clcg2_45"; }
				void np2clcg2_45::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_53::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_53::get_name() const { return "np2clcg2_53"; }
				void np2clcg2_53::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_57::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_57::get_name() const { return "np2clcg2_57"; }
				void np2clcg2_57::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_60::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_60::get_name() const { return "np2clcg2_60"; }
				void np2clcg2_60::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_64::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_64::get_name() const { return "np2clcg2_64"; }
				void np2clcg2_64::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_67::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_67::get_name() const { return "np2clcg2_67"; }
				void np2clcg2_67::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}
				Uint16 np2clcg2_70::raw16() {
					lcg_A.advance(); lcg_B.advance();
					return lcg_A.state + lcg_B.state;
				}
				std::string np2clcg2_70::get_name() const { return "np2clcg2_70"; }
				void np2clcg2_70::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state);
					lcg_A.sanitize(); lcg_B.sanitize();
				}

				Uint16 np2clcg3_40::raw16() {
					lcg_A.advance(); lcg_B.advance(); lcg_C.advance();
					return lcg_A.state + lcg_B.state + lcg_C.state;
				}
				std::string np2clcg3_40::get_name() const { return "np2clcg3_40"; }
				void np2clcg3_40::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state); walker->handle(lcg_C.state);
					lcg_A.sanitize(); lcg_B.sanitize(); lcg_C.sanitize();
				}
				Uint16 np2clcg3_49::raw16() {
					lcg_A.advance(); lcg_B.advance(); lcg_C.advance();
					return lcg_A.state + lcg_B.state + lcg_C.state;
				}
				std::string np2clcg3_49::get_name() const { return "np2clcg3_49"; }
				void np2clcg3_49::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state); walker->handle(lcg_C.state);
					lcg_A.sanitize(); lcg_B.sanitize(); lcg_C.sanitize();
				}
				Uint16 np2clcg3_58::raw16() {
					lcg_A.advance(); lcg_B.advance(); lcg_C.advance();
					return lcg_A.state + lcg_B.state + lcg_C.state;
				}
				std::string np2clcg3_58::get_name() const { return "np2clcg3_58"; }
				void np2clcg3_58::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state); walker->handle(lcg_C.state);
					lcg_A.sanitize(); lcg_B.sanitize(); lcg_C.sanitize();
				}
				Uint16 np2clcg3_67::raw16() {
					lcg_A.advance(); lcg_B.advance(); lcg_C.advance();
					return lcg_A.state + lcg_B.state + lcg_C.state;
				}
				std::string np2clcg3_67::get_name() const { return "np2clcg3_67"; }
				void np2clcg3_67::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg_A.state); walker->handle(lcg_B.state); walker->handle(lcg_C.state);
					lcg_A.sanitize(); lcg_B.sanitize(); lcg_C.sanitize();
				}

				np2lcg32_varqual::np2lcg32_varqual(Uint32 _m, Uint32 _a, Uint32 _c) {
					m = _m;
					if (!m) issue_error("np2lcg32_varqual: LCG modulus is zero");
					if (_a >= m) issue_error("np2lcg32_varqual: multiplier out of range");
					if (_c >= m) issue_error("np2lcg32_varqual: additive out of range");
					a = _a % m;
					c = _c % m;
				}
				Uint32 np2lcg32_varqual::raw32() {
					Uint64 x2 = x;
					while (true) {
						x2 *= a;
						x2 += c;
						x2 %= m;
						x = x2;
						return x;
					}
				}
				std::string np2lcg32_varqual::get_name() const {
					std::ostringstream str;
					str << "np2lcg(16," << m << "," << a << "," << c << ")";
					return str.str();
				}
				void np2lcg32_varqual::walk_state(StateWalkingObject *s) {
					s->handle(x);
					if (x >= m) x %= m;
				}

				np2lcg16_varqual::np2lcg16_varqual(Uint32 _m, Uint32 _a, Uint32 _c) {
					m = _m;
					if (!m) issue_error("np2lcg16_varqual: LCG modulus is zero");
					//if (m < 65536) issue_error("np2lcg16_varqual: LCG modulus is less than 65536");
					if (_a >= m) issue_error("np2lcg16_varqual: multiplier out of range");
					if (_c >= m) issue_error("np2lcg16_varqual: additive out of range");
					//mm = m & 0xFFFF0000;
					a = _a % m;
					c = _c % m;
				}
				Uint16 np2lcg16_varqual::raw16() {
					Uint64 x2 = x;
					while (true) {
						x2 *= a;
						x2 += c;
						x2 %= m;
						x = x2;
						//if (x < mm) return x;
						return x;
					}
				}
				std::string np2lcg16_varqual::get_name() const {
					std::ostringstream str;
					str << "np2lcg(16," << m << "," << a << "," << c << ")";
					return str.str();
				}
				void np2lcg16_varqual::walk_state(StateWalkingObject *s) {
					s->handle(x);
					if (x >= m) x %= m;
				}
				np2lcg8_varqual::np2lcg8_varqual(Uint32 _m, Uint32 _a, Uint32 _c) {
					m = _m;
					if (!m) issue_error("np2lcg8_varqual: LCG modulus is zero");
					//if (m < 256) issue_error("np2lcg8_varqual: LCG modulus is less than 256");
					if (_a >= m) issue_error("np2lcg8_varqual: multiplier out of range");
					if (_c >= m) issue_error("np2lcg8_varqual: additive out of range");
					//mm = m & 0xFFFFFF00;
					a = _a % m;
					c = _c % m;
				}
				Uint8 np2lcg8_varqual::raw8() {
					Uint64 x2 = x;
					while (true) {
						x2 *= a;
						x2 += c;
						x2 %= m;
						x = x2;
						//if (x < mm) return x;
						return x;
					}
				}
				std::string np2lcg8_varqual::get_name() const {
					std::ostringstream str;
					str << "np2lcg(8," << m << "," << a << "," << c << ")";
					return str.str();
				}
				void np2lcg8_varqual::walk_state(StateWalkingObject *s) {
					s->handle(x);
					if (x >= m) x %= m;
				}


				//similar to lcg16of32, but with a longer period
				Uint16 lcg16of32_extended::raw16() {
					state = state * 1103515245 + add;
					if (!(state & 0x7fff)) add += 2;
					return Uint16(state >> 16);
				}
				std::string lcg16of32_extended::get_name() const {return "lcg16of32_extended";}
				void lcg16of32_extended::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(add);
					add |= 1;
				}

				//similar to lcg32, but with a longer period
				Uint32 lcg32_extended::raw32() {
					state = state * 1103515245 + add;
					if (!(state & 0x7fff)) add += 2;
					return state;
				}
				std::string lcg32_extended::get_name() const {return "lcg32_extended";}
				void lcg32_extended::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(add);
					add |= 1;
				}
				Uint32 clcg32of95_varqual::raw32() {
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = Uint64(lcg2) * 1579544716;
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return lcg2 + Uint32(lcg1 >> outshift);
				}
				std::string clcg32of95_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(32," << (32 + outshift + 31) << ")";
					return str.str();
				}
				void clcg32of95_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					lcg2 &= 0x7FffFFff;
					if (!lcg2) lcg2 = 1;
				}
				Uint16 clcg16of95_varqual::raw16() {
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = Uint64(lcg2) * 1579544716;
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					//return Uint16(lcg2 >> 12) + Uint16(lcg1 >> outshift);
					return Uint16(lcg2) + Uint16(lcg1 >> outshift);
					//return Uint16(((Uint64(lcg2) << ) + lcg1) >> outshift);
				}
				std::string clcg16of95_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(16," << (16 + outshift + 31) << ")";
					return str.str();
				}
				void clcg16of95_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					lcg2 &= 0x7FffFFff;
					if (!lcg2) lcg2 = 1;
				}
				Uint8 clcg8of95_varqual::raw8() {
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint8(lcg2) + Uint8(lcg1 >> outshift);
				}
				std::string clcg8of95_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(8," << (8 + outshift + 31) << ")";
					return str.str();
				}
				void clcg8of95_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					lcg2 &= 0x7FffFFff;
					if (!lcg2) lcg2 = 1;
				}
				Uint32 clcg32of108_varqual::raw32() {
					p2mlcg *= 0x5ED963A741C64E6Dull;
					p2mlcg += 12345;
					np2lcg.advance();
					return Uint32(np2lcg.state) + Uint32(p2mlcg >> outshift);
				}
				std::string clcg32of108_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(32," << (32 + outshift + NP2LCG_BITS) << ")";
					return str.str();
				}
				void clcg32of108_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(p2mlcg);
					walker->handle(np2lcg.state);
					np2lcg.sanitize();
				}
				Uint16 clcg16of108_varqual::raw16() {
					p2mlcg *= 0x5ED963A741C64E6Dull;
					p2mlcg += 12345;
					np2lcg.advance();
					return Uint16(np2lcg.state) + Uint16(p2mlcg >> outshift);
				}
				std::string clcg16of108_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(16," << (16 + outshift + NP2LCG_BITS) << ")";
					return str.str();
				}
				void clcg16of108_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(p2mlcg);
					walker->handle(np2lcg.state);
					np2lcg.sanitize();
				}
				Uint8 clcg8of108_varqual::raw8() {
					/*
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = Uint64(lcg2) * 1579544716;
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					//return Uint16(lcg2 >> 12) + Uint16(lcg1 >> outshift);
					return Uint16(lcg2) + Uint16(lcg1 >> outshift);
					//return Uint16(((Uint64(lcg2) << ) + lcg1) >> outshift);
					*/
					p2mlcg *= 0x5ED963A741C64E6Dull;
					p2mlcg += 12345;
					np2lcg.advance();
					return Uint8(np2lcg.state) + Uint8(p2mlcg >> outshift);
				}
				std::string clcg8of108_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(8," << (8 + outshift + NP2LCG_BITS) << ")";
					return str.str();
				}
				void clcg8of108_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(p2mlcg);
					walker->handle(np2lcg.state);
					np2lcg.sanitize();
				}

				void pcg32::seed(Uint64 s) { state = 0; raw32(); state += s; raw32(); }
				Uint32 pcg32::raw32() {
					Uint64 oldstate = state;
					state = state * 0x5851f42d4c957f2dULL + inc;
					Uint32 xorshifted = Uint32(((oldstate >> 18u) ^ oldstate) >> 27u);
					Uint32 rot = Uint32(oldstate >> 59u);
					return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
				}
				std::string pcg32::get_name() const {
					if (inc == 0xda3e39cb94b95bdbULL) return "pcg32";
					std::ostringstream str;
					str << "pcg32(" << std::hex << inc << ")";
					return str.str();
				}
				void pcg32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				void pcg32_norot::seed(Uint64 s) { state = 0; raw32(); state += s; raw32(); }
				Uint32 pcg32_norot::raw32() {
					Uint64 oldstate = state;
					state = state * 0x5851f42d4c957f2dULL + inc;
					Uint32 xorshifted = Uint32(((oldstate >> 18u) ^ oldstate) >> 27u);
					return xorshifted;
					//Uint32 rot = Uint32(oldstate >> 59u);
					//return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
				}
				std::string pcg32_norot::get_name() const {
					if (inc == 0xda3e39cb94b95bdbULL) return "pcg32_norot";
					std::ostringstream str;
					str << "pcg32_norot(" << std::hex << inc << ")";
					return str.str();
				}
				void pcg32_norot::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				void cmrg32of192::seed(Uint64 s) {
					n1m0 = Uint32(s);
					n1m1 = n1m2 = 1;
					n2m0 = Uint32(s >> 32);
					n2m1 = n2m2 = 1;
				}
				Uint32 cmrg32of192::raw32() {
					Uint64 n1 = (Uint64(n1m1) * 1403580 - Uint64(n1m2) * 810728) % ((1ull << 32) - 209);
					Uint64 n2 = (Uint64(n2m0) * 527612 - Uint64(n2m2) * 1370589) % ((1ull << 32) - 22853);
					n1m2 = n1m1;
					n1m1 = n1m0;
					n1m0 = Uint32(n1);
					n2m2 = n2m1;
					n2m1 = n2m0;
					n2m0 = Uint32(n2);
					return Uint32(n1 + n2);
				}
				std::string cmrg32of192::get_name() const {return "cmrg32of192";}
				void cmrg32of192::walk_state(StateWalkingObject *walker) {
					walker->handle(n1m0);
					walker->handle(n1m1);
					walker->handle(n1m2);
					walker->handle(n2m0);
					walker->handle(n2m1);
					walker->handle(n2m2);
				}

				Uint32 xsh_lcg_bad::raw32() {
					Uint64 tmp = x1 ^ (x1 << 11);
					x1 = x2;
					x2 = x3;
					x3 = x0;
					x0 = (x0 >> 19) ^ tmp ^ (tmp >> 8);
					lcg = (lcg * 279470273) % 4294967291;
					return Uint32(x0 ^ lcg);
				}
				std::string xsh_lcg_bad::get_name() const { return "xsh_lcg_bad"; }
				void xsh_lcg_bad::seed(Uint64 s) {
					x1 = s;
					x0 = x2 = x3 = 0xFFffFFffFFffFFffull;//changed to prevent the bad all-zeroes case
					lcg = 2233445566;
					for (int i = 0; i < 64; i++) raw32();
				}
				void xsh_lcg_bad::walk_state(StateWalkingObject *walker) {
					walker->handle(x0);
					walker->handle(x1);
					walker->handle(x2);
					walker->handle(x3);
					walker->handle(lcg);
				}


				Uint32 xlcg32of64_varqual::raw32() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					state = (state ^ X) * M;
					return Uint32(state >> outshift);
				}
				std::string xlcg32of64_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void xlcg32of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint16 xlcg16of64_varqual::raw16() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					state = (state ^ X) * M;
					return Uint16(state >> outshift);
				}
				std::string xlcg16of64_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void xlcg16of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint8 xlcg8of64_varqual::raw8() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					state = (state ^ X) * M;
					return Uint8(state >> outshift);
				}
				std::string xlcg8of64_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void xlcg8of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 xlcg32of128_varqual::raw32() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					Uint64 a = Uint32(low) * Uint64(M);
					Uint64 b = (low >> 32) * Uint64(M);
					low = a + (b << 32);
					b += a >> 32;
					high = high * M + (b >> 32);
					low ^= X;
					if (outshift >= 64) return Uint16(high >> (outshift - 64));
					if (outshift > 48) return Uint16((low >> outshift) | (high << (64 - outshift)));
					return Uint16(low >> outshift);
				}
				std::string xlcg32of128_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void xlcg32of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
				}
				Uint16 xlcg16of128_varqual::raw16() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					Uint64 a = Uint32(low) * Uint64(M);
					Uint64 b = (low >> 32) * Uint64(M);
					low = a + (b << 32);
					b += a >> 32;
					high = high * M + (b >> 32);
					low ^= X;
					if (outshift >= 64) return Uint16(high >> (outshift - 64));
					if (outshift > 48) return Uint16((low >> outshift) | (high << (64 - outshift)));
					return Uint16(low >> outshift);
				}
				std::string xlcg16of128_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void xlcg16of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
				}
				Uint8 xlcg8of128_varqual::raw8() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					Uint64 a = Uint32(low) * Uint64(M);
					Uint64 b = (low >> 32) * Uint64(M);
					low = a + (b << 32);
					b += a >> 32;
					high = high * M + (b >> 32);
					low ^= X;
					if (outshift >= 64) return Uint8(high >> (outshift - 64));
					if (outshift > 56) return Uint8((low >> outshift) | (high << (64 - outshift)));
					return Uint8(low >> outshift);
				}
				std::string xlcg8of128_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void xlcg8of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
				}

				Uint32 cxlcg32of96_varqual::raw32() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					lcg1 = (lcg1 ^ X) * M;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return lcg2 + Uint32(lcg1 >> outshift);
				}
				std::string cxlcg32of96_varqual::get_name() const {
					std::ostringstream str;
					str << "cxlcg(32," << (32 + outshift + 31) << ")";
					return str.str();
				}
				void cxlcg32of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}
				Uint16 cxlcg16of96_varqual::raw16() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					lcg1 = (lcg1 ^ X) * M;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint16(lcg2) + Uint16(lcg1 >> outshift);
				}
				std::string cxlcg16of96_varqual::get_name() const {
					std::ostringstream str;
					str << "cxlcg(16," << (16 + outshift + 31) << ")";
					return str.str();
				}
				void cxlcg16of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}
				Uint8 cxlcg8of96_varqual::raw8() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					lcg1 = (lcg1 ^ X) * M;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint8(lcg2) + Uint8(lcg1 >> outshift);
				}
				std::string cxlcg8of96_varqual::get_name() const {
					std::ostringstream str;
					str << "cxlcg(8," << (8 + outshift + 31) << ")";
					return str.str();
				}
				void cxlcg8of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}

				bigbadlcg64X::bigbadlcg64X(int discard_bits_, int shift_bits_) : discard_bits(discard_bits_) {
					int max_discard_bits = MAX_N * 64 - 64;
					if (discard_bits_ < 0 || discard_bits_ > max_discard_bits) issue_error("bigbadlcg64 - discard_bits out of range (0 <= discard_bits <= ?960?)");
					n = (discard_bits_ + 63) / 64 + 1;
					shift_i = shift_bits_ / 64;
					shift_b = shift_bits_ & 63;
					if (shift_bits_ < 1 || shift_bits_ > discard_bits_) issue_error("bigbadlcg64 - shift_bits out of range (1 <= shift_bits <= discard_bits)");
				}
				//bigbadlcgX::~bigbadlcgX() { delete[] state; }
				Uint64 bigbadlcg64X::raw64() {
					static const Uint64 K = 0xA3EC647659359ACDull;
					Uint64 rv = state[n - 1];
					if (discard_bits & 63) {
						int b = discard_bits & 63;
						rv = (rv << (64 - b)) | (state[n - 2] >> b);
					}
					Uint64 olda[MAX_N];
					for (int i = 0; i < n; i++) olda[i] = state[i];
					bool carry = false;
					if (shift_b) {
						if (true) {
							Uint64 old = olda[0] << shift_b;
							state[shift_i] += old;
							carry = state[shift_i] < old;
						}
						for (int i = shift_i + 1; i < n; i++) {
							Uint64 old = (olda[i - shift_i] << shift_b) | (olda[i - shift_i - 1] >> (64 - shift_b));
							state[i] += old;
							bool c1 = state[i] < old;
							state[i] += carry ? 1 : 0;
							bool c2 = carry && !state[i];
							carry = c1 || c2;
						}
					}
					else {
						if (true) {
							Uint64 old = olda[0];
							state[shift_i] += old;
							bool carry = state[shift_i] < old;
						}
						for (int i = shift_i + 1; i < n; i++) {
							Uint64 old = olda[i - shift_i];
							state[i] += old;
							bool c1 = state[i] < old;
							state[i] += carry ? 1 : 0;
							bool c2 = carry && !state[i];
							carry = c1 || c2;
						}
					}
					state[0] += K;
					state[1] += (state[0] < K) ? 1 : 0;
					if (state[1] == 0 && state[0] < K) {
						for (int i = 2; i < n && !state[i - 1]; i++) state[i]++;
					}
					return rv;
				}
				std::string bigbadlcg64X::get_name() const {
					std::ostringstream str;
					str << "bblcg(64," << (discard_bits + OUTPUT_BITS) << "," << (shift_i * 64 + shift_b) << ")";
					return str.str();
				}
				void bigbadlcg64X::walk_state(StateWalkingObject *walker) {
					for (int i = 0; i < n; i++) walker->handle(state[i]);
				}
				bigbadlcg32X::bigbadlcg32X(int discard_bits_, int shift_) : base_lcg(discard_bits_, shift_) {}
				Uint32 bigbadlcg32X::raw32() { return Uint32(base_lcg.raw32()); }
				std::string bigbadlcg32X::get_name() const {
					std::ostringstream str;
					str << "bblcg(32," << (base_lcg.discard_bits + OUTPUT_BITS) << "," << (base_lcg.shift_i * 64 + base_lcg.shift_b) << ")";
					return str.str();
				}
				void bigbadlcg32X::walk_state(StateWalkingObject *walker) { base_lcg.walk_state(walker); }
				bigbadlcg16X::bigbadlcg16X(int discard_bits_, int shift_) : base_lcg(discard_bits_, shift_) {}
				Uint16 bigbadlcg16X::raw16() { return Uint16(base_lcg.raw64()); }
				std::string bigbadlcg16X::get_name() const {
					std::ostringstream str;
					str << "bblcg(16," << (base_lcg.discard_bits + OUTPUT_BITS) << "," << (base_lcg.shift_i * 64 + base_lcg.shift_b) << ")";
					return str.str();
				}
				void bigbadlcg16X::walk_state(StateWalkingObject *walker) { base_lcg.walk_state(walker); }
				bigbadlcg8X::bigbadlcg8X(int discard_bits_, int shift_) : base_lcg(discard_bits_, shift_) {}
				Uint8 bigbadlcg8X::raw8() { return Uint8(base_lcg.raw32()); }
				std::string bigbadlcg8X::get_name() const {
					std::ostringstream str;
					str << "bblcg(8," << (base_lcg.discard_bits + OUTPUT_BITS) << "," << (base_lcg.shift_i * 64 + base_lcg.shift_b) << ")";
					return str.str();
				}
				void bigbadlcg8X::walk_state(StateWalkingObject *walker) { base_lcg.walk_state(walker); }

			}
		}
	}
}
