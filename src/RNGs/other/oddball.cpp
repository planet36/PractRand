#include <string>
#include <sstream>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>

#include "PractRand/RNGs/other/oddball.h"
//#include "PractRand/test_helpers.h"

namespace PractRand {
	using namespace Internals;
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				Uint16 mmr16::raw16() {
					Uint16 old = a;
					a = b * 0x69ad;
					b = rotate16(b, 7) ^ c;
					c = rotate16(c, 5) + old;
					return old;
				}
				std::string mmr16::get_name() const { return "mmr16"; }
				void mmr16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 mmr32::raw32() {
					Uint32 old = a;
					a = b * 0xAC4969AD;
					b = rotate16(b, 13) ^ c;
					c = rotate16(c, 9) + old;
					return old;
				}
				std::string mmr32::get_name() const { return "mmr32"; }
				void mmr32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint16 garthy16::raw16() {
					if (!counter) scale += 2;
					scale += 2;
					Uint16 temp = value * scale;
					value += ((temp << 7) | (temp >> 9)) ^ counter++;
					return value;
				}
				std::string garthy16::get_name() const { return "garthy16"; }
				void garthy16::walk_state(StateWalkingObject *walker) {
					walker->handle(value); walker->handle(counter); walker->handle(scale);
					scale |= 1;
				}
				Uint32 garthy32::raw32() {
					if (!counter) scale += 2;
					scale += 2;
					Uint32 temp = value * scale;
					value += ((temp << 13) | (temp >> 19)) ^ counter++;
					return value;
				}
				std::string garthy32::get_name() const { return "garthy32"; }
				void garthy32::walk_state(StateWalkingObject *walker) {
					walker->handle(value); walker->handle(counter); walker->handle(scale);
					scale |= 1;
				}

				Uint16 binarymult16::raw16() {
					//with Bays-Durham shuffle (size 16) fails @ 32 GB
					Uint16 old = a;
					a = b * (c | 1);
					b = c ^ (old >> 7);
					c ^= old + d++;
					return a;
				}
				std::string binarymult16::get_name() const { return "binarymult16"; }
				void binarymult16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}
				Uint32 binarymult32::raw32() {
					Uint32 old = a;
					a = b * (c | 1);
					b = c ^ (old >> 13);
					c ^= old + d++;
					return a;
				}
				std::string binarymult32::get_name() const { return "binarymult32"; }
				void binarymult32::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}

				Uint16 rxmult16::raw16() {
					if (!a) { c++; if (!c) { c = 1; d += 2; } }
					a = a * 0x9ad + d;
					b = (((b << 7) | (b >> 9)) + a) ^ c;
					Uint16 tmp = b * 5245;
					tmp ^= tmp >> 8;
					return tmp + a;
				}
				std::string rxmult16::get_name() const { return "rxmult16"; }
				void rxmult16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
					d |= 1;
				}

				Uint64 multish2x64::raw64() {
					Uint64 old = ~a;
					a = (a * 0xa536c4b9) + b;
					b += (old << 21) | (old >> 43);
					return old;
				}
				std::string multish2x64::get_name() const { return "multish2x64"; }
				void multish2x64::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b);
				}
				Uint32 multish3x32::raw32() {
					Uint32 old = a;
					a = (b * 0xa536c4b9) + c++;
					b = ((b << 7) | (b >> 25)) + old;
					return old;
				}
				std::string multish3x32::get_name() const { return "multish3x32"; }
				void multish3x32::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 multish4x16::raw16() {
					Uint16 old = a;
					if (!c++) d++;
					a = (b^d) * 0x96b9 + c;
					b = ((b << 5) | (b >> 11)) ^ old;
					return old;
				}
				std::string multish4x16::get_name() const { return "multish4x16"; }
				void multish4x16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}

				Uint16 mwrca16::raw16() {
					Uint16 old = a * 0xa395;
					a ^= rotate16(b, 8);
					b = old;
					return a + b;
				}
				std::string mwrca16::get_name() const { return "mwrca16"; }
				void mwrca16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}
				Uint32 mwrca32::raw32() {
					Uint32 old = a * 0x5925a395;
					a ^= rotate32(b, 16);
					b = old;
					return a + b;
				}
				std::string mwrca32::get_name() const { return "mwrca32"; }
				void mwrca32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}
				Uint64 mwrca64::raw64() {
					Uint64 old = a * 0x5925a395864b139dull;
					a ^= rotate64(b, 32);
					b = old;
					return a + b;
				}
				std::string mwrca64::get_name() const { return "mwrca64"; }
				void mwrca64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}

				Uint16 mwrcc16::raw16() {
					Uint16 old = a * 0xa395;
					a += counter++;
					a ^= rotate16(b, 8);
					b = old;
					return a;
				}
				std::string mwrcc16::get_name() const { return "mwrcc16"; }
				void mwrcc16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint32 mwrcc32::raw32() {
					Uint32 old = a * 0x5925a395;
					a += counter++;
					a ^= rotate32(b, 16);
					b = old;
					return a;
				}
				std::string mwrcc32::get_name() const { return "mwrcc32"; }
				void mwrcc32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}

				Uint16 mwrcca16::raw16() {
					Uint16 old = a * 0xa395;
					a += counter++;
					a ^= rotate16(b, 8);
					b = old;
					return a + b;
				}
				std::string mwrcca16::get_name() const { return "mwrcca16"; }
				void mwrcca16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint32 mwrcca32::raw32() {
					Uint32 old = a * 0x5925a395;
					a += counter++;
					a ^= rotate32(b, 16);
					b = old;
					return a + b;
				}
				std::string mwrcca32::get_name() const { return "mwrcca32"; }
				void mwrcca32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint64 mwrcca64::raw64() {
					Uint64 old = a * 0x5925a395864b139dull;
					a += counter++;
					a ^= rotate64(b, 32);
					b = old;
					return a + b;
				}
				std::string mwrcca64::get_name() const { return "mwrcca64"; }
				void mwrcca64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}

				Uint16 old_mwlac16::raw16() {
					Uint16 oa;
					oa = a;
					a = (b * 0x9785) ^ (a >> 7);
					b = c + (oa >> 2);
					c = d;
					d += ~oa;
					return c;
				}
				std::string old_mwlac16::get_name() const { return "old_mwlac16"; }
				void old_mwlac16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}
				Uint16 mwlac_varA::raw16() {
					Uint16 oa;
					oa = a * 0x9785;//   1001011110000101
					a = b ^ rotate16(a, 7);
					b += c;
					c = oa;
					return c;
				}
				std::string mwlac_varA::get_name() const { return "mwlac_varA"; }
				void mwlac_varA::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varB::raw16() {
					Uint16 oa;
					oa = a * 0x9785;//   1001011110000101
					b = rotate(b, 13);
					a = b ^ rotate16(a, 7);
					b += c;
					c = oa;
					return c;
				}
				std::string mwlac_varB::get_name() const { return "mwlac_varB"; }
				void mwlac_varB::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varC::raw16() {
					a *= 0x9785;//   1001011110000101
					b = rotate16(b, 5);
					c = rotate16(c, 13);
					b ^= a;
					a ^= c;
					c += b;
					return b;
				}
				std::string mwlac_varC::get_name() const { return "mwlac_varC"; }
				void mwlac_varC::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varD::raw16() {
					a += b;
					b -= c;
					c += a;
					a *= 0x9785;//   1001011110000101
					b = rotate16(b, 7);
					c = rotate16(c, 4);
					return a;
				}
				std::string mwlac_varD::get_name() const { return "mwlac_varD"; }
				void mwlac_varD::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varE::raw16() {
					c ^= a;
					a += b;
					b -= c;
					a += c;
					c *= 0x9785;//   1001011110000101
					//shift:	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
					//						37	39-	41	34-	34	38-	39	38-	35	20
					b = rotate16(b, 6);
					return a;
				}
				std::string mwlac_varE::get_name() const { return "mwlac_varE"; }
				void mwlac_varE::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varF::raw16() {
					a += b;
					b ^= c;
					c += a;
					a = rotate16(a, 7);
					b = rotate16(b, 5);
					c *= 0x9785;//   1001011110000101
					Uint16 tmp = a;
					a = b; b = c; c = tmp;
					return a;
				}
				std::string mwlac_varF::get_name() const { return "mwlac_varF"; }
				void mwlac_varF::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}

				Uint32 mwc64x::raw32() {
					Uint32 c = state >> 32;
					Uint32 x = Uint32(state);
					state = x * Uint64(4294883355U) + c;
					return x ^ c;
				}
				std::string mwc64x::get_name() const { return "mwc64x"; }
				void mwc64x::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint64 cxm64_varqual::raw64() {
					const Uint64 K = 0x6595a395a1ec531b;
					Uint64 tmp = high >> 32;
					low += K;
					high += K + ((low < K) ? 1 : 0);
					tmp ^= high ^ 0;//(Uint64)this;
					for (int i = 1; i < num_mult; i++) {
						tmp *= K;
						tmp ^= tmp >> 32;
					}
					tmp *= K;
					return tmp + low;
				}
				std::string cxm64_varqual::get_name() const {
					std::ostringstream str;
					str << "cxm" << num_mult << "n64";
					return str.str();
				}
				void cxm64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
					//walker->handle(num_mult);
				}


				Uint32 mo_Cmfr32::raw32() {
					state = ~(2911329625u * state); state = rotate32(state, 17);
					return state;
				}
				std::string mo_Cmfr32::get_name() const { return "mo_Cmfr32"; }
				void mo_Cmfr32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_Cmr32::raw32() {
					state = 4031235431u * state; state = rotate32(state, 15);
					return state;
				}
				std::string mo_Cmr32::get_name() const { return "mo_Cmr32"; }
				void mo_Cmr32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_Cmr32of64::raw32() {
					state = 38217494031235431ull * state; state = rotate64(state, 37);
					return Uint32(state);
				}
				std::string mo_Cmr32of64::get_name() const { return "mo_Cmr32of64"; }
				void mo_Cmr32of64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}

				Uint32 murmlac32::raw32() {
					Uint32 tmp = state1;
					for (int i = 0; i < rounds; i++) {
						tmp *= 4031235431u;
						tmp ^= tmp >> 16;
					}
					state1 += state2;
					state2 = tmp;
					return state1;
				}
				std::string murmlac32::get_name() const {
					std::ostringstream str;
					str << "murmlac32(" << rounds << ")";
					return str.str();
				}
				void murmlac32::walk_state(StateWalkingObject *walker) {
					walker->handle(state1); walker->handle(state2);
				}

				Uint64 mulcr64::raw64() {
					Uint64 rv = a * count;
					a = rotate64(a, 24) + b;
					count += 2;
					b = rotate64(b, 37) ^ rv;
					return rv;
				}
				std::string mulcr64::get_name() const { return "mulcr64"; }
				void mulcr64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count);
					count |= 1;
				}
				Uint32 mulcr32::raw32() {
					Uint32 rv = a * 2911329625u;
					a = b ^ count++;
					b = rotate32(b, 11) + rv;
					return rv;
				}
				std::string mulcr32::get_name() const { return "mulcr32"; }
				void mulcr32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count);
				}
				Uint16 mulcr16::raw16() {
					Uint16 rv = a * 2911329625u;
					a = b ^ count++;
					b = rotate16(b, 6) + rv;
					return rv;
				}
				std::string mulcr16::get_name() const { return "mulcr16"; }
				void mulcr16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count);
				}

				Uint8 mulcrax8::raw8() {
					Uint8 rv = a * count_low;
					a = b ^ count_high;
					b = rotate8(b, 3) + rv;
					count_low += 2;
					if (count_low == 1) {
						count_high++;
						if (!count_high) count_high = 1;
					}
					return rv;
				}
				std::string mulcrax8::get_name() const { return "mulcrax8"; }
				void mulcrax8::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count_low);
					walker->handle(count_high);
					count_low |= 1;
					if (!count_high) count_high = 1;
				}
				Uint16 mulcrax16::raw16() {
					Uint16 rv = a * count_low;
					a = b ^ count_high;
					b = rotate16(b, 5) + rv;
					count_low += 2;
					if (count_low == 1) {
						count_high++;
						if (!count_high) count_high = 1;
					}
					return rv;
				}
				std::string mulcrax16::get_name() const { return "mulcrax16"; }
				void mulcrax16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count_low);
					walker->handle(count_high);
					count_low |= 1;
					if (!count_high) count_high = 1;
				}
				Uint32 mulcrax32::raw32() {
					Uint32 rv = a * count_low;
					a = b ^ count_high;
					b = rotate32(b, 13) + rv;
					count_low += 2;
					if (count_low == 1) {
						count_high++;
						if (!count_high) count_high = 1;
					}
					return rv;
				}
				std::string mulcrax32::get_name() const { return "mulcrax32"; }
				void mulcrax32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count_low);
					walker->handle(count_high);
					count_low |= 1;
					if (!count_high) count_high = 1;
				}

				Uint16 mulcrx16::raw16() {
					Uint16 rv = a;
					a = rotate16(a * count_low, 6) + count_high;
					count_low += 2;
					if (count_low == 1) {
						count_high++;
						if (!count_high) count_high = 1;
					}
					return rv ^ a;
				}
				std::string mulcrx16::get_name() const { return "mulcrx16"; }
				void mulcrx16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(count_low);
					walker->handle(count_high);
					count_low |= 1;
					if (!count_high) count_high = 1;
				}
				Uint32 mulcrx32::raw32() {
					Uint32 rv = a;
					a = rotate32(a * count_low, 13) + count_high;
					count_low += 2;
					if (count_low == 1) {
						count_high++;
						if (!count_high) count_high = 1;
					}
					return rv ^ a;
				}
				std::string mulcrx32::get_name() const { return "mulcrx32"; }
				void mulcrx32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(count_low);
					walker->handle(count_high);
					count_low |= 1;
					if (!count_high) count_high = 1;
				}
				Uint64 mulcrx64::raw64() {
					Uint64 rv = a;
					a = rotate64(a * count_low, 25) + count_high;
					count_low += 2;
					if (count_low == 1) {
						count_high++;
						if (!count_high) count_high = 1;
					}
					return rv ^ a;
				}
				std::string mulcrx64::get_name() const { return "mulcrx64"; }
				void mulcrx64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(count_low);
					walker->handle(count_high);
					count_low |= 1;
					if (!count_high) count_high = 1;
				}
				Uint16 varrotA::raw16() {
					Uint16 tmp = rotate16(a + b, counter & 15);
					a = b + counter++;
					b ^= c + (c << 3);
					c = tmp;
					b = rotate16(b, 4);
					return tmp;
				}
				std::string varrotA::get_name() const { return "varrotA"; }
				void varrotA::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				Uint16 varrotB::raw16() {
					Uint16 tmp = a + counter++;
					a = b + rotate16(c, 5);
					b = c ^ tmp;
					c = rotate16(tmp, a & 15);
					return b;
				}
				std::string varrotB::get_name() const { return "varrotB"; }
				void varrotB::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				Uint32 varrotC::raw32() {
					a = rotate(a, c & 31) + b;
					b += c;
					c ^= a;
					b = rotate(b, 5);
					return a;
				}
				std::string varrotC::get_name() const { return "varrotC"; }
				void varrotC::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 varrotD::raw16() {
					a += b;
					b += c;
					c ^= a;
					a = rotate(a, b >> 12);
					b = rotate(b, c >> 12);
					return a;
				}
				std::string varrotD::get_name() const { return "varrotD"; }
				void varrotD::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 varrotE::raw16() {
					a += b;
					b += c;
					c = rotate(c, counter & 15) ^ (a + counter);
					b = rotate(b, 5);
					counter++;
					return a;
				}
				std::string varrotE::get_name() const { return "varrotE"; }
				void varrotE::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				Uint16 varrotF::raw16() {
					Uint16 old = a + b;
					a = b + c;
					b = c + counter++;
					c = old;
					a = rotate(a, b >> 12);
					b = rotate(b, c >> 12);
					return a;
				}
				std::string varrotF::get_name() const { return "varrotF"; }
				void varrotF::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				Uint16 lxwm16::raw16() {
					Uint32 tmp = a;
					tmp *= 0xB6C5ul;
					a = (tmp >> 16) ^ b;
					b = c;
					c = Uint16(tmp);
					return b;
				}
				std::string lxwm16::get_name() const { return "lxwm16"; }
				void lxwm16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 lxwm32::raw32() {
					Uint64 tmp = a;
					tmp *= 0x79b9B6C5ull;
					a = (tmp >> 32) ^ b;
					b = Uint32(tmp);
					return a;
				}
				std::string lxwm32::get_name() const { return "lxwm32"; }
				void lxwm32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}
				Uint16 mrsf16::raw16() {
					Uint16 old = a;
					a = b * 0xB635;
					b = rotate16(b, 10) ^ old;
					return old+a;//*/
				}
				std::string mrsf16::get_name() const { return "mrsf16"; }
				void mrsf16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}
				void mrsf16::seed(Uint64 s) {
					if (!s) s = (1ull << 32) - 1;
					a = Uint16(s);
					b = Uint16(s >> 16);
					for (int i = 0; i < 16; i++) raw16();
				}
				Uint32 hierarchyA::raw32() {
					const Uint64 K = 0x9e3779b97f4a7c15ull;
					c = c * b + 137;
					Uint32 tmp = c ^ rotate32(c + a, 21);
					if (!c) {
						set_frame();
					}

					tmp *= K;
					tmp ^= tmp >> 16;
					tmp *= K;
					return tmp;
				}
				void hierarchyA::set_frame() {
					const Uint64 K = 0x9e3779b97f4a7c15ull;
					Uint64 old;
					old = top_a * K;
					top_a = top_b + top_c++;
					top_b = rotate64(top_b, 21) ^ old;
					a = top_a + old;
					old = top_a * K;
					top_a = top_b + top_c++;
					top_b = rotate64(top_b, 21) ^ old;
					b = top_a + old;
					old = top_a * K;
					top_a = top_b + top_c++;
					top_b = rotate64(top_b, 21) ^ old;
					c = top_a + old;

					b |= 1;
					c = Internals::fast_forward_lcg32(Uint32(Sint32(-253)), 0, b, 137);
				}
				std::string hierarchyA::get_name() const { return "hierarchyA"; }
				void hierarchyA::walk_state(StateWalkingObject *walker) {
					walker->handle(top_a); walker->handle(top_b); walker->handle(top_c);
					if (walker->get_properties() & walker->FLAG_SEEDER) {
						set_frame();
						Uint16 how_far = 0;
						walker->handle(how_far);
						how_far %= 253;
						c = Internals::fast_forward_lcg32(how_far, c, b, 137);
					}
					else {
						walker->handle(a); walker->handle(b); walker->handle(c);
					}
					b |= 1;
				}
			}
		}
	}
}