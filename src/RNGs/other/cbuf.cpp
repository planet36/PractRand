#include <string>
#include <sstream>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>

#include "PractRand/RNGs/other/mt19937.h"
#include "PractRand/RNGs/other/cbuf.h"

using namespace PractRand;
using namespace PractRand::Internals;

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				Uint8 lfsr_medium::raw8() {
					if (used < SIZE) return cbuf[used++];
					for (int i = 0; i < LAG; i++) {
						cbuf[i] ^= cbuf[i+(SIZE-LAG)] ^ table1[cbuf[i+1]] ^ table2[cbuf[i+2]];
					}
					for (int i = LAG; i < SIZE-2; i++) {
						cbuf[i] ^= cbuf[i-LAG] ^ table1[cbuf[i+1]] ^ table2[cbuf[i+2]];
					}
					cbuf[SIZE-2] ^= cbuf[SIZE-2-LAG] ^ table1[cbuf[SIZE-1]] ^ table2[cbuf[0]];
					cbuf[SIZE-1] ^= cbuf[SIZE-1-LAG] ^ table1[cbuf[0]] ^ table2[1];
					used = 1;
					return cbuf[0];
				}
				std::string lfsr_medium::get_name() const {return "lfsr_medium";}
				void lfsr_medium::walk_state(StateWalkingObject *walker) {
					for (int i = 0; i < SIZE; i++) walker->handle(cbuf[i]);
					walker->handle(used);
					if (used >= SIZE) used = 0;
				}
				lfsr_medium::lfsr_medium() {
					used = 0;
					Uint8 vartaps = 1+2;//255 - 16;
					for (Uint32 i = 0; i < 256; i++) {
						Uint8 low = 0;
						Uint8 high = 0;
						for (int b = 0; b < 8; b++) {
							if ((vartaps >> b) & 1) {
								low ^= i >> b;
								if (b) high ^= i << (8-b);
							}
						}
						table1[i] = low;
						table2[i] = high;
					}
				}

				//Mitchell-Moore: LFib32(Uint32, 55, 24, ADD)
				Uint32 mm32::raw32() {
					Uint32 tmp;
					tmp = cbuf[index1] += cbuf[index2];
					if ( ++index1 == 55 ) index1 = 0;
					if ( ++index2 == 55 ) index2 = 0;
					return tmp;
				}
				std::string mm32::get_name() const {return "mm32";}
				void mm32::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					for (int i = 0; i < 55; i++) walker->handle(cbuf[i]);
					if (index1 >= 55) index1 %= 55;
					index2 = index1 - 24;
					if (index2 >= 55) index2 += 55;//it's an unsigned value
				}
				//Mitchell-Moore modified: LFib16(Uint32, 55, 24, ADD) >> 16
				Uint16 mm16of32::raw16() {
					Uint32 tmp;
					tmp = cbuf[index1] += cbuf[index2];
					if ( ++index1 == 55 ) index1 = 0;
					if ( ++index2 == 55 ) index2 = 0;
					return tmp >> 16;
				}
				std::string mm16of32::get_name() const {return "mm16of32";}
				void mm16of32::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					for (int i = 0; i < 55; i++) walker->handle(cbuf[i]);
					if (index1 >= 55) index1 %= 55;
					index2 = index1 - 24;
					if (index2 >= 55) index2 += 55;//it's an unsigned value
				}
				//Mitchell-Moore modified: LFib32(Uint32, 55, 24, ADD) >> 16
				Uint32 mm32_awc::raw32() {
					Uint32 tmp1, tmp2, tmp3;
					tmp1 = cbuf[index1];
					tmp2 = cbuf[index2];
					tmp3 = tmp1 + tmp2 + carry;
					cbuf[index1] = tmp3;
					carry = (tmp3 - tmp1) >> 31;
					if ((tmp3 == tmp1) && tmp2) carry = 1;
					if (++index1 == 55) index1 = 0;
					if (++index2 == 55) index2 = 0;
					return tmp3;
				}
				std::string mm32_awc::get_name() const { return "mm32_awc"; }
				void mm32_awc::walk_state(StateWalkingObject *walker) {
					walker->handle(carry);
					walker->handle(index1);
					for (int i = 0; i < 55; i++) walker->handle(cbuf[i]);
					if (index1 >= 55) index1 %= 55;
					index2 = index1 - 24;
					if (index2 >= 55) index2 += 55;//it's an unsigned value
					carry &= 1;
				}
				//Mitchell-Moore modified: LFib32(Uint32, 55, 24, ADC) >> 16
				Uint16 mm16of32_awc::raw16() {
					Uint32 tmp1, tmp2, tmp3;
					tmp1 = cbuf[index1];
					tmp2 = cbuf[index2];
					tmp3 = tmp1 + tmp2 + carry;
					cbuf[index1] = tmp3;
					carry = (tmp3 - tmp1) >> 31;
					if ((tmp3 == tmp1) && tmp2) carry = 1;
					if (++index1 == 55) index1 = 0;
					if (++index2 == 55) index2 = 0;
					return Uint16(tmp3);
				}
				std::string mm16of32_awc::get_name() const { return "mm16of32_awc"; }
				void mm16of32_awc::walk_state(StateWalkingObject *walker) {
					walker->handle(carry);
					walker->handle(index1);
					for (int i = 0; i < 55; i++) walker->handle(cbuf[i]);
					if (index1 >= 55) index1 %= 55;
					index2 = index1 - 24;
					if (index2 >= 55) index2 += 55;//it's an unsigned value
					carry &= 1;
				}

				//used by Marsaglia in KISS4691 (2010)
				Uint32 mwc4691::raw32() {
					index = (index < 4691-1) ? index + 1 : 0;
					Uint32 x, t;
					x = cbuf[index];
					t = (x << 13) + carry + x;
					carry = (x>>19) + (t<=x);
					if (!t && !x) carry--;
					cbuf[index] = t;
					return t;
				}
				std::string mwc4691::get_name() const {return "mwc4691";}
				void mwc4691::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(carry);
					for (int i = 0; i < 4691; i++) walker->handle(cbuf[i]);
				}
				
				//
				Uint32 cbuf_accum::raw32() {
					Uint32 tmp = cbuf[--index];
					accum = ((accum << 11) | (accum >> 21)) + ~tmp;
					cbuf[index] = accum;
					if (!index) index = L;
					return accum;
				}
				std::string cbuf_accum::get_name() const { return "cbuf_accum"; }
				void cbuf_accum::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum);
					for (int i = 0; i < L; i++) walker->handle(cbuf[i]);
					if (index >= L) index %= L;
					if (!index) index = L;
				}
				Uint32 cbuf_accum_big::raw32() {
					Uint32 tmp = cbuf[--index];
					accum = ((accum << 11) | (accum >> 21)) + ~tmp;
					cbuf[index] = accum;
					if (!index) index = L;
					return accum;
				}
				std::string cbuf_accum_big::get_name() const { return "cbuf_accum_big"; }
				void cbuf_accum_big::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum);
					for (int i = 0; i < L; i++) walker->handle(cbuf[i]);
					if (index >= L) index %= L;
					if (!index) index = L;
				}
				Uint32 cbuf_2accum_small::raw32() {
					Uint32 tmp = cbuf[--index] + accum2;
					accum2 += accum1;
					enum { SHIFT = 11 };// 3,11 for small, 12,11 for medium
					accum1 = ((accum1 << SHIFT) | (accum1 >> (32-SHIFT))) ^ tmp;
					//		1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
					//	3	22	27	28	28	29	30	31	37	43	42	39	43	37	37	37	16	36				42			37	31	30	29	28	28	26	23
					//	4	23	28	28	29	29	31	32	44?							40	17	40							>40	31	30	29	28	28	26	24
					//	5	24	28	29	29	30	31	33									17
					//	6	25	29	30	30	30	31	33									17
					//	7	27	30	31	31	31	33	34									17
					//	8	34	38	40	41												18
					cbuf[index] = accum2;
					if (!index) index = L;
					return accum2;
				}
				std::string cbuf_2accum_small::get_name() const { return "cbuf_2accum_small"; }
				void cbuf_2accum_small::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum1);
					walker->handle(accum2);
					for (int i = 0; i < L; i++) walker->handle(cbuf[i]);
					if (index >= L) index %= L;
					if (!index) index = L;
				}
				Uint32 cbuf_2accum::raw32() {
					Uint32 tmp = cbuf[--index] + accum2;
					accum2 += accum1;
					accum1 = ((accum1 << 11) | (accum1 >> 21)) ^ tmp;
					cbuf[index] = accum2;
					if (!index) index = L;
					return accum2;
				}
				std::string cbuf_2accum::get_name() const { return "cbuf_2accum"; }
				void cbuf_2accum::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum1);
					walker->handle(accum2);
					for (int i = 0; i < L; i++) walker->handle(cbuf[i]);
					if (index >= L) index %= L;
					if (!index) index = L;
				}
				Uint32 dual_cbuf_small::raw32() {
					Uint32 tmp1, tmp2;
					tmp1 = cbuf1[--index1];
					tmp2 = cbuf2[--index2];
					cbuf1[index1] = tmp1 + tmp2;
					cbuf2[index2] = ((tmp1 << 11) | (tmp1 >> 21)) ^ tmp2;
					if (!index1) index1 = L1;
					if (!index2) index2 = L2;
					return tmp1 + tmp2;
				}
				std::string dual_cbuf_small::get_name() const { return "dual_cbuf_small"; }
				void dual_cbuf_small::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					walker->handle(index2);
					for (int i = 0; i < L1; i++) walker->handle(cbuf1[i]);
					for (int i = 0; i < L2; i++) walker->handle(cbuf2[i]);
					if (index1 > L1) index1 %= L1;
					if (!index1) index1 = L1;
					if (index2 > L2) index2 %= L2;
					if (!index2) index2 = L2;
				}
				Uint32 dual_cbuf::raw32() {
					Uint32 tmp1, tmp2;
					tmp1 = cbuf1[--index1];
					tmp2 = cbuf2[--index2];
					cbuf1[index1] = tmp1 + tmp2;
					cbuf2[index2] = ((tmp1 << 11) | (tmp1 >> 21)) ^ tmp2;
					if (!index1) index1 = L1;
					if (!index2) index2 = L2;
					return tmp1 + tmp2;
				}
				std::string dual_cbuf::get_name() const { return "dual_cbuf"; }
				void dual_cbuf::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					walker->handle(index2);
					for (int i = 0; i < L1; i++) walker->handle(cbuf1[i]);
					for (int i = 0; i < L2; i++) walker->handle(cbuf2[i]);
					if (index1 > L1) index1 %= L1;
					if (!index1) index1 = L1;
					if (index2 > L2) index2 %= L2;
					if (!index2) index2 = L2;
				}
				Uint32 dual_cbuf_big::raw32() {
					Uint32 tmp1, tmp2;
					tmp1 = cbuf1[--index1];
					tmp2 = cbuf2[--index2];
					cbuf1[index1] = tmp1 + tmp2;
					cbuf2[index2] = ((tmp1 << 11) | (tmp1 >> 21)) ^ tmp2;
					if (!index1) index1 = L1;
					if (!index2) index2 = L2;
					return tmp1 + tmp2;
				}
				std::string dual_cbuf_big::get_name() const { return "dual_cbuf_big"; }
				void dual_cbuf_big::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					walker->handle(index2);
					for (int i = 0; i < L1; i++) walker->handle(cbuf1[i]);
					for (int i = 0; i < L2; i++) walker->handle(cbuf2[i]);
					if (index1 > L1) index1 %= L1;
					if (!index1) index1 = L1;
					if (index2 > L2) index2 %= L2;
					if (!index2) index2 = L2;
				}
				Uint32 dual_cbufa_small::raw32() {
					Uint32 tmp1, tmp2;
					tmp1 = cbuf1[--index1];
					tmp2 = cbuf2[--index2];
					accum = ((accum << 11) | (accum >> 21)) + tmp1;
					cbuf1[index1] = tmp1 ^ tmp2;
					cbuf2[index2] = accum;
					if ( !index1 ) index1 = L1;
					if ( !index2 ) index2 = L2;
					return accum;
				}
				std::string dual_cbufa_small::get_name() const { return "dual_cbufa_small"; }
				void dual_cbufa_small::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					walker->handle(index2);
					walker->handle(accum);
					for (int i = 0; i < L1; i++) walker->handle(cbuf1[i]);
					for (int i = 0; i < L2; i++) walker->handle(cbuf2[i]);
					if (index1 > L1) index1 %= L1;
					if (index2 > L2) index2 %= L2;
					if ( !index1 ) index1 = L1;
					if ( !index2 ) index2 = L2;
				}
				Uint32 dual_cbuf_accum::raw32() {
					Uint32 tmp1, tmp2;
					tmp1 = cbuf1[--index1];
					tmp2 = cbuf2[--index2];
					accum = ((accum << 11) | (accum >> 21)) + tmp1;
					cbuf1[index1] = tmp1 ^ tmp2;
					cbuf2[index2] = accum;
					if (!index1) index1 = L1;
					if (!index2) index2 = L2;
					return accum;
				}
				std::string dual_cbuf_accum::get_name() const { return "dual_cbuf_accum"; }
				void dual_cbuf_accum::walk_state(StateWalkingObject *walker) {
					walker->handle(index1);
					walker->handle(index2);
					walker->handle(accum);
					for (int i = 0; i < L1; i++) walker->handle(cbuf1[i]);
					for (int i = 0; i < L2; i++) walker->handle(cbuf2[i]);
					if (index1 > L1) index1 %= L1;
					if (index2 > L2) index2 %= L2;
					if (!index1) index1 = L1;
					if (!index2) index2 = L2;
				}

				Uint16 cbuf3tap_small::raw16() {
					index = (index + 1) & (LENGTH - 1);
					Uint16 tmp1 = cbuf[index];
					Uint16 tmp2 = cbuf[index ^ (LENGTH / 2)];
					Uint16 tmp3 = cbuf[(index + 3) & (LENGTH - 1)];
					cbuf[index] = rotate16(tmp2 + tmp3, 12) ^ rotate16(tmp1 + tmp2, 5);
					return tmp1;
				}
				std::string cbuf3tap_small::get_name() const { return "cbuf3tap_small"; }
				void cbuf3tap_small::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					index &= LENGTH - 1;
					for (int i = 0; i < LENGTH; i++) walker->handle(cbuf[i]);
				}
				Uint16 cbuf3tap::raw16() {
					index = (index + 1) & (LENGTH - 1);
					Uint16 tmp1 = cbuf[index];
					Uint16 tmp2 = cbuf[index ^ (LENGTH / 2)];
					Uint16 tmp3 = cbuf[(index + 3) & (LENGTH - 1)];
					cbuf[index] = rotate16(tmp2 + tmp3, 12) ^ rotate16(tmp1 + tmp2, 5);
					return tmp1;
				}
				std::string cbuf3tap::get_name() const { return "cbuf3tap"; }
				void cbuf3tap::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					index &= LENGTH - 1;
					for (int i = 0; i < LENGTH; i++) walker->handle(cbuf[i]);
				}
				Uint16 cbuf3tap_big::raw16() {
					index = (index + 1) & (LENGTH - 1);
					Uint16 tmp1 = cbuf[index];
					Uint16 tmp2 = cbuf[index ^ (LENGTH / 2)];
					Uint16 tmp3 = cbuf[(index + 3) & (LENGTH - 1)];
					cbuf[index] = rotate16(tmp2 + tmp3, 12) ^ rotate16(tmp1 + tmp2, 5);
					return tmp1;
				}
				std::string cbuf3tap_big::get_name() const { return "cbuf3tap_big"; }
				void cbuf3tap_big::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					index &= LENGTH - 1;
					for (int i = 0; i < LENGTH; i++) walker->handle(cbuf[i]);
				}

				Uint16 cbufa2tap_small::raw16() {
					index &= LENGTH - 1;
					accum += accum << 3;
					Uint16 tmp1 = cbuf[index];
					Uint16 tmp2 = cbuf[index ^ (LENGTH / 2)];
					accum = rotate16(accum, 8) ^ (tmp1 + tmp2);
					cbuf[index] = accum;
					index--;
					return accum;
				}
				std::string cbufa2tap_small::get_name() const { return "cbufa2tap_small"; }
				void cbufa2tap_small::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum);
					for (int i = 0; i < LENGTH; i++) walker->handle(cbuf[i]);
				}
				Uint16 cbufa2tap::raw16() {
					// L	2	4	8	16	32	64
					// Q	22	26	33	35	40	>41
					index &= LENGTH - 1;
					accum += accum << 3;
					Uint16 tmp1 = cbuf[index];
					Uint16 tmp2 = cbuf[index ^ (LENGTH / 2)];
					accum = rotate16(accum, 8) ^ (tmp1 + tmp2);
					cbuf[index] = accum;
					index--;
					return accum;
				}
				std::string cbufa2tap::get_name() const { return "cbufa2tap"; }
				void cbufa2tap::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum);
					for (int i = 0; i < LENGTH; i++) walker->handle(cbuf[i]);
				}
				Uint16 cbufa2tap_big::raw16() {
					index &= LENGTH - 1;
					accum += accum << 3;
					Uint16 tmp1 = cbuf[index];
					Uint16 tmp2 = cbuf[index ^ (LENGTH / 2)];
					accum = rotate16(accum, 8) ^ (tmp1 + tmp2);
					cbuf[index] = accum;
					index--;
					return accum;
				}
				std::string cbufa2tap_big::get_name() const { return "cbufa2tap_big"; }
				void cbufa2tap_big::walk_state(StateWalkingObject *walker) {
					walker->handle(index);
					walker->handle(accum);
					for (int i = 0; i < LENGTH; i++) walker->handle(cbuf[i]);
				}



				Uint32 ranrot32small::raw32() {
					/*if (position < LAG1) return buffer[position++];
					for (unsigned long i = 0; i < LAG2; i++) {
						buffer[i] =
							((buffer[i + LAG1 - LAG1] << ROT1) | (buffer[i + LAG1 - LAG1] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i + LAG1 - LAG2] << ROT2) | (buffer[i + LAG1 - LAG2] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					for (unsigned long i = LAG2; i < LAG1; i++) {
						buffer[i] =
							((buffer[i - 0] << ROT1) | (buffer[i - 0] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i - LAG2] << ROT2) | (buffer[i - LAG2] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					position = 0;
					return buffer[position++];*/
					if (position) return buffer[--position];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] =
							((buffer[i + 0] << ROT1) | (buffer[i + 0] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i - (LAG1 - LAG2)] << ROT2) | (buffer[i - (LAG1 - LAG2)] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] =
							((buffer[i + 0] << ROT1) | (buffer[i + 0] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i + LAG2] << ROT2) | (buffer[i + LAG2] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot32small::get_name() const { return "ranrot32small"; }
				void ranrot32small::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot32::raw32() {
					if (position) return buffer[--position];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] =
							((buffer[i + LAG1 - LAG1] << ROT1) | (buffer[i + LAG1 - LAG1] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i + LAG2 - LAG1] << ROT2) | (buffer[i + LAG2 - LAG1] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] =
							((buffer[i + 0] << ROT1) | (buffer[i + 0] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i + LAG2] << ROT2) | (buffer[i + LAG2] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot32::get_name() const { return "ranrot32"; }
				void ranrot32::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot32big::raw32() {
					if (position) return buffer[--position];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] =
							((buffer[i + LAG1 - LAG1] << ROT1) | (buffer[i + LAG1 - LAG1] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i + LAG2 - LAG1] << ROT2) | (buffer[i + LAG2 - LAG1] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] =
							((buffer[i + 0] << ROT1) | (buffer[i + 0] >> (sizeof(buffer[0]) * 8 - ROT1))) +
							((buffer[i + LAG2] << ROT2) | (buffer[i + LAG2] >> (sizeof(buffer[0]) * 8 - ROT2)));
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot32big::get_name() const { return "ranrot32big"; }
				void ranrot32big::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot3tap32small::func(Uint32 a, Uint32 b, Uint32 c) {
					return rotate(a, ROT1) + rotate(b, ROT2) + rotate(c, ROT3);//30 @ 7/5, 36 @ 11/7, 36 @ 17/9
				}
				Uint32 ranrot3tap32small::raw32() {
					if (position) return buffer[--position];
					/*Uint32 old = buffer[LAG1 - 1];
					for (unsigned long i = 0; i < LAG2; i++) {
						buffer[i] = old = func(buffer[i + LAG1 - LAG1], buffer[i + LAG1 - LAG2], old);
					}
					for (unsigned long i = LAG2; i < LAG1; i++) {
						buffer[i] = old = func(buffer[i], buffer[i - LAG2], old);
					}*/
					Uint32 old = buffer[0];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2 - LAG1], old);
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2], old);
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot3tap32small::get_name() const { return "ranrot3tap32small"; }
				void ranrot3tap32small::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot3tap32::func(Uint32 a, Uint32 b, Uint32 c) {
					return rotate(a, ROT1) + rotate(b, ROT2) + rotate(c, ROT3);//30 @ 7/5, 36 @ 11/7, 36 @ 17/9
				}
				Uint32 ranrot3tap32::raw32() {
					if (position) return buffer[--position];
					Uint32 old = buffer[0];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2 - LAG1], old);
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2], old);
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot3tap32::get_name() const { return "ranrot3tap32"; }
				void ranrot3tap32::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot3tap32big::func(Uint32 a, Uint32 b, Uint32 c) {
					return rotate(a, ROT1) + rotate(b, ROT2) + rotate(c, ROT3);//30 @ 7/5, 36 @ 11/7, 36 @ 17/9
				}
				Uint32 ranrot3tap32big::raw32() {
					if (position) return buffer[--position];
					Uint32 old = buffer[0];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2 - LAG1], old);
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2], old);
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot3tap32big::get_name() const { return "ranrot3tap32big"; }
				void ranrot3tap32big::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot32hetsmall::func(Uint32 a, Uint32 b, Uint32 c) {
					a = rotate(a, ROT1); b = rotate(b, ROT2); c = rotate(c, ROT3);
					return (a + b) ^ c;
				}
				Uint32 ranrot32hetsmall::raw32() {
					if (position) return buffer[--position];
					Uint32 old = buffer[0];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2 - LAG1], old);
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2], old);
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot32hetsmall::get_name() const { return "ranrot32hetsmall"; }
				void ranrot32hetsmall::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot32het::func(Uint32 a, Uint32 b, Uint32 c) {
					a = rotate(a, ROT1); b = rotate(b, ROT2); c = rotate(c, ROT3);
					return (a + b) ^ c;
				}
				Uint32 ranrot32het::raw32() {
					if (position) return buffer[--position];
					Uint32 old = buffer[0];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2 - LAG1], old);
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2], old);
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot32het::get_name() const { return "ranrot32het"; }
				void ranrot32het::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 ranrot32hetbig::func(Uint32 a, Uint32 b, Uint32 c) {
					a = rotate(a, ROT1); b = rotate(b, ROT2); c = rotate(c, ROT3);
					return (a + b) ^ c;
				}
				Uint32 ranrot32hetbig::raw32() {
					if (position) return buffer[--position];
					Uint32 old = buffer[0];
					for (long i = LAG1 - 1; i >= LAG1 - LAG2; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2 - LAG1], old);
					}
					for (long i = LAG1 - LAG2 - 1; i >= 0; i--) {
						buffer[i] = old = func(buffer[i], buffer[i + LAG2], old);
					}
					position = LAG1;
					return buffer[--position];
				}
				std::string ranrot32hetbig::get_name() const { return "ranrot32hetbig"; }
				void ranrot32hetbig::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}

				Uint16 fibmul16of32::raw16() {
					if (position < LAG1) return Uint16(buffer[position++] >> 16);
					for (unsigned long i = 0; i < LAG2; i++) {
						buffer[i] = buffer[i] * buffer[i+LAG1-LAG2];
					}
					for (unsigned long i = LAG2; i < LAG1; i++) {
						buffer[i] = buffer[i] * buffer[i-LAG2];
					}
					position = 0;
					return Uint16(buffer[position++] >> 16);
				}
				std::string fibmul16of32::get_name() const {return "fibmul16of32";}
				void fibmul16of32::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
					for (int i = 0; i < LAG1; i++) buffer[i] |= 1;
				}
				Uint32 fibmul32of64::raw32() {
					if (position < LAG1) return Uint32(buffer[position++] >> 32);
					for (unsigned long i = 0; i < LAG2; i++) {
						buffer[i] = buffer[i+LAG1-LAG1] * buffer[i+LAG1-LAG2];
					}
					for (unsigned long i = LAG2; i < LAG1; i++) {
						buffer[i] = buffer[i] * buffer[i-LAG2];
					}
					position = 0;
					return Uint32(buffer[position++] >> 32);
				}
				std::string fibmul32of64::get_name() const {return "fibmul32of64";}
				void fibmul32of64::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
					for (int i = 0; i < LAG1; i++) buffer[i] |= 1;
				}
				Uint16 fibmulmix16::raw16() {
					if (position < LAG1) return buffer[position++];
					Uint16 prev = buffer[LAG1 - 1];
					enum {
						SH1  = 0
						,SH2 = 0
						,SH3 = 5
					};
					for (unsigned long i = 0; i < LAG2; i++) {
						Uint16 a = buffer[i + LAG1 - LAG1];
						Uint16 b = buffer[i + LAG1 - LAG2];
						a = rotate16(a, SH1); b = rotate16(b, SH2); prev += rotate16(prev, SH3);
						Uint16 product = a * (b | 1);
						//product ^= product >> 8; // fixes BRank failure
						prev ^= product;
						//prev ^= prev >> 8; // does NOT fix BRank failure
						buffer[i] = prev;
					}
					for (unsigned long i = LAG2; i < LAG1; i++) {
						Uint16 a = buffer[i];
						Uint16 b = buffer[i - LAG2];
						a = rotate16(a, SH1); b = rotate16(b, SH2); prev += rotate16(prev, SH3);
						Uint16 product = a * (b | 1);
						//product ^= product >> 8; // fixes BRank failure
						prev ^= product;
						//prev ^= prev >> 8; // does NOT fix BRank failure
						buffer[i] = prev;
					}
					position = 0;
					return buffer[position++];
				}
				std::string fibmulmix16::get_name() const { return "fibmulmix16"; }
				void fibmulmix16::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}
				Uint32 fibmulmix32::raw32() {
					if (position < LAG1) return buffer[position++];
					Uint32 prev = buffer[LAG1 - 1];
					enum {
						SH1 = 3
						, SH2 = 5
						, SH3 = 19
					};
					for (unsigned long i = 0; i < LAG2; i++) {
						Uint32 a = buffer[i + LAG1 - LAG1];
						Uint32 b = buffer[i + LAG1 - LAG2];
						a = rotate32(a, SH1); b = rotate32(b, SH2); prev += rotate32(prev, SH3);
						Uint32 product = a * (b | 1);
						//product ^= product >> 16; // fixes BRank failure
						prev ^= product;
						//prev ^= prev >> 16; // does NOT fix BRank failure
						buffer[i] = prev;
					}
					for (unsigned long i = LAG2; i < LAG1; i++) {
						Uint32 a = buffer[i];
						Uint32 b = buffer[i - LAG2];
						a = rotate32(a, SH1); b = rotate32(b, SH2); prev += rotate32(prev, SH3);
						Uint32 product = a * (b | 1);
						//product ^= product >> 16; // fixes BRank failure
						prev ^= product;
						//prev ^= prev >> 16; // does NOT fix BRank failure
						buffer[i] = prev;
					}
					position = 0;
					return buffer[position++];
				}
				std::string fibmulmix32::get_name() const { return "fibmulmix32"; }
				void fibmulmix32::walk_state(StateWalkingObject *walker) {
					walker->handle(position);
					for (int i = 0; i < LAG1; i++) walker->handle(buffer[i]);
					if (position >= LAG1) position %= LAG1;
				}


				Uint32 mt19937_unhashed::raw32() {
					return implementation.untempered_raw32();
				}
				std::string mt19937_unhashed::get_name() const {return "mt19937_unhashed";}
				void mt19937_unhashed::walk_state(StateWalkingObject *walker) {
					implementation.walk_state(walker);
				}

				void chacha_weakenedA::refill() {
					for (int i = 0; i < 9; i++) buffer[i] = state[i];
					for (int i = 0; i < quality; i++) {
						// chacha/salsa use 4x4 matrix of 32 bit values, here I'll use a 3x3 matrix of 8 bit values to reduce quality
						Uint8 tmp, offset;
						offset = 0;
						tmp = buffer[0 + offset];
						buffer[0 + offset] = buffer[1 + offset] + buffer[2 + offset];
						buffer[1 + offset] = buffer[2 + offset] ^ tmp;
						buffer[2 + offset] = rotate8(tmp, 3);
						offset = 3;
						tmp = buffer[0 + offset];
						buffer[0 + offset] = buffer[1 + offset] + buffer[2 + offset];
						buffer[1 + offset] = buffer[2 + offset] ^ tmp;
						buffer[2 + offset] = rotate8(tmp, 4);
						offset = 6;
						tmp = buffer[0 + offset];
						buffer[0 + offset] = buffer[1 + offset] + buffer[2 + offset];
						buffer[1 + offset] = buffer[2 + offset] ^ tmp;
						buffer[2 + offset] = rotate8(tmp, 5);
						//buffer[0] += buffer[3]; buffer[3] += buffer[6]; buffer[7] += buffer[4]; buffer[4] += buffer[1];
						i++;
						if (i >= quality) break;
						offset = 0;
						tmp = buffer[3 + offset];
						buffer[3 + offset] = buffer[6 + offset] + buffer[2 + offset];
						buffer[6 + offset] = buffer[0 + offset] ^ tmp;
						buffer[0 + offset] = rotate8(tmp, 3);
						offset = 1;
						tmp = buffer[3 + offset];
						buffer[3 + offset] = buffer[6 + offset] + buffer[2 + offset];
						buffer[6 + offset] = buffer[0 + offset] ^ tmp;
						buffer[0 + offset] = rotate8(tmp, 4);
						offset = 2;
						tmp = buffer[3 + offset];
						buffer[3 + offset] = buffer[6 + offset] + buffer[2 + offset];
						buffer[6 + offset] = buffer[0 + offset] ^ tmp;
						buffer[0 + offset] = rotate8(tmp, 5);
					}
					/*
						9 outputs per buffer, revised version
							Quality			3	4	5	6	7	8	9	10	11	12	13	14	15	16
							PRstd			10	11	12	12	15	18	18	22	25	28	30	36	39	43+
							gj tiny/sm/st	DDD	CCD	CCC	BBC	BBC	AAB	AAB	47B	46A	135	-27	---	---	---
							gj big/huge												8	2	-1	--
							gj tera/TT															-?
							diehard 0/1/2	2	2	2	2	2	2	2	2	22	-~2	112	---	---	---
							diehard a										3	-	-	-	-	-
							TestU01 SC/C/BC	F	F	F	D	A	5	6	3!	~!	-6	-9!	---	--	

						PRstd			10	10	15	16	25	26	38						9 outputs per buffer (canonical for this PRNG)
						PRstd			10	13	15	16	18	27	35						8 outputs per buffer
					*/
				}
				Uint8 chacha_weakenedA::advance_and_refill_helper() {
					if (!++state[0]) if (!++state[4]) if (!++state[8]) if (!++state[3]) if (!++state[7]) if (!++state[2]) ++state[6];//cycle length is 2**56 times however many outputs are read from each buffer
					refill();
					bufpos = 1;
					return buffer[0];
				}
				std::string chacha_weakenedA::get_name() const { return std::string("chacha_weakenedA(") + std::to_string(quality) + ")"; }
				void chacha_weakenedA::walk_state(StateWalkingObject *walker) {
					for (int i = 0; i < 9; i++) walker->handle(state[i]);
					for (int i = 0; i < 9; i++) walker->handle(buffer[i]);
					walker->handle(bufpos);
					if (bufpos > 8) bufpos = 8;
					//walker->handle(quality);
				}

				void chacha_weakenedB::refill() {
					for (int i = 0; i < 8; i++) buffer[i] = state[i];
					for (int i = 0; i < quality; i++) {
						// matrix is 2x2x2 this time ; sure that's an extra dimension, but it's still smaller
						buffer[0] += buffer[1]; buffer[1] += buffer[2]; buffer[2] += buffer[3]; buffer[3] ^= buffer[0]; buffer[0] = rotate8(buffer[0], 3);
						buffer[4] += buffer[5]; buffer[5] += buffer[6]; buffer[6] += buffer[7]; buffer[7] ^= buffer[4]; buffer[5] = rotate8(buffer[5], 3);
						if (++i >= quality) break;
						buffer[0] += buffer[2]; buffer[2] += buffer[4]; buffer[4] ^= buffer[6]; buffer[6] += buffer[0]; buffer[0] = rotate8(buffer[0], 5);
						buffer[1] += buffer[3]; buffer[3] += buffer[5]; buffer[5] ^= buffer[7]; buffer[7] += buffer[1]; buffer[5] = rotate8(buffer[5], 4);
						if (++i >= quality) break;
						buffer[0] += buffer[4]; buffer[4] ^= buffer[1]; buffer[1] += buffer[5]; buffer[5] += buffer[0]; buffer[1] = rotate8(buffer[1], 4);
						buffer[2] += buffer[6]; buffer[6] ^= buffer[3]; buffer[3] += buffer[7]; buffer[7] += buffer[2]; buffer[3] = rotate8(buffer[3], 3);
					}
					for (int i = 0; i < 8; i++) buffer[i] += state[i];
					/*
						Quality		2	3	4	5	6	7	8	9	10	11	12
						PRstd		10	12	11	12	15	18	22	27	35	42	
					*/
				}
				Uint8 chacha_weakenedB::advance_and_refill_helper() {
					if (!++state[0]) if (!++state[1]) if (!++state[2]) if (!++state[3]) if (!++state[4]) if (!++state[5]) ++state[6];//cycle length is 2**56 times however many outputs are read from each buffer
					refill();
					bufpos = 1;
					return buffer[0];
				}
				std::string chacha_weakenedB::get_name() const { return std::string("chacha_weakenedB(") + std::to_string(quality) + ")"; }
				void chacha_weakenedB::walk_state(StateWalkingObject *walker) {
					for (int i = 0; i < 8; i++) walker->handle(state[i]);
					for (int i = 0; i < 8; i++) walker->handle(buffer[i]);
					walker->handle(bufpos);
					if (bufpos > 8) bufpos = 8;
					//walker->handle(quality);
				}
			}
		}
	}
}
