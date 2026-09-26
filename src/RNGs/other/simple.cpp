#include "PractRand/RNGs/other/simple.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <bit>
#include <cstdint>
#include <sstream>
#include <string>

namespace PractRand::RNGs::Polymorphic::NotRecommended {
	using namespace Internals;
				uint16_t xsalta16x3::raw16() {//slightly more complex output function
					uint16_t tmp = 0, old = 0;
					tmp = a + c + ((b >> 11) | (b << 5));
					old = a;
					a = b ^ (b >> 3);
					b = c ^ (c << 7);
					c = c ^ (c >> 8) ^ old;
					return tmp;
				}
				std::string xsalta16x3::get_name() const { return "xsalta16x3"; }
				void xsalta16x3::walk_state(StateWalkingObject* walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
					if (!a && !b && !c) a = 1;
				}
				uint16_t xsaltb16x3::raw16() {//radically different output function
					uint16_t tmp = 0, old = 0;
					tmp = (a & b) | (b & c) | (c & a);
					old = a;
					a = b ^ (b >> 3);
					b = c ^ (c << 7);
					c = c ^ (c >> 8) ^ old;
					return c + ((tmp << 5) | (tmp >> 11));
				}
				std::string xsaltb16x3::get_name() const { return "xsaltb16x3"; }
				void xsaltb16x3::walk_state(StateWalkingObject* walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
					if (!a && !b && !c) a = 1;
				}
				uint16_t xsaltc16x3::raw16() {//deviating from the standard LFSR state function
					uint16_t old = 0;
					old = a;
					a = b ^ (b >> 5);
					b = c + (c << 3);
					c = c ^ (c >> 7) ^ old;
					return a + b;
				}
				std::string xsaltc16x3::get_name() const { return "xsaltc16x3"; }
				void xsaltc16x3::walk_state(StateWalkingObject* walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
					if (!a && !b && !c) a = 1;
				}

				uint32_t xorshift32::raw32() {
					a ^= a << 13;
					a ^= a >> 17;
					a ^= a << 5;
					return a;
				}
				std::string xorshift32::get_name() const { return "xorshift32"; }
				void xorshift32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					if (!a) a = 1;
				}
				uint64_t xorshift64::raw64() {
					a ^= a << 13;
					a ^= a >> 7;
					a ^= a << 17;
					return a;
				}
				std::string xorshift64::get_name() const { return "xorshift64"; }
				void xorshift64::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					if (!a) a = 1;
				}
				void xorshift64of128::xrs(int bits) {
					if (bits < 64) {
						low ^= low >> bits;
						low ^= high << (64 - bits);
						high ^= high >> bits;
					}
					else if (bits > 64) { low ^= high >> (bits - 64); }
					else { low ^= high; }
				}
				void xorshift64of128::xls(int bits) {
					if (bits < 64) {
						high ^= high << bits;
						high ^= low >> (64 - bits);
						low ^= low << bits;
					}
					else if (bits > 64) { low ^= high >> (bits - 64); }
					else { low ^= high; }
				}
				uint64_t xorshift64of128::raw64() {
					xls(16);
					xrs(53);
					xls(47);
					return low;
				}
				std::string xorshift64of128::get_name() const { return "xorshift64of128"; }
				void xorshift64of128::walk_state(StateWalkingObject* walker) {
					walker->handle(low); walker->handle(high);
					if (!high && !low) low = 1;
				}

				std::string xorshift32of128::get_name() const { return "xorshift32of128"; }
				void xorshift32of128::walk_state(StateWalkingObject* walker) { impl.walk_state(walker); }
				std::string xorshift32of64::get_name() const { return "xorshift32of64"; }
				void xorshift32of64::walk_state(StateWalkingObject* walker) { impl.walk_state(walker); }
				std::string xorshift16of32::get_name() const { return "xorshift16of32"; }
				void xorshift16of32::walk_state(StateWalkingObject* walker) { impl.walk_state(walker); }

				uint32_t xorshift32x4::raw32() {
					uint32_t tmp = x ^ (x << 15);
					x = y;
					y = z;
					z = w;
					tmp ^= tmp >> 4;
					w ^= (w >> 21) ^ tmp;
					return w;
				}
				std::string xorshift32x4::get_name() const { return "xorshift32x4"; }
				void xorshift32x4::walk_state(StateWalkingObject* walker) {
					walker->handle(x);
					walker->handle(y);
					walker->handle(z);
					walker->handle(w);
					if (!(x || y || z || w)) x = 1;
				}

				uint32_t xorwow32plus32::raw32() {
					i += 362437;
					x ^= x << 13;
					x ^= x >> 17;
					x ^= x << 5;
					return x + i;
				}
				std::string xorwow32plus32::get_name() const { return "xorwow32plus32"; }
				void xorwow32plus32::walk_state(StateWalkingObject* walker) {
					walker->handle(x);
					if (!x) x = 1;
					walker->handle(i);
				}
				uint32_t xorwow32of96::raw32() {
					a += 362437;
					return a + impl.raw32();
				}
				std::string xorwow32of96::get_name() const { return "xorwow32of96"; }
				void xorwow32of96::walk_state(StateWalkingObject* walker) {
					impl.walk_state(walker);
					walker->handle(a);
				}
				uint32_t xorwow32x6::raw32() {
					uint32_t tmp = x;
					x = y;
					y = z;
					z = w ^ (w << 1);
					w = v ^ (v >> 7);
					v ^= (v << 4) ^ tmp;
					d += 362437;
					return v + d;
				}
				std::string xorwow32x6::get_name() const { return "xorwow32x6"; }
				void xorwow32x6::walk_state(StateWalkingObject* walker) {
					walker->handle(x);
					walker->handle(y);
					walker->handle(z);
					walker->handle(w);
					walker->handle(v);
					walker->handle(d);
					if (!(x || y || z || w || v)) x = 1;
				}
				uint64_t xorshift128plus::raw64() {
					uint64_t a = state0;
					uint64_t b = state1;

					state0 = b;
					a ^= a << 23;
					a ^= a >> 18;
					a ^= b;
					a ^= b >> 5;
					state1 = a;

					return a + b;
				}
				std::string xorshift128plus::get_name() const { return "xorshift128plus"; }
				void xorshift128plus::walk_state(StateWalkingObject* walker) {
					walker->handle(state0);
					walker->handle(state1);
				}
				uint64_t xoroshiro128plus::raw64() {
					uint64_t result = state0 + state1;
					uint64_t tmp = state0 ^ state1;
					state0 = std::rotl(state0, 55) ^ tmp ^ (tmp << 14);
					state1 = std::rotl(tmp, 36);
					return result;
				}
				std::string xoroshiro128plus::get_name() const { return "xoroshiro128plus"; }
				void xoroshiro128plus::walk_state(StateWalkingObject* walker) {
					walker->handle(state0);
					walker->handle(state1);
				}
				uint64_t xoroshiro128plus_2p64::raw64() {
					uint64_t result = state0 + state1;
					static constexpr uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };

					uint64_t s0 = 0;
					uint64_t s1 = 0;
					for (const auto i : JUMP) {
						for (int b = 0; b < 64; b++) {
							if (i & 1ULL << b) {
								s0 ^= state0;
								s1 ^= state1;
							}
							uint64_t tmp = state0 ^ state1;
							state0 = std::rotl(state0, 55) ^ tmp ^ (tmp << 14);
							state1 = std::rotl(tmp, 36);
						}
					}
					state0 = s0;
					state1 = s1;
					return result;
				}
				std::string xoroshiro128plus_2p64::get_name() const { return "xoroshiro128plus_2p64"; }
				void xoroshiro128plus_2p64::walk_state(StateWalkingObject* walker) {
					walker->handle(state0);
					walker->handle(state1);
				}
				void tinyMT::next_state() {
					uint32_t x = 0, y = 0;
					y = state[3];
					x = (state[0] & 0x7fFFffFF) ^ state[1] ^ state[2];
					x ^= x << 1;
					y ^= (y >> 1) ^ x;
					state[0] = state[1];
					state[1] = state[2];
					state[2] = x ^ (y << 10);
					state[3] = y;
					state[1] ^= -int32_t(y & 1) & state_param1;
					state[2] ^= -int32_t(y & 1) & state_param2;
				}
				uint32_t tinyMT::raw32() {
					next_state();
					uint32_t a = 0, b = 0;
					a = state[3];
					b = state[0] + (state[2] >> 8);
					a ^= b;
					a ^= -int32_t(b & 1) & out_param;
					return a;
				}
				std::string tinyMT::get_name() const { return "tinyMT"; }
				void tinyMT::walk_state(StateWalkingObject* walker) {
					walker->handle(state[0]);
					walker->handle(state[1]);
					walker->handle(state[2]);
					walker->handle(state[3]);
					state_param1 = 0x8f7011ee; // 0x8f7011ee, 0xfc78ff1f, 0x3793fdff is one of many good triples
					state_param2 = 0xfc78ff1f;
					out_param = 0x3793fdff;
				}


				uint32_t sapparot::raw32() {
					uint32_t tmp = 0;
					tmp = a + 0x9e3779b9;
					tmp = (tmp << 7) | (tmp >> 25);
					a = b ^ (~tmp) ^ (tmp << 3);
					a = (a << 7) | (a >> 25);
					b = tmp;
					return a ^ b;
				}
				std::string sapparot::get_name() const { return "sapparot"; }
				void sapparot::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
				}

				uint16_t sap16of48::raw16() {
					uint16_t tmp = 0;
					tmp = a + 0x79b9 + c;
					tmp = (tmp << 5) | (tmp >> 11);
					a = b ^ (~tmp) ^ (tmp << 3);
					a = (a << 5) | (a >> 11);
					b = tmp;
					c = (c + a) ^ b;
					return b;
				}
				std::string sap16of48::get_name() const { return "sap16of48"; }
				void sap16of48::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t sap32of96::raw32() {
					uint32_t tmp = 0;
					tmp = a + 0x9e3779b9 + c;
					tmp = (tmp << 7) | (tmp >> 25);
					a = b ^ (~tmp) ^ (tmp << 3);
					a = (a << 7) | (a >> 25);
					b = tmp;
					c = (c + a) ^ b;
					return b;
				}
				std::string sap32of96::get_name() const { return "sap32of96"; }
				void sap32of96::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				uint32_t flea32x1::raw32() {
					constexpr int SHIFT1 = 15;
					constexpr int SHIFT2 = 27;
					uint32_t e = a[d % SIZE];
					a[d % SIZE] = std::rotl(b, SHIFT1);
					b = c + std::rotl(d, SHIFT2);
					c = d + a[i++ % SIZE];
					d = e + c;
					return b;
				}
				std::string flea32x1::get_name() const { return "flea32x1"; }
				void flea32x1::walk_state(StateWalkingObject* walker) {
					for (auto& z : a) walker->handle(z);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					walker->handle(i);
				}

				/*
				sfc version 1:
				tmp = a ^ counter++;
				a = b + (b << SHIFT1);
				b = ((b << SHIFT2) | (b >> (WORD_BITS - SHIFT2))) + tmp;
				steps:
				load a			load b			load counter
				*	b<<SHIFT1		b<<<SHIFT2		counter + 1		a ^ counter
				+ b				+ a ^ counter	store counter
				store a			store b
				overall:
				quality suffers, but it's fast
				note that (b + (b << SHIFT1)) can be a single LEA on x86 for small shifts
				constants:
				16: 2,5 ?
				32: 5,12 ?
				64: 7,41 ?
				sfc version 2:
				code:
				tmp = a ^ b;
				a = b + (b << SHIFT1);
				b = ((b << SHIFT2) | (b >> (WORD_BITS - SHIFT2))) + tmp + counter++;
				steps:
				load a			load b				load counter
				*	b << SHIFT1		b <<< SHIFT2		a ^ b			counter + 1
				+ b				+((a^b)+counter)***
				store a			store b
				suggested values for SHIFT1,SHIFT2 by output avalanche:
				16 bit: *2,5 - 2,6 - 3,5 - 3,6
				32 bit: 4,9 - 5,9 - *5,12 - 6,10 - 9,13
				64 bit: 6,40 - 7,17 - *7,41 - 10,24
				overall:
				does okay @ 32 & 64 bit, but poorly at 16 bit
				note that (b + (b << SHIFT1)) can be a single LEA on x86 for small shifts
				sfc version 3:
				code
				tmp = a + b + counter++;
				a = b ^ (b >> SHIFT1);
				b = ((b << SHIFT2) | (b >> (WORD_BITS - SHIFT2))) + tmp;
				steps:
				load a			load b			load counter
				b >> SHIFT1		b <<< SHIFT2	a+b+counter***		counter + 1
				*	^ b				+(a+b+counter)	store counter
				store a			store b
				overall:
				good statistical properties
				slower than earlier versions on my CPU, but still fast
				refence values for SHIFT1,SHIFT2:
				16 bit: 2,5
				32 bit: 5,12
				64 bit: 7,41
				sfc version 4:
				code
				Word old = a + b + counter++;
				a = b ^ (b >> SHIFT2);
				b = c + (c << SHIFT3);
				c = old + std::rotl(c,SHIFT1);
				return old;
				steps:
				load a			load b			load c				load counter
				*	b >> SHIFT2		b << SHIFT3		c <<< SHIFT1		counter + 1		a+b+counter***
				^ b				+ b				+(a+b+counter)		store counter
				store a			store b			store c
				overall:
				very good behavior on statistical tests
				uses an extra word / register - not as nice for inlining
				adequate speed
				refence values for SHIFT1,SHIFT2,SHIFT3:
				8 bit:  3,2,1
				16 bit: 7,5,2
				32 bit: 25,8,3
				64 bit: 25,12,3
				*/
				uint16_t sfc_v1_16::raw16() {
					uint16_t tmp = a ^ counter++;
					a = b + (b << 2);
					constexpr int BARREL_SHIFT = 5;
					b = std::rotl(b, BARREL_SHIFT) + tmp;
					return a;
				}
				std::string sfc_v1_16::get_name() const { return "sfc_v1_16"; }
				void sfc_v1_16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				uint32_t sfc_v1_32::raw32() {
					uint32_t tmp = a ^ counter++;
					a = b + (b << 5);
					constexpr int BARREL_SHIFT = 12;
					b = std::rotl(b, BARREL_SHIFT) + tmp;
					return a;
				}
				std::string sfc_v1_32::get_name() const { return "sfc_v1_32"; }
				void sfc_v1_32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				uint16_t sfc_v2_16::raw16() {
					uint16_t tmp = a ^ b;
					a = b + (b << 2);
					constexpr int BARREL_SHIFT = 5;
					b = std::rotl(b, BARREL_SHIFT) + tmp + counter++;
					return tmp;
				}
				std::string sfc_v2_16::get_name() const { return "sfc_v2_16"; }
				void sfc_v2_16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				uint32_t sfc_v2_32::raw32() {
					uint32_t tmp = a ^ b;
					a = b + (b << 5);
					constexpr int BARREL_SHIFT = 12;
					b = std::rotl(b, BARREL_SHIFT) + tmp + counter++;
					return tmp;
				}
				std::string sfc_v2_32::get_name() const { return "sfc_v2_32"; }
				void sfc_v2_32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				uint16_t sfc_v3_16::raw16() {
					uint16_t tmp = a + b + counter++;
					a = b ^ (b >> 2);
					constexpr int BARREL_SHIFT = 5;
					b = std::rotl(b, BARREL_SHIFT) + tmp;
					return tmp;
				}
				std::string sfc_v3_16::get_name() const { return "sfc_v3_16"; }
				void sfc_v3_16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				uint32_t sfc_v3_32::raw32() {
					uint32_t tmp = a + b + counter++;
					a = b ^ (b >> 5);
					constexpr int BARREL_SHIFT = 12;
					b = std::rotl(b, BARREL_SHIFT) + tmp;
					return tmp;
				}
				std::string sfc_v3_32::get_name() const { return "sfc_v3_32"; }
				void sfc_v3_32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				uint16_t jsf16::raw16() {
					uint16_t e = a - ((b << 13) | (b >> 3));
					a = b ^ ((c << 9) | (c >> 7));
					b = c + d;
					c = d + e;
					d = e + a;
					return d;
				}
				std::string jsf16::get_name() const { return "jsf16"; }
				/*seed:
					a = uint16_t(s);
					b = uint16_t(s >> 16);
					c = uint16_t(s >> 32);
					d = uint16_t(s >> 48);
					if (!(a|b) && !(c|d)) d++;
					for (int i = 0; i < 20; i++) raw16();
					*/
				void jsf16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					if (!(a | b) && !(c | d)) d++;
				}

				uint32_t tyche::raw32() {
					a += b;
					d ^= a;
					d = std::rotl(d, 16);
					c += d;
					b ^= c;
					b = std::rotl(b, 12);
					a += b;
					d ^= a;
					d = std::rotl(d, 8);
					c += d;
					b ^= c;
					b = std::rotl(b, 7);
					return b;
				}
				void tyche::seed(uint64_t s) { seed(s, 0); }
				void tyche::seed(uint64_t s, uint32_t idx) {
					a = uint32_t(s >> 32);
					b = uint32_t(s);
					c = 2654435769;
					d = 1367130551 ^ idx;
					for (int i = 0; i < 20; i++) raw32();
				}
				std::string tyche::get_name() const { return "tyche"; }
				void tyche::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					if (!(a | b) && !(c | d)) d++;
				}
				uint16_t tyche16::raw16() {
					a += b;
					d ^= a;
					d = std::rotl(d, 8);
					c += d;
					b ^= c;
					b = std::rotl(b, 6);
					a += b;
					d ^= a;
					d = std::rotl(d, 5);
					c += d;
					b ^= c;
					b = std::rotl(b, 3);
					return b;
				}
				std::string tyche16::get_name() const { return "tyche16"; }
				void tyche16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					if (!(a | b) && !(c | d)) d++;
				}

				uint32_t simpleA::raw32() {
					constexpr int BARREL_SHIFT1 = 19;
					uint32_t tmp = b ^ std::rotl(a, BARREL_SHIFT1);
					a = ~b + c;
					b = c;
					c += tmp;
					return b;
				}
				std::string simpleA::get_name() const { return "simpleA"; }
				void simpleA::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t simpleB::raw16() {
					constexpr int BARREL_SHIFT1 = 3;
					constexpr int BARREL_SHIFT2 = 5;
					uint16_t tmp = std::rotl(a, BARREL_SHIFT1) ^ ~b;
					a = b + c;
					b = c;
					c = tmp + std::rotl(c, BARREL_SHIFT2);
					return tmp;
				}
				std::string simpleB::get_name() const { return "simpleB"; }
				void simpleB::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t simpleC::raw16() {
					constexpr int BARREL_SHIFT1 = 3;
					constexpr int BARREL_SHIFT2 = 5;
					uint16_t tmp = std::rotl(a, BARREL_SHIFT1) ^ ~b;
					a = b + std::rotl(c, BARREL_SHIFT2);
					b = c ^ (c >> 2);
					c += tmp;
					return tmp;
				}
				std::string simpleC::get_name() const { return "simpleC"; }
				void simpleC::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t simpleD::raw32() {
					constexpr int BARREL_SHIFT1 = 19;
					uint32_t old = a;
					uint32_t tmp = b ^ std::rotl(a, BARREL_SHIFT1);
					a = b + c;
					b = c ^ old;
					c = old + tmp;
					return tmp;
				}
				std::string simpleD::get_name() const { return "simpleD"; }
				void simpleD::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t simpleE::raw32() {
					uint32_t old = a + b;
					a = b ^ c;
					b = c + old;
					c = old + ((c << 13) | (c >> 19));
					return old;
				}
				std::string simpleE::get_name() const { return "simpleE"; }
				void simpleE::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t simpleF::raw16() {
					uint16_t old = a ^ b;
					a = b ^ (c & d);
					b = c + old;
					c = ~d;
					d = old + ((d << 5) | (d >> 11));
					return c;
				}
				std::string simpleF::get_name() const { return "simpleF"; }
				void simpleF::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t simpleG::raw32() {
					uint32_t old = a ^ (b >> 7);
					a = b + c + d;
					b = c ^ d;
					c = d + old;
					d = old;
					return a + c;
				}
				std::string simpleG::get_name() const { return "simpleG"; }
				void simpleG::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
				}

				/*
				static uint64_t shift_array64(uint64_t vec[2], unsigned long bits) {
					bits -= 64;
					if (!(bits % 64)) return vec[bits / 64];
					return (vec[bits / 64] << (bits & 63)) | (vec[1 + bits / 64] >> (64 - (bits & 63)));
				}
				*/
				uint32_t trivium_weakenedA::raw32() {
					uint32_t tmp_a = uint32_t(b >> 2) ^ uint32_t(b >> 17);
					uint32_t tmp_b = uint32_t(a >> 2) ^ uint32_t(a >> 23);
					uint32_t new_a = tmp_a ^ uint32_t(a >> 16) ^ (uint32_t(b >> 16) & uint32_t(b >> 18));
					uint32_t new_b = tmp_b ^ uint32_t(b >> 13) ^ (uint32_t(a >> 22) & uint32_t(a >> 24));
					a <<= 32; a |= new_a;
					b <<= 32; b |= new_b;

					return tmp_a ^ tmp_b;
				}
				std::string trivium_weakenedA::get_name() const { return "trivium_weakenedA"; }
				void trivium_weakenedA::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
				}
				uint16_t trivium_weakenedB::raw16() {
					uint16_t tmp_a = uint16_t(c >> 2) ^ uint16_t(c >> 15);
					uint16_t tmp_b = uint16_t(a >> 2) ^ uint16_t(a >> 13);
					uint16_t tmp_c = uint16_t(b >> 5) ^ uint16_t(b >> 10);
					uint16_t new_a = tmp_a ^ uint16_t(a >> 9) ^ (uint16_t(c >> 14) & uint16_t(c >> 16));
					uint16_t new_b = tmp_b ^ uint16_t(b >> 7) ^ (uint16_t(a >> 12) & uint16_t(a >> 14));
					uint16_t new_c = tmp_c ^ uint16_t(c >> 11) ^ (uint16_t(b >> 9) & uint16_t(b >> 11));
					a <<= 16; a |= new_a;
					b <<= 16; b |= new_b;
					c <<= 16; c |= new_c;

					return tmp_a ^ tmp_b ^ tmp_c;
				}
				std::string trivium_weakenedB::get_name() const { return "trivium_weakenedB"; }
				void trivium_weakenedB::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				uint32_t mo_Lesr32::raw32() {
					state = (state << 7) - state; state = std::rotl(state, 23);
					return state;
				}
				std::string mo_Lesr32::get_name() const { return "mo_Lesr32"; }
				void mo_Lesr32::walk_state(StateWalkingObject* walker) {
					walker->handle(state);
				}
				uint32_t mo_ResrRers32::raw32() {
					a = std::rotl(a, 21) - a; a = std::rotl(a, 26);
					b = std::rotl(b, 20) - std::rotl(b, 9);
					return a ^ b;
				}
				std::string mo_ResrRers32::get_name() const { return "mo_ResrRers32"; }
				void mo_ResrRers32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
				}
				uint32_t mo_Rers32of64::raw32() {
					state = std::rotl(state, 8) - std::rotl(state, 29);
					return uint32_t(state);
				}
				std::string mo_Rers32of64::get_name() const { return "mo_Rers32of64"; }
				void mo_Rers32of64::walk_state(StateWalkingObject* walker) {
					walker->handle(state);
				}
				uint32_t mo_Resr32of64::raw32() {
					state = std::rotl(state, 21) - state; state = std::rotl(state, 20);
					return uint32_t(state);
				}
				std::string mo_Resr32of64::get_name() const { return "mo_Resr32of64"; }
				void mo_Resr32of64::walk_state(StateWalkingObject* walker) {
					walker->handle(state);
				}
				uint32_t mo_Resdra32of64::raw32() {
					state = std::rotl(state, 42) - state;
					state += std::rotl(state, 14);
					return uint32_t(state);
				}
				std::string mo_Resdra32of64::get_name() const { return "mo_Resdra32of64"; }
				void mo_Resdra32of64::walk_state(StateWalkingObject* walker) {
					walker->handle(state);
				}
				uint32_t murmlacish::raw32() {
					uint32_t tmp = state1 + (state1 << 3);
					tmp ^= tmp >> 8;
					state1 = std::rotl(state1, 11) + state2;
					state2 += state3 ^ (state3 >> 7) ^ tmp;
					state3 += tmp + (tmp << 3);
					return state1;
				}
				std::string murmlacish::get_name() const {
					std::ostringstream str;
					str << "murmlacish";
					return str.str();
				}
				void murmlacish::walk_state(StateWalkingObject* walker) {
					walker->handle(state1);
					walker->handle(state2);
					walker->handle(state3);
				}

				uint16_t gjishA::raw16() {
					b += a; c = std::rotl(c, 4); a ^= b;
					c += b; a = std::rotl(a, 7); b ^= c;
					a += c; b = std::rotl(b, 11); c ^= a;
					return a;
				}
				std::string gjishA::get_name() const { return "gjishA"; }
				void gjishA::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t gjishB::raw16() {
					b += a; c = std::rotl(c, 4); a ^= b;
					c += b; a = std::rotl(a, 7); b ^= c;
					c += counter++;
					a += c; b = std::rotl(b, 11); c ^= a;
					return a;
				}
				std::string gjishB::get_name() const { return "gjishB"; }
				void gjishB::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				uint32_t gjishC::raw32() {
					b += a; c = std::rotl(c, 21); a ^= b;
					c += b; a = std::rotl(a, 13); b ^= c;
					a += c; b = std::rotl(b, 0); c ^= a;
					return a;
				}
				std::string gjishC::get_name() const { return "gjishC"; }
				void gjishC::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t gjishD::raw32() {
					//			4,10,17		4,10,19		5,7,11		5,7,12		5,7,13		5,7,14		5,7,15		5,8,10		5,8,11		5,8,12		5,8,13		5,8,14		5,8,15
					//--big		3			5			1			1			1			1			4			6			2			4			3			1			4
					//--huge	78						32			10			31			23									67						79			62

					//			6,8,15		6,8,17		6,10,13		6,11,20		6,11,18		6,11,15		6,11,14		6,11,19		6,15,29		6,9,21		6,9,17
					//--big		1			4			7			10			6			2			4			5			7			7			10
					//--huge	41															81

					//			7,9,13		7,9,14	(7,14,9)	7,9,15	(7,15,9)	7,9,16		7,10,15		7,10,14		7,11,17
					//--big		3			1		1			1		3			7			4			4			3
					//--huge	50			43		35			40		52												76

					//			9,13,25		8,13,29
					//--big		7			2
					//--huge				65

					//			21,13,0		21,13,3		21,13,5		19,11,6		18,11,5		23,13,7		23,13,6		21,13,7		12,19,7		12,19,5		14,21,6
					//--big		~1			~9			~9			~5			6			4			7			1?			4			6			11
					//--huge	61																					82
					b += a; c = std::rotl(c, 5); a ^= b;
					c += b; a = std::rotl(a, 8); b ^= c;
					a += c; b = std::rotl(b, 16); c ^= a;
					return a;
				}
				std::string gjishD::get_name() const { return "gjishD"; }
				void gjishD::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				uint16_t ara16::raw16() {
					a += std::rotl(static_cast<uint16_t>(b + c), 3);
					b += std::rotl(static_cast<uint16_t>(c + a), 5);
					c += std::rotl(static_cast<uint16_t>(a + b), 7);
					return a;
				}
				std::string ara16::get_name() const { return "ara16"; }
				void ara16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t ara32::raw32() {
					a += std::rotl(b + c, 7);
					b += std::rotl(c + a, 11);
					c += std::rotl(a + b, 15);
					return a;
				}
				std::string ara32::get_name() const { return "ara32"; }
				void ara32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t arx16::raw16() {
					a ^= std::rotl(static_cast<uint16_t>(b + c), 3);
					b ^= std::rotl(static_cast<uint16_t>(c + a), 5);
					c ^= std::rotl(static_cast<uint16_t>(a + b), 7);
					return a;
				}
				std::string arx16::get_name() const { return "arx16"; }
				void arx16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint32_t arx32::raw32() {
					a ^= std::rotl(b + c, 7);
					b ^= std::rotl(c + a, 11);
					c ^= std::rotl(a + b, 15);
					return a;
				}
				std::string arx32::get_name() const { return "arx32"; }
				void arx32::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t hara16::raw16() {
					a += std::rotl(static_cast<uint16_t>(b + c), 3);
					b ^= std::rotl(static_cast<uint16_t>(c + a), 5);
					c += std::rotl(static_cast<uint16_t>(a + b), 7);
					return a;
				}
				std::string hara16::get_name() const { return "hara16"; }
				void hara16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t harx16::raw16() {
					a ^= std::rotl(static_cast<uint16_t>(b + c), 3);
					b += std::rotl(static_cast<uint16_t>(c + a), 5);
					c ^= std::rotl(static_cast<uint16_t>(a + b), 7);
					return a;
				}
				std::string harx16::get_name() const { return "harx16"; }
				void harx16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				uint16_t learx16::raw16() {
					a ^= std::rotl(static_cast<uint16_t>(b + c), 3);
					b ^= std::rotl(static_cast<uint16_t>(c + (c << 3)), 5);
					c ^= std::rotl(static_cast<uint16_t>(a + (a << 3)), 7);
					return a;
				}
				std::string learx16::get_name() const { return "learx16"; }
				void learx16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				uint16_t hlearx16::raw16() {
					a ^= std::rotl(static_cast<uint16_t>(b + c), 3);
					b ^= std::rotl(static_cast<uint16_t>(c + (c << 3)), 5);
					c += std::rotl(static_cast<uint16_t>(a + (a << 3)), 7);
					return a;
				}
				std::string hlearx16::get_name() const { return "hlearx16"; }
				void hlearx16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				uint16_t alearx16::raw16() {
					a ^= std::rotl(static_cast<uint16_t>(b + c), 3);
					b ^= std::rotl(static_cast<uint16_t>(c + (c << 3)), 5);
					c ^= std::rotl(static_cast<uint16_t>(a + (a << 3)), 7) + b;
					return a;
				}
				std::string alearx16::get_name() const { return "alearx16"; }
				void alearx16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				uint16_t arac16::raw16() {
					a += std::rotl(static_cast<uint16_t>(b + c), 3) + counter++;
					b += std::rotl(static_cast<uint16_t>(c + a), 5);
					c += std::rotl(static_cast<uint16_t>(a + b), 7);
					return a;
				}
				std::string arac16::get_name() const { return "arac16"; }
				void arac16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				uint16_t arxc16::raw16() {
					a ^= std::rotl(static_cast<uint16_t>(b + c), 3) + counter++;
					b ^= std::rotl(static_cast<uint16_t>(c + a), 5);
					c ^= std::rotl(static_cast<uint16_t>(a + b), 7);
					return a;
				}
				std::string arxc16::get_name() const { return "arxc16"; }
				void arxc16::walk_state(StateWalkingObject* walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
}
