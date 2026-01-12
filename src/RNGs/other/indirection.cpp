#include <string>
#include <sstream>
#include <cstdlib>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include <vector>

#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/arbee.h"
//#include "PractRand/test_helpers.h"

using namespace PractRand::Internals;

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				Uint8 rc4::raw8() {
					b += arr[a];
					Uint8 tmp = arr[b];
					arr[b] = arr[a];
					arr[a] = tmp;
					return arr[Uint8(arr[a++] + arr[b])];
				}
				std::string rc4::get_name() const {return "rc4";}
				void rc4::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					if (walker->is_clumsy() && !walker->is_read_only()) {
						Uint64 seed;
						walker->handle(seed);
						PractRand::RNGs::Raw::arbee seeder(seed);
						for (int i = 0; i < 256; i++) arr[i] = i;
						for (int i = 0; i < 256; i++) {
							Uint8 ai = i, bi = seeder.raw8();
							Uint8 tmp = arr[ai]; arr[ai] = arr[bi]; arr[bi] = tmp;
						}
					}
					else {
						for (int i = 0; i < 256; i++) walker->handle(arr[i]);
					}
				}

				Uint8 rc4_weakenedA::raw8() {
					b += arr[a];
					Uint8 tmp = arr[b];
					arr[b] = arr[a];
					arr[a++] = tmp;
					return tmp;
				}
				std::string rc4_weakenedA::get_name() const { return "rc4_weakenedA"; }
				Uint8 rc4_weakenedB::raw8() {
					b += arr[a];
					Uint8 tmp = arr[b];
					arr[b] = arr[a];
					arr[a++] = tmp;
					return tmp + arr[b];
				}
				std::string rc4_weakenedB::get_name() const { return "rc4_weakenedB"; }
				Uint8 rc4_weakenedC::raw8() {
					b += arr[a];
					Uint8 tmp = arr[b];
					arr[b] = arr[a];
					arr[a++] = tmp;
					return arr[tmp];
				}
				std::string rc4_weakenedC::get_name() const { return "rc4_weakenedC"; }
				Uint8 rc4_weakenedD::raw8() {
					b += arr[a];
					Uint8 tmp = arr[b];
					arr[b] = arr[a];
					arr[a++] = tmp;
					return arr[tmp] + b;
				}
				std::string rc4_weakenedD::get_name() const { return "rc4_weakenedD"; }

				Uint8 ibaa8x::raw8() {
					if (left) {
						return table[--left];
					}
					const int half_size = 1<<(table_size_L2-1);
					const int mask = (1<<table_size_L2)-1;
					Uint8 *base = &table[mask+1];
					for (int i = 0; i <= mask; i++) {
						Uint8 x, y;
						x = base[i];
						a = ((a << 5) | (a >> 3)) + base[(i+half_size) & mask];
						y = base[x & mask] + a + b;
						base[i] = y;
						b = base[(y >> table_size_L2) & mask] + x;
						table[i] = b;
					}
					left = mask;
					return table[mask];
				}
				std::string ibaa8x::get_name() const {
					std::ostringstream tmp;
					int table_size = 1 << table_size_L2;
					tmp << "ibaa8x(" << table_size_L2 << ")";
					return tmp.str();
				}
				void ibaa8x::walk_state(StateWalkingObject *walker) {
					walker->handle(a);walker->handle(b);
					int table_size = 1<<table_size_L2;
					for (int i = 0; i < table_size * 2; i++) walker->handle(table[i]);
					walker->handle(left);
					if (left >= table_size) left = 0;
				}
				ibaa8x::ibaa8x(int table_size_L2_) : table_size_L2(table_size_L2_) {
					table = new Uint8[2 << table_size_L2];
				}
				ibaa8x::~ibaa8x() {delete[] table;}

				Uint16 ibaa16x::raw16() {
					if (left) {
						return table[--left];
					}
					const int half_size = 1<<(table_size_L2-1);
					const int mask = (1<<table_size_L2)-1;
					Uint16 *base = &table[mask+1];
					for (int i = 0; i <= mask; i++) {
						Uint16 x, y;
						x = base[i];
						a = ((a << 11) | (a >> 5)) + base[(i+half_size) & mask];
						y = base[x & mask] + a + b;
						base[i] = y;
						b = base[(y >> table_size_L2) & mask] + x;
						table[i] = b;
					}
					left = mask;
					return table[mask];
				}
				std::string ibaa16x::get_name() const {
					std::ostringstream tmp;
					int table_size = 1<<table_size_L2;
					tmp << "ibaa16x(" << table_size_L2 << ")";
					return tmp.str();
				}
				void ibaa16x::walk_state(StateWalkingObject *walker) {
					walker->handle(a);walker->handle(b);
					int table_size = 1<<table_size_L2;
					for (int i = 0; i < table_size * 2; i++) walker->handle(table[i]);
					walker->handle(left);
					if (left >= table_size) left = 0;
				}
				ibaa16x::ibaa16x(int table_size_L2_) : table_size_L2(table_size_L2_) {
					table = new Uint16[2 << table_size_L2];
				}
				ibaa16x::~ibaa16x() {delete[] table;}

				Uint32 ibaa32x::raw32() {
					if (left) {
						return table[--left];
					}
					const int half_size = 1<<(table_size_L2-1);
					const int mask = (1<<table_size_L2)-1;
					Uint32 *base = &table[mask+1];
					for (int i = 0; i <= mask; i++) {
						Uint32 x, y;
						x = base[i];
						a = ((a << 19) | (a >> 13)) + base[(i+half_size) & mask];
						y = base[x & mask] + a + b;
						base[i] = y;
						b = base[(y >> table_size_L2) & mask] + x;
						table[i] = b;
					}
					left = mask;
					return table[mask];
				}
				std::string ibaa32x::get_name() const {
					std::ostringstream tmp;
					int table_size = 1<<table_size_L2;
					tmp << "ibaa32x(" << table_size_L2 << ")";
					return tmp.str();
				}
				void ibaa32x::walk_state(StateWalkingObject *walker) {
					walker->handle(a);walker->handle(b);
					unsigned long table_size = 1<<table_size_L2;
					for (unsigned long i = 0; i < table_size * 2; i++) walker->handle(table[i]);
					walker->handle(left);
					if (left >= table_size) left = 0;
				}
				ibaa32x::ibaa32x(int table_size_L2_) : table_size_L2(table_size_L2_) {
					table = new Uint32[2 << table_size_L2];
				}
				ibaa32x::~ibaa32x() {delete[] table;}

				
				#define ind32(mm,x)  (*(Uint32 *)(((Uint8 *)(mm)) + ((x) & ((MASK)<<2))))
				#define rngstep32(mix,a,b,mm,m,m2,r,x) \
				{ \
				  x = *m;  \
				  a = (a^(mix)) + *(m2++); \
				  *(m++) = y = ind32(mm,x) + a + b; \
				  *(r++) = b = ind32(mm,y>>table_size_L2) + x; \
				}
				Uint32 isaac32x::raw32() {
					if (left) {
						return table[--left];
					}
					const int HALF_SIZE = 1<<(table_size_L2-1);
					const int MASK = (1<<table_size_L2)-1;
					Uint32 *base = &table[MASK+1];
					Uint32 *m, *m2, *mend, *r;
					Uint32 x, y;
					m = base;
					r = table;
					b += ++c;
					if (table_size_L2 != 2) {
						for (m = base, mend = m2 = m+HALF_SIZE; m<mend; )
						{
							rngstep32( a<<13, a, b, base, m, m2, r, x);
							rngstep32( a>> 6, a, b, base, m, m2, r, x);
							rngstep32( a<< 2, a, b, base, m, m2, r, x);
							rngstep32( a>>16, a, b, base, m, m2, r, x);
						}
						for (m2 = base; m2<mend; )
						{
							rngstep32( a<<13, a, b, base, m, m2, r, x);
							rngstep32( a>> 6, a, b, base, m, m2, r, x);
							rngstep32( a<< 2, a, b, base, m, m2, r, x);
							rngstep32( a>>16, a, b, base, m, m2, r, x);
						}
					}
					else {
						for (m = base, mend = m2 = m+HALF_SIZE; m<mend; )
						{
							rngstep32( a<<13, a, b, base, m, m2, r, x);
							rngstep32( a>> 6, a, b, base, m, m2, r, x);
						}
						for (m2 = base; m2<mend; )
						{
							rngstep32( a<< 2, a, b, base, m, m2, r, x);
							rngstep32( a>>16, a, b, base, m, m2, r, x);
						}
					}
					left = MASK;
					return table[left];
				}
				std::string isaac32x::get_name() const {
					std::ostringstream tmp;
					int table_size = 1<<table_size_L2;
					tmp << "isaac32x(" << table_size_L2 << ")";
					return tmp.str();
				}
				void isaac32x::walk_state(StateWalkingObject *walker) {
					walker->handle(a);walker->handle(b);walker->handle(c);
					unsigned long table_size = 1<<table_size_L2;
					for (unsigned long i = 0; i < table_size*2; i++) walker->handle(table[i]);
					walker->handle(left);
					if (left >= table_size) left = 0;
				}
				isaac32x::isaac32x(int table_size_L2_) : table_size_L2(table_size_L2_) {
					if (table_size_L2 < 2) PractRand::issue_error("invalid table size for isaac32_small");
					table = new Uint32[2 << table_size_L2];
				}
				isaac32x::~isaac32x() {delete[] table;}
				#undef ind32
				#undef rngstep32


				
				#define ind16(mm,x)  (*(Uint16 *)((Uint8 *)(mm) + ((x) & ((MASK)<<1))))
				#define rngstep16(mix,a,b,mm,m,m2,r,x) \
				{ \
				  x = *m;  \
				  a = (a^(mix)) + *(m2++); \
				  *(m++) = y = ind16(mm,x) + a + b; \
				  *(r++) = b = ind16(mm,y>>table_size_L2) + x; \
				}
				Uint16 isaac16x::raw16() {
					if (left) {
						return table[--left];
					}
					const int HALF_SIZE = 1<<(table_size_L2-1);
					const int MASK = (1<<table_size_L2)-1;
					Uint16 *base = &table[MASK+1];
					Uint16 *m, *m2, *mend, *r;
					Uint16 x, y;
					m = base;
					r = table;
					b += ++c;
					if (table_size_L2 != 2) {
						for (m = base, mend = m2 = m+HALF_SIZE; m<mend; )
						{
							//13, 6, 2, 16 -> 7, 3, 2, 5
							rngstep16( a<<7, a, b, base, m, m2, r, x);
							rngstep16( a>>3, a, b, base, m, m2, r, x);
							rngstep16( a<<2, a, b, base, m, m2, r, x);
							rngstep16( a>>5, a, b, base, m, m2, r, x);
						}
						for (m2 = base; m2<mend; )
						{
							rngstep16( a<<7, a, b, base, m, m2, r, x);
							rngstep16( a>>3, a, b, base, m, m2, r, x);
							rngstep16( a<<2, a, b, base, m, m2, r, x);
							rngstep16( a>>5, a, b, base, m, m2, r, x);
						}
					}
					else {
						for (m = base, mend = m2 = m+HALF_SIZE; m<mend; )
						{
							rngstep16( a<<7, a, b, base, m, m2, r, x);
							rngstep16( a>>3 , a, b, base, m, m2, r, x);
						}
						for (m2 = base; m2<mend; )
						{
							rngstep16( a<<2, a, b, base, m, m2, r, x);
							rngstep16( a>>5, a, b, base, m, m2, r, x);
						}
					}
					left = MASK;
					return table[left];
				}
				std::string isaac16x::get_name() const {
					std::ostringstream tmp;
					int table_size = 1<<table_size_L2;
					tmp << "isaac16x(" << table_size_L2 << ")";
					return tmp.str();
				}
				void isaac16x::walk_state(StateWalkingObject *walker) {
					walker->handle(a);walker->handle(b);walker->handle(c);
					unsigned long table_size = 1<<table_size_L2;
					for (unsigned long i = 0; i < table_size*2; i++) walker->handle(table[i]);
					walker->handle(left);
					if (left >= table_size) left = 0;
				}
				isaac16x::isaac16x(int table_size_L2_) : table_size_L2(table_size_L2_) {
					if (table_size_L2 < 2) PractRand::issue_error("invalid table size for isaac16_small");
					table = new Uint16[2 << table_size_L2];
				}
				isaac16x::~isaac16x() {delete[] table;}
				#undef ind16
				#undef rngstep16


				Uint8 efiix8x::raw8() {
					Uint8 iterated = iteration_table  [i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a ^ b;
					a = b + i;
					b = c + indirect;
					c = old + rotate8(c, 3);
					i++;
					return b ^ iterated;//*/

					//1+1: , 1+2: 38-39, 2+2: 38, 1+4: ?, 2+4: 41-42, 
					/*Uint8 iterated = iteration_table  [i & iteration_table_size_m1  ];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated ^ a;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a + i++;
					a = b + iterated;
					b = c ^ indirect;
					c = old + rotate( c, 3 );
					return b;//*/
					
					//"^b" - 1+1: 38, 1+2: >36
					//"^a" - 1+1: 38
					/*Uint8 iterated = iteration_table  [i & iteration_table_size_m1  ] ^ a;
					Uint8 indirect = indirection_table[c & indirection_table_size_m1] + i;
					indirection_table[c & indirection_table_size_m1] = iterated;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a + b;
					a = b + iterated;
					b = c + indirect;
					c = old ^ rotate( c, 3 );
					i++; return old;//*/

					//1+1: 36?, 2+2: 25, 4+4: 33, 8+8: 35, 
					/*Uint8 iterated = iteration_table  [i & iteration_table_size_m1  ] ^ i;
					Uint8 indirect = indirection_table[c & indirection_table_size_m1] + a;
					indirection_table[c & indirection_table_size_m1] = iterated;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a ^ b;
					a = b + indirect;
					b = c + iterated;
					c = old + rotate( c, 3 );
					return b;//*/
				}
				void efiix8x::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(i);
					for (int x = 0; x <= iteration_table_size_m1  ; x++) walker->handle(iteration_table[  x]);
					for (int x = 0; x <= indirection_table_size_m1; x++) walker->handle(indirection_table[x]);
				}
				std::string efiix8x::get_name() const {
					std::ostringstream tmp;
					int itlog2, indlog2;
					for (itlog2 = 0; (1 << itlog2) <= iteration_table_size_m1; itlog2++);
					for (indlog2 = 0; (1 << indlog2) <= indirection_table_size_m1; indlog2++);
					tmp << "efiix8x(" << itlog2 << "," << indlog2 << ")";
					return tmp.str();
				}
				efiix8x::efiix8x(int iteration_table_size_L2, int indirection_table_size_L2) {
					if (iteration_table_size_L2 > 8) PractRand::issue_error("iteration table size log2 too large for efiix8_varqual");
					if (indirection_table_size_L2 > 8) PractRand::issue_error("indirection table size log2 too large for efiix8_varqual");
					int iteration_table_size = 1 << iteration_table_size_L2;
					int indirection_table_size = 1 << indirection_table_size_L2;
					iteration_table_size_m1 = iteration_table_size - 1;
					indirection_table_size_m1 = indirection_table_size - 1;
					iteration_table   = new Uint8[iteration_table_size  ];
					indirection_table = new Uint8[indirection_table_size];
				}
				efiix8x::~efiix8x() {
					delete[] iteration_table;
					delete[] indirection_table;
				}


				Uint8 efiix4x::rotate4(Uint8 value, int bits) {
					value &= 15;
					bits &= 3;
					return (value << bits) | (value >> (4-bits));
				}
				Uint8 efiix4x::raw4() {
					/*
						more recent testing (PR v0.96):

						efiix4, optimal			0		1		2		3		4
							(0,X)				22		26		34		50		82
							(1,X)				26		30		38		54		86
							(2,X)				34		38		46		62		94
							(3,X)				50		54		62

						efiix4, PRstd			0		1		2		3		4
							(0,X)				21		25		25		26		29
							(1,X)				24		28		28		24		29
							(2,X)				31		34		36		31		35
							(3,X)				42		41		47

						efiix4, gjrand			0		1		2		3		4
							(0,X)				CCCC??	9CCC??	--1AC?	--15??	---13?
							(1,X)				19CC??	-3AC??	---19?	1456??	1135??
							(2,X)				---5C?	---15?	----~4	---13?	----1?
							(3,X)				------	-----_	------	-----1	----??

						efiix4, RaBiGeTe		0		1		2		3		4	(only tried 20, 24, and 27 - 8 megabits, 128 megabits, and 1 gigabit)
							(0,X)				20		20		24		24		27
							(1,X)				22		24		27		22		27
							(2,X)				27		pass	pass	27		pass?
							(3,X)				pass	pass?	pass?	pass?	pass?

						efiix4, TestU01			0		1		2		3		4
							(0,X)				14/?/?	6/120/?	3/27/?	~0/14/?	0/2/5
							(1,X)				2/42/?	1/26/?	2/22/?	0/14/?	0/2/1
							(2,X)				0/0/6	0/0/~0	0/1/~0	pass	pass
							(3,X)				pass	pass	pass	pass	pass

					*/
					Uint8 old = a ^ b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//variable rotation variant:
					//	iter. TSL2					indirection TSL2
					//				0			1			2			3			4			SUM 0-4	SUM 0-2
					//	0			23			23			24			26			28
					//	1			24			28			28			30			33
					//	2			31			36			37			39			43
					//	3			>43			
					//	4			
					//	SUM 0-2
					/*Uint8 old = a ^ b;
					Uint8 iterated = rotate4(iteration_table[i & iteration_table_size_m1], b & 3);
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/
					
					//other variable rotation variant:
					//	iter. TSL2					indirection TSL2
					//				0			1			2			3			4			SUM 0-4	SUM 0-2
					//	0			20			24			25			25			29
					//	1			25			26			28			30			32
					//	2			31			35			37			37			>38
					//	3			>38
					//	4			
					//	SUM 0-2
					/*Uint8 old = a ^ b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = rotate4(indirection_table[c & indirection_table_size_m1], b);
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i & iteration_table_size_m1] = indirect;// rotate4(indirect, b & 3);
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//and another:
					//	iter. TSL2					indirection TSL2
					//				0			1			2			3			4			SUM 0-4	SUM 0-2
					//	0			
					//	1			
					//	2			
					//	3						
					//	4			
					//	SUM 0-2
					/*
					Uint8 old = a ^ b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = rotate4(b, iterated) + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//yet another:
					//	iter. TSL2					indirection TSL2
					//				0			1			2			3			4			SUM 0-4	SUM 0-2
					//	0			
					//	1			
					//	2			
					//	3						
					//	4			
					//	SUM 0-2
					/*Uint8 old = rotate4(a+i, b);
					Uint8 iterated = iteration_table[i & iteration_table_size_m1] ^= a;
					Uint8 indirect = indirection_table[c & indirection_table_size_m1] += b;
					a = b + iterated;
					b = c + indirect;
					c = old ^ rotate4(c, 2);
					i++;
					return iterated;//*/


					//this is another I like
					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			24			24			25			28			122			69
					//	2			24			27			28			30			32			141			79
					//	4			31			33			36			38			37			175			100
					//	8			43			38			>42			>42			>41			209+		124+
					//	16			
					//	SUM 1-4		76			84			88			93			97
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i++ & iteration_table_size_m1] = indirect + b;
					Uint8 rv = a + indirect;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					return rv;//*/

					//this version is the best of those that have only two vars in the mixing pool
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			17			22			21			23			24			107			60
					//	2			21			24			24			27			29			125			69
					//	4			27			29			32			36			36			160			88
					//	8			36			39			40			45			37			197			115
					//	16			41			>46			44			>46			41			>219		>131
					//	SUM 1-4		65			75			77			86			89			392			217
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = b + indirect;
					Uint8 tmp = (a + b) ^ indirect;
					a = b + i;
					b = rotate4(b, 1) + tmp;
					return iterated + tmp;//*/

					//it's as if quality runs in a ceiling sometimes, a ceiling not tied to statesize or anything else obvious
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			18			22			21			23			24			108			61
					//	2			19			24			24			27			28			122			67
					//	4			27			28			32			36			35			158			87
					//	8			40			39			41			~46			40			~206		120
					//	16			
					//	SUM 1-4		64			74			77			86			87			388			215
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = b + indirect;
					Uint8 tmp = (a + b) ^ indirect;
					a = b + (b << 2) + i;
					b = rotate4(b, 1) + tmp;
					return iterated + tmp;//*/

					//or... this version does decently for only 2 vars in the mixing, though the structure is different enough I'm not 100% I'd call that a mixing pool
					//				1			2			4			8			16			SUM 1-16	SUM 1-4		SUM 4-16
					//	1			18,18		19,20		21,21		22,24		23,24		103,110		58,62		66,69
					//	2			20,21		23,24		22,25		24,25		28,28		117,125		65,72		74,78
					//	4			27,23		30,28		29,35		35,36		40,39		161,161		86,86		104,110
					//	8			38,25		40,33		42,45		43,40		38,43		201,186		120,103		123,128
					//	16			44,33		44,39		>46,47		>46,45		>46,>46		>228,211+	>134,119	>140,>138
					//	SUM 1-8		103,92		112,105		114,125		124,125		129,134
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = b + i;
					iteration_table[_iter] = a;
					Uint8 tmp = a + indirect;
					a = b + iterated;
					b = rotate4(b, 1) + tmp;
					return iterated ^ tmp;//*/
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = b + i;
					iteration_table[_iter] = a;
					Uint8 tmp = a ^ indirect;
					a = b + iterated;
					b = rotate4(b, 2) + tmp;
					return iterated ^ tmp;//*/

					//similar, but back to 3 vars mixing pool
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			
					//	2			
					//	4			
					//	8			
					//	16			
					//	SUM 1-4		
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = b + i;
					iteration_table[_iter] = a;
					Uint8 tmp = a + b;
					a = b + indirect;
					b = c + iterated;
					c = rotate4(c, 1) + tmp;
					return iterated ^ tmp;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			24			24			25			27			121			69
					//	2			24			27			28			28			29			136			79
					//	4			31			31			38			39			41			180			100
					//	8			35			38			>42			>42			>42			202+		116+
					//	16			>41			>41
					//	SUM 1-4		
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i++ & iteration_table_size_m1] = indirect + b;
					Uint8 rv = c + iterated;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					return rv;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			19			22			24			25			27
					//	2			24			25			26			27			27
					//	4			29			35			38			36			38
					//	8			30			42			>42			>42			>42
					//	16			37			>42
					//	SUM 1-4		
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1] ^ b;
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i++ & iteration_table_size_m1] = indirect;
					Uint8 rv = c + iterated;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					return rv;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			20/20/21	22/24/24	23/25/20	26/24/22	27/26/26	118/119/113	65/69/65
					//	2			23/20/23	21/21/21	22/22/22	26/26/25	31/30/30	123/119/121	66/63/66
					//	4			30/27/30	33/32/34	36/37/34	41/31/38	38/35/37	178/162/173	99/96/98
					//	8			33/32/35	38/34/37	>41/39/40	>41//		38//					113+/105/112
					//	16			>41//		>41//
					//	SUM 1-4		73/67/74	76/77/79	81/84/76	93/81/85	96/91/93
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					iterated ^= a;
					a = b + c;
					b = c + i;
					c = indirect + rotate4(iterated, 1);// 1/2/3
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21/21/21	24/24/24	24/25/26	27/25/26	26/25/26	122/120/123	69/70/71
					//	2			24/24/24	27/28/27	28/29/26	28/22/28	28/22/28	135/125/133	79/81/77
					//	4			27/31/29	33/35/34	35/37/35	37/37/34	36/38/34	168/178/166	95/103/98
					//	8			27/32/29	33/40/35	41/43/39	>42/>42/>40	>42/>42/>40				101/115/103
					//	16			//38		//>40
					//	SUM 1-4		72/76/74	84/87/85	87/91/87	92/84/88	90/85/88
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 old = a + b;//37
					a = b + c;
					b = c + i;
					c = indirect + rotate4(old, 3);
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			20			24			25			27
					//	2			24			27			28			30			31
					//	4			28			34			37			38			37
					//	8			29			35			43			>42			>42
					//	16			39			>42
					//	SUM 1-4		
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 old = a;//38
					a = b + c;
					b = c + i;
					c = indirect + rotate4(old, 1);
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			24			24			26			28
					//	2			23			27			29			31			33
					//	4			24			32			38			42			41
					//	8			25			33			42
					//	16			32			41
					//	SUM 1-4		
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 old = iterated;
					a = b + c;
					b = c + i;
					c = indirect + rotate4(old, 1);
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			
					//	2															31
					//	4									37			31			34
					//	8			32			34
					//	16			
					//	SUM 1-4		
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 old = iterated + a;
					a = b + c;
					b = c + i;
					c = indirect ^ rotate4(old, 2);
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			
					//	2			
					//	4			29			34			34			35			36
					//	8			30			40
					//	16			
					//	SUM 1-4		
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 old = a + b;
					a = rotate4(a, 1) + b;//+1
					b = c + i;
					c = indirect ^ rotate4(old, 2);//^2
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			
					//	2			
					//	4			30/28/30	34/31/34	37/37/35	41/31/39	38/34/38
					//	8			32/32/34	39/34/37
					//	16			
					//	SUM 1-4		
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 old = a + iterated;
					a = b + c;
					b = c + i;
					c = indirect ^ rotate4(old, 3);
					return indirect + a;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			17/			17/			19/			22/			24/
					//	2			21/			24/			23/			26/			26/
					//	4			27/27/27	31/27/31	32/33/23	34/29/30	35/35/35
					//	8			33/			35/			35/			36/			36/
					//	16			39/34/		41/38		41/40		42/35		39/39
					//	SUM 1-4		
					/*int _ind = b & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 tmp = a ^ indirect;
					a = b + i;
					b = rotate4(b, 2) + tmp;
					return iterated + tmp;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			17/			22/			21/			23/			24/
					//	2			21/			24/			23/			26/			26/
					//	4			25/26/		29/27/		33/33/		35/29/		34/35/
					//	8			28/28/		35/33/		38/37/		38/32/		36/36/
					//	16			32/32/		37/36/		40/41/		41/35/		39/39/
					//	SUM 1-4		
					/*int _ind = b & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated + b;
					iteration_table[_iter] = a;
					Uint8 tmp = (a + b) ^ indirect;
					a = b + i;
					b = rotate4(b, 2) + tmp;
					return iterated + tmp;//*/

					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			17/18/18	17/19/16	18/19/18	22/20/22	24/22/24	98/98/98	52/56/52
					//	2			21/21/20	24/22/24	20/24/24	24/24/27	28/27/29	117/118/124	65/67/68
					//	4			27/27/27	29/26/31	28/34/27	34/34/35	38/39/39	156/160/159	84/87/85
					//	8			33/31/33	37/33/35	38/39/31	40/36/37	41/40/41	189/179/177	108/103/99
					//	16			39/34/36	38/40/37	40/>42/38	43/40/40	44/43/44	204/200+/194	117/117+/111
					//	SUM 1-4		65/66/65	70/67/71	66/77/69	80/78/84	90/88/92?	371/376/381	201/210/205
					/*int _ind = a & indirection_table_size_m1;
					int _iter = i++  & iteration_table_size_m1;
					Uint8 iterated = iteration_table[_iter];
					Uint8 indirect = indirection_table[_ind];
					indirection_table[_ind] = iterated;
					iteration_table[_iter] = a;
					Uint8 tmp = a ^ indirect;
					a = b + i;
					b = rotate4(b, 3) + tmp;
					return iterated + tmp;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			24			25			27			30			127			70
					//	2			24			26			28			25			28			131			78
					//	4			31			36			38			30			36			171			105
					//	8			45			41
					//	16			
					//	SUM 1-4		76			86			91			82			94			429
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			22			25			27
					//	2			24			27			28			24
					//	4			28			31
					//	8			
					//	16			
					//	SUM 1-4		
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = b + i;
					b = c + indirect;
					c = old ^ rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			20			23			26			28
					//	2			24			27			29			29			32
					//	4			30			21			22			25			29
					//	8			
					//	16			
					//	SUM 1-4		
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated ^= a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			21			21			23			26			28
					//	2			24			27			29			29			32
					//	4			30			21			22			25			29
					//	8			
					//	16			
					//	SUM 1-4		
					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated ^ a;
					iteration_table[i & iteration_table_size_m1] = indirect;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/


					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1] + i;
					indirection_table[c & indirection_table_size_m1] = iterated + a;
					iteration_table[i++ & iteration_table_size_m1] = indirect;
					a = b + indirect;

					//b = rotate4(b, 1) + c;//4+2: 36
					//c = old + rotate4(c, 3);
					//b = rotate4(b, 2) + c;//4+2: 34
					//c = old + rotate4(c, 1);
					//b = rotate4(b, 2) + c;//4+2: 36
					//c = old + rotate4(c, 3);
					//b = rotate4(b, 3) + c;//4+2: 36
					//c = old + rotate4(c, 3);
					//b = rotate4(b, 3) + c;//4+2: 31
					//c = old + rotate4(c, 2);
					//b = rotate4(b, 3) + c;//4+2: 35
					//c = old + rotate4(c, 1);

					//b = rotate4(b, 1) + c;//4+2: 
					//c = old ^ rotate4(c, 3);
					//b = rotate4(b, 2) + c;//4+2: 
					//c = old ^ rotate4(c, 1);
					//b = rotate4(b, 2) + c;//4+2: 
					//c = old ^ rotate4(c, 3);
					//b = rotate4(b, 3) + c;//4+2: 
					//c = old ^ rotate4(c, 3);
					//b = rotate4(b, 3) + c;//4+2: 
					//c = old ^ rotate4(c, 2);
					//b = rotate4(b, 3) + c;//4+2: 36
					//c = old ^ rotate4(c, 1);

					//b = rotate4(b, 1) ^ c;//4+2: 34
					//c = old + rotate4(c, 3);
					//b = rotate4(b, 2) ^ c;//4+2: 34
					//c = old + rotate4(c, 1);
					//b = rotate4(b, 2) ^ c;//4+2: 34
					//c = old + rotate4(c, 3);
					//b = rotate4(b, 3) ^ c;//4+2: 32
					//c = old + rotate4(c, 3);
					//b = rotate4(b, 3) ^ c;//4+2: 35
					//c = old + rotate4(c, 2);
					//b = rotate4(b, 3) ^ c;//4+2: 35
					//c = old + rotate4(c, 1);

					return b ^ iterated;//*/


					/*Uint8 old = a + b;
					Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated + a;//4+2: 36
					iteration_table[i & iteration_table_size_m1] = indirect;
					//indirection_table[c & indirection_table_size_m1] = iterated + b;//4+2: 35
					//iteration_table[i & iteration_table_size_m1] = indirect;
					//indirection_table[c & indirection_table_size_m1] = iterated + old;//4+2: 32
					//iteration_table[i & iteration_table_size_m1] = indirect;
					//indirection_table[c & indirection_table_size_m1] = iterated ^ old;//4+2: 30
					//iteration_table[i & iteration_table_size_m1] = indirect;

					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 25
					//iteration_table[i & iteration_table_size_m1] = indirect + b;
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 36
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a,1) + rotate4(b, 3) + rotate4(c, 0);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 36
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 1) + rotate4(b, 2) + rotate4(c, 0);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 35
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 1) + rotate4(b, 1) + rotate4(c, 0);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 35
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 3) + rotate4(b, 1) + rotate4(c, 0);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 35
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 3) + rotate4(b, 2) + rotate4(c, 0);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 33
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 3) + rotate4(b, 3) + rotate4(c, 0);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 35
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 0) + rotate4(b, 1) + rotate4(c, 2);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 30
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 0) + rotate4(b, 2) + rotate4(c, 2);
					//indirection_table[c & indirection_table_size_m1] = iterated;//4+2: 34
					//iteration_table[i & iteration_table_size_m1] = indirect + rotate4(a, 0) + rotate4(b, 3) + rotate4(c, 2);

					//indirection_table[c & indirection_table_size_m1] = iterated + a;//4+2: 33
					//iteration_table[i & iteration_table_size_m1] = indirect + b;
					//indirection_table[c & indirection_table_size_m1] = iterated + old;//4+2: 36
					//iteration_table[i & iteration_table_size_m1] = indirect + b;
					//indirection_table[c & indirection_table_size_m1] = iterated ^ old;//4+2: 35
					//iteration_table[i & iteration_table_size_m1] = indirect + b;
					//indirection_table[c & indirection_table_size_m1] = iterated + old;//4+2: 36
					//iteration_table[i & iteration_table_size_m1] = indirect ^ b;
					//indirection_table[c & indirection_table_size_m1] = iterated + old;//4+2: 30
					//iteration_table[i & iteration_table_size_m1] = indirect + b + c;
					a = b + i;
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			18			25			26			26			29			124			69
					//	2			19			25			27			27			31			129			71
					//	4			21			31			38			36			40			166			90
					//	8			25			35			>43			43			>43			192???		>103
					//	16			33			40			47			48?						218???		120
					//	SUM 1-4		58			81			91			89			100			419			230
					//	SUM 1-16	116			156			182+		178+		195+???		826			453+
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + b;
					a = b + i;
					b = c + indirect;
					c = old ^ rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//			002		022		042			102		122		142
					//1+1 --	22		22		21			21		23		21
					//1+2 --	18		18		17			20		19		20
					//1+4 --	21		21		21			22		22		22
					//1+8 --	23		23		23			23		23		23
					//1+16--	26		26		26			26		26		26
					//2+1 --	24		26		24			24		25		24
					//2+2 --	27		29		27			27		29		27
					//2+4 --	31		31		31			~32		33		32
					//2+8 --	29		33		29			30		34		28
					//2+16--	32		36		32			32		35		32
					//4+1 --	24		23		24			31		33		31
					//4+2 --	33		34		33			34		36		34
					//4+4 --	36		36		36			37		39		37
					//4+8 --	38		39		38			38		40		38
					//4+16--	40		42					39		42
					//8+1 --	26		25		26			40		43		??
					//8+2 --	38		38		38			43		??		??
					//16+1--	34		34		34			>42		??		??
					//SUM 1-4	384+98	397+97	382+98		397+133	417+139	395			original:	376	542
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + b;
					a = rotate4(b,1) + (i >> 0);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return indirect;//*/



					//			002		012		022		032		042			102		112		122		132		142
					//1+1 --	21		22
					//1+2 --	24		25
					//1+4 --	24		24
					//1+8 --	25		25
					//1+16--	27		27
					//2+1 --	24		25
					//2+2 --	27		28
					//2+4 --	28		29
					//2+8 --	29		29
					//2+16--	29		30
					//4+1 --	28		29
					//4+2 --	32		32
					//4+4 --	38		39
					//4+8 --	39		40
					//4+16--	41		41
					//8+1 --	29		30
					//8+2 --	38		38
					//16+1--	37		38
					//SUM 1-4				original:	376	542
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					iterated += c;
					Uint8 old = a + b;
					a = rotate4(b, 0) + (i >> 1);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return iterated;//*/
					

					//			002		012		022		032		042			102		112		122		132		142
					//1+1 --	
					//1+2 --	
					//1+4 --	
					//1+8 --	
					//1+16--	
					//2+1 --	
					//2+2 --	
					//2+4 --	
					//2+8 --	
					//2+16--	
					//4+1 --			32
					//4+2 --			35
					//4+4 --			38
					//4+8 --			37
					//4+16--			42
					//8+1 --			>42
					//8+2 --			43
					//16+1--			>42
					//SUM 1-4				original:	376	542
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					iterated += b;
					Uint8 old = a + b;
					a = rotate4(b, 0) + (i >> 1);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return iterated;//*/
					

					//			002		012		022		032		042			102		112		122		132		142
					//1+1 --	
					//1+2 --	
					//1+4 --	
					//1+8 --	
					//1+16--	
					//2+1 --	
					//2+2 --	
					//2+4 --	
					//2+8 --	
					//2+16--	
					//4+1 --			
					//4+2 --			
					//4+4 --	37
					//4+8 --			
					//4+16--			
					//8+1 --			
					//8+2 --			
					//16+1--			
					//SUM 1-4				original:	376	542
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					iterated += a;
					Uint8 old = a + b;
					a = rotate4(b, 0) + (i >> 0);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return iterated;//*/


					//			002		012		022		032		042			102		112		122		132		142
					//1+1 --					. . 22 21
					//1+2 --					. . 22 24
					//1+4 --					. . 22 25
					//1+8 --					. . 24 26
					//1+16--					. . 25 27
					//2+1 --					. . 24 24
					//2+2 --					. . 27 27
					//2+4 --					. . 25 27
					//2+8 --					. . 27 28
					//2+16--	. . .			. . 29 28
					//4+1 --	. 31 31			. . 31 30
					//4+2 --	35 35 35		. . 34 33
					//4+4 --	38 38 38		. . 37 36
					//4+8 --	36 37 38		. . 39 38
					//4+16--	40 >39 40		. . 42 39
					//8+1 --	41				. . . >41
					//8+2 --	>42
					//16+1--	>42
					//SUM 1-4				original:	376	542
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					iterated ^= c;
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + b;
					a = rotate4(b, 1) + (i >> 0);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return iterated + old;//*/


					//			002		022		042			102		122		142
					//1+1 --	21		22		24			23		23		24
					//1+2 --	24		25		26			26		26		28
					//1+4 --	28		29		30			28		28		28
					//1+8 --	25		29		30			30		30		30
					//1+16--	28		32		32			31		31		32
					//2+1 --	24		26		27			25		25		28
					//2+2 --	27		29		31			29		29		29
					//2+4 --	32		33		34			30		30		30
					//2+8 --	29		33		33			32		32		32
					//2+16--	33		36		37			34		34		34
					//4+1 --	24		23		34			33		33		35
					//4+2 --	32		32		37			37		37		37
					//4+4 --	35		38		40			37		37		38
					//4+8 --	35		35		35			38		38		39
					//4+16--	35		41		41			37		39		40
					//8+1 --	27		25		37			46		??		>43
					//8+2 --	34		37		40			>45		??		??
					//16+1--	35		34		>42			>45		??		??
					//SUM 1-4	397	528	422	559	450	612		433	612	433		444			original:	376	542
					//SUM 1-4	384		397		383			397		417		395			original:	376	542
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + b;
					a = rotate4(b,1) + (i >> 4);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return a + indirect;//*/

					
					//			002		022		042			102		122		142
					//1+1 --	21		23
					//1+2 --	24		25
					//1+4 --	24		28
					//1+8 --	25		29
					//1+16--	28		32
					//2+1 --	24		26					24
					//2+2 --	27		28					28
					//2+4 --	28		30					29
					//2+8 --	30		32					29
					//2+16--	32		34					30
					//4+1 --	31		33					31
					//4+2 --	33		33					35
					//4+4 --	36		40					35
					//4+8 --	38		41					36
					//4+16--	37		43					37
					//8+1 --	43		>42
					//8+2 --	
					//16+1--	
					//SUM 1-4	
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					int ind = (a) & indirection_table_size_m1;
					Uint8 indirect = indirection_table[ind];
					indirection_table[ind] = iterated;
					iteration_table[i++ & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + b;
					Uint8 rv = a + indirect;
					a = rotate4(b,1) + (i >> 0);
					b = c + indirect;
					c = old + rotate4(c, 2);
					return rv;//*/


					//			102		122		142
					//1+1 --	23		20
					//1+2 --	22		21
					//1+4 --			22
					//1+8 --			26
					//1+16--			31
					//2+1 --			26		
					//2+2 --			28
					//2+4 --			29
					//2+8 --			31
					//2+16--			34
					//4+1 --	31		32
					//4+2 --	34		36
					//4+4 --	38		39
					//4+8 --	38		41
					//4+16--	39		39
					//8+1 --			>39
					//8+2 --			>39
					//16+1--			>39
					//SUM 1-4			
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					int ind = (a)& indirection_table_size_m1;
					Uint8 indirect = indirection_table[ind];
					indirection_table[ind] = iterated;
					iteration_table[i++ & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + (b << 1);
					Uint8 rv = a + indirect;
					a = b + (i >> 2);
					b = c + indirect;
					c = old + rotate4(c, 2);
					return rv;//*/


					//			002		022		042			102		122		142
					//1+1 --	21		23
					//1+2 --	22		26
					//1+4 --	25		29
					//1+8 --	26		27
					//1+16--	27		28
					//2+1 --	24		25
					//2+2 --	27		28
					//2+4 --	28		29		29
					//2+8 --	29		31		31
					//2+16--	33		34		34
					//4+1 --	31		33		34
					//4+2 --	35		35		36
					//4+4 --	38		39		40
					//4+8 --	35		35		42
					//4+16--	41		41		46
					//8+1 --	>39		>38		46
					//8+2 --			
					//16+1--			
					//SUM 1-4	
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					int ind = (a + i) & indirection_table_size_m1;
					Uint8 indirect = indirection_table[ind];
					indirection_table[ind] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + b;
					a = rotate4(b,0) + (i >> 4);
					b = c + indirect;
					c = old + rotate4(c, 2);
					i++;
					return a + iterated;//*/


					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-4
					//	1			20			25			26			25			27			123			71
					//	2			23			25			27			30			32			137			75
					//	4			24			34			39			39			42			178			97
					//	8			27			39			>45												112+
					//	16			35
					//	SUM 1-4		67			84			92			94			101
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[a & indirection_table_size_m1];
					indirection_table[a & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = indirect + b;
					Uint8 old = a + c;
					a = b + i;
					b = c ^ indirect;
					c = old + rotate4(c, 2);
					i++;
					return b ^ iterated;//*/

					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16	SUM 1-8
					//	1			12			13			16			21			21			83			62
					//	2			13			14			18			24			25			94			69
					//	4			17			19			22			28			34			120			86
					//	8			25			27			32			39			>=45		>=168		123
					//	16			41			43
					//	SUM 1-8		67			73			88			112			>=125		
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated;
					iteration_table[i++ & iteration_table_size_m1] = indirect + b;//a or b makes no real difference here?
					c = rotate4(c, 3);
					a += b + i;
					b ^= c + indirect;
					c += a;
					a = rotate4(a, 2);
					return a ^ iterated;//*/


					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16
					//	1			14			21			26			21			26			108
					//	2			14			23			25			25			28			115
					//	4			14			28			33			32			35			142
					//	8			14			29			36			35			38			152
					//	16			16			39			>44(46?)	40			43			>182
					//	SUM 1-16	72			140			>164		153			170
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[b & indirection_table_size_m1];
					indirection_table[b & indirection_table_size_m1] = iterated;
					iteration_table[i & iteration_table_size_m1] = a;
					Uint8 old = a + i;
					a = b + c;
					b = rotate4(b, 2) ^ (b + c);
					c = indirect + i;
					i++;
					return a ^ iterated;//*/



					//	iter. tsize					indirection table size
					//				1			2			4			8			16			SUM 1-16
					//	1			~16			25			25			27			29			122
					//	2			24			28			30			30			33			145
					//	4			31			34			38			39			41			183
					//	8			31			35			39			39			41			185
					//	16			35			38			39			40			43			195
					//	SUM 1-16	137			160			171			175			187
					/*Uint8 iterated = iteration_table[i & iteration_table_size_m1];
					Uint8 indirect = indirection_table[b & indirection_table_size_m1];
					indirection_table[b & indirection_table_size_m1] = iterated;
					Uint8 old = a + c;
					iteration_table[i & iteration_table_size_m1] = old;
					a += rotate4(a, 2);
					a ^= b + i;
					b ^= rotate4(b, 3);
					b += c + indirect;
					c = rotate4(c, 1);
					c += old;
					i++;
					return indirect ^ iterated;//*/


					//8 - 1+1: , 1+2: 38-39, 2+2: 38, 1+4: ?, 2+4: 41-42, 
					/*Uint8 iterated = iteration_table  [i & iteration_table_size_m1  ];
					Uint8 indirect = indirection_table[c & indirection_table_size_m1];
					indirection_table[c & indirection_table_size_m1] = iterated ^ a;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a + i++;
					a = b + iterated;
					b = c ^ indirect;
					c = old + rotate4( c, 2 );
					return b;//*/
					
					//8 "^b" - 1+1: 38, 1+2: >36
					/*
					"^b" / "^a" on iterated, at 4 bit, shift 3
								1		2		4		8		16
						1		19/21	23/26	24/28	27/29	31/31
						2		26/23	29/29	30/30	32/31	34/35
						4		31/32	32/33	34/36	36/36	
						8		39/42-43
					"^b" / "^a" on iterated, at 4 bit, shift 2
								1		2		4		8		16
						1		19/21	23/18	27/23	29/28	32/31
						2		23/26	24/24	28/28	29/29	32/
						4		32/26	29/30	33/34	34/36
						8		34/34	36/
					"^b" / "^a" on iterated, at 4 bit, shift 1
								1		2		4		8		16
						1		21/21	26/26	26/26	30/29	31/31
						2		25/26	28/28	29/28	32/31	34/35
						4		32/32	31/33	33/36	36/37
						8		37/42?
					*/
					/*Uint8 iterated = iteration_table  [i & iteration_table_size_m1  ] ^ a;
					Uint8 indirect = indirection_table[c & indirection_table_size_m1] + i;
					indirection_table[c & indirection_table_size_m1] = iterated;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a + b;
					a = b + iterated;
					b = c + indirect;
					c = old ^ rotate4( c, 3 );
					i++; return old;//*/

					//8 - 1+1: 36?, 2+2: 25, 4+4: 33, 8+8: 35, 
					/*Uint8 iterated = iteration_table  [i & iteration_table_size_m1  ] ^ i;
					Uint8 indirect = indirection_table[c & indirection_table_size_m1] + a;
					indirection_table[c & indirection_table_size_m1] = iterated;
					iteration_table  [i & iteration_table_size_m1  ] = indirect;
					Uint8 old = a ^ b;
					a = b + indirect;
					b = c + iterated;
					c = old + rotate4( c, 2 );
					return b;//*/
				}
				Uint8 efiix4x::raw8() {
					Uint8 rv = raw4() & 15; return rv | ((raw4() & 15) << 4);
				}
				void efiix4x::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(i);
					for (int x = 0; x <= iteration_table_size_m1  ; x++) walker->handle(iteration_table[  x]);
					for (int x = 0; x <= indirection_table_size_m1; x++) walker->handle(indirection_table[x]);
				}
				int efiix4x::get_native_output_size() const {
					return -1;
				}
				std::string efiix4x::get_name() const {
					std::ostringstream tmp;
					int itlog2, indlog2;
					for (itlog2 = 0; (1 << itlog2) <= iteration_table_size_m1; itlog2++);
					for (indlog2 = 0; (1 << indlog2) <= indirection_table_size_m1; indlog2++);
					tmp << "efiix4x(" << itlog2 << "," << indlog2 << ")";
					return tmp.str();
				}
				efiix4x::efiix4x(int iteration_table_size_L2, int indirection_table_size_L2) {
					if (iteration_table_size_L2 > 8) PractRand::issue_error("iteration table size log2 too large for efiix4_varqual");
					if (indirection_table_size_L2 > 4) PractRand::issue_error("indirection table size log2 too large for efiix4_varqual");
					int iteration_table_size = 1 << iteration_table_size_L2;
					int indirection_table_size = 1 << indirection_table_size_L2;
					iteration_table_size_m1 = iteration_table_size - 1;
					indirection_table_size_m1 = indirection_table_size - 1;
					iteration_table   = new Uint8[iteration_table_size  ];
					indirection_table = new Uint8[indirection_table_size];
				}
				efiix4x::~efiix4x() {
					delete[] iteration_table;
					delete[] indirection_table;
				}

				genindA::genindA(int size_L2) {
					if (size_L2 > 16) issue_error("genindA - size too large");
					if (size_L2 < 0) issue_error("genindA - size too small");
					shift = size_L2;
					table_size_mask = (1 << shift) - 1;
					table = new Uint16[table_size_mask + 1];
				}
				genindA::~genindA() {
					delete[] table;
				}
				Uint16 genindA::raw16() {
					int i1 = a >> (16 - shift);
					int i2 = i & table_size_mask;
					Uint16 o = table[i2];
					table[i2] = a;
					a += table[i1] + o + i++;
					return o;
				}
				void genindA::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(i);
					for (int x = 0; x <= table_size_mask; x++) walker->handle(table[x]);
				}
				std::string genindA::get_name() const {
					std::stringstream buf;
					buf << "genindA(" << shift << ")";
					return buf.str();
				}
				genindB::genindB(int size_L2) {
					if (size_L2 > 16) issue_error("genindB - size too large");
					if (size_L2 < 0) issue_error("genindB - size too small");
					shift = size_L2;
					table_size_mask = (1 << shift) - 1;
					table = new Uint16[table_size_mask + 1];
				}
				genindB::~genindB() {
					delete[] table;
				}
				Uint16 genindB::raw16() {
					int i1 = a >> (16 - shift);
					int i2 = i++ & table_size_mask;
					Uint16 &t1 = table[i1];
					Uint16 &t2 = table[i2];
					Uint16 old = a ^ i;
					a ^= t2 + b;
					b = rotate16(b, 5) + old;
					t1 = t2;
					t2 = old;
					return a;
				}
				void genindB::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(i);
					//for (int x = 0; x < (1 << TABLE_SIZE_L2); x++) walker->handle(table[x]);
					for (int x = 0; x <= table_size_mask; x++) walker->handle(table[x]);
				}
				std::string genindB::get_name() const {
					std::stringstream buf;
					buf << "genindB(" << shift << ")";
					return buf.str();
				}
				genindC::~genindC() {
					delete[] table;
				}
				genindC::genindC(int size_L2) {
					if (size_L2 > 16) issue_error("genindC - size too large");
					if (size_L2 < 1) issue_error("genindC - size too small");
					table_size_L2 = size_L2;
					table = new Uint16[1ull << table_size_L2];
					left = -1;
				}
				Uint16 genindC::refill() {
					int size = 1 << table_size_L2;
					int half_size = 1 << (table_size_L2 - 1);
					int mask = half_size - 1;
					Uint16 *table2 = table + half_size;
					for (int i = 0; i < half_size; i++) {
						Uint16 o = table[i];
						table[i] += a;
						a = table2[o & mask] + rotate16(a, 5);
					}
					for (int i = 0; i < half_size; i++) {
						Uint16 o = table2[i];
						table2[i] += a;
						a = table[o & mask] + rotate16(a, 5);
					}
					left = (1 << table_size_L2) - 1;
					return table[left--];
				}
				void genindC::walk_state(StateWalkingObject *walker) {
					int size = 1 << table_size_L2;
					walker->handle(left);
					if (left >= size) left = -1;
					walker->handle(a);
					for (int x = 0; x < size ; x++) walker->handle(table[x]);
				}
				std::string genindC::get_name() const {
					std::stringstream buf;
					buf << "genindC(" << table_size_L2 << ")";
					return buf.str();
				}
				genindD::~genindD() {
					delete[] table;
				}
				genindD::genindD(int size_L2) {
					if (size_L2 > 16) issue_error("genindD - size too large");
					if (size_L2 < 0) issue_error("genindD - size too small");
					table_size_L2 = size_L2;
					table = new Uint16[1ull << table_size_L2];
					mask = (1ull << table_size_L2) - 1;
				}
				Uint16 genindD::raw16() {
					int i1 = i++ & mask;
					int i2 = a & mask;
					Uint16 tmp = table[i1] ^ a;
					a += tmp;
					table[i1] = table[i2];
					table[i2] = tmp;
					a = rotate(a, 5);
					return tmp;
				}
				void genindD::walk_state(StateWalkingObject *walker) {
					walker->handle(i);
					walker->handle(a);
					for (int i = 0; i <= mask; i++) walker->handle(table[i]);
				}
				std::string genindD::get_name() const {
					std::stringstream buf;
					buf << "genindD(" << table_size_L2 << ")";
					return buf.str();
				}

				genindE::~genindE() {
					delete[] table1;
					delete[] table2;
				}
				genindE::genindE(int size_L2) {
					table_size_L2 = size_L2;
					if (size_L2 < 1) issue_error("genindE - size too small");
					if (size_L2 > 16) issue_error("genindE - size too large");
					table1 = new Uint16[1ull << (table_size_L2)];
					table2 = new Uint16[1ull << (table_size_L2)];
					mask = (1ull << table_size_L2) - 1;
				}
				Uint16 genindE::raw16() {
					int X = i++;
					int Y = a & mask;
					Uint16 A = rotate(table2[X], 3) + rotate(table2[Y], 0);
					Uint16 B = rotate(table1[X], 0) ^ rotate(a, 2);
					table1[X] = a;
					a = B + A;
					if (i > mask) {
						i = 0;
						Uint16 *tmp = table1;
						table1 = table2;
						table2 = tmp;
					}
					return B;
				}
				void genindE::walk_state(StateWalkingObject *walker) {
					walker->handle(i);
					i &= mask;
					walker->handle(a);
					for (int i = 0; i <= mask; i++) walker->handle(table1[i]);
					for (int i = 0; i <= mask; i++) walker->handle(table2[i]);
				}
				std::string genindE::get_name() const {
					std::stringstream buf;
					buf << "genindE(" << table_size_L2 << ")";
					return buf.str();
				}

				genindF::~genindF() {
					delete[] table1;
					delete[] table2;
				}
				genindF::genindF(int size_L2) {
					table_size_L2 = size_L2;
					if (size_L2 < 1) issue_error("genindF - size too small");
					if (size_L2 > 16) issue_error("genindF - size too large");
					table1 = new Uint16[1ull << (table_size_L2)];
					table2 = new Uint16[1ull << (table_size_L2)];
					mask = (1ull << (table_size_L2)) - 1;
				}
				Uint16 genindF::raw16() {
					a ^= table2[i++ & mask];
					int i1 = a & mask;
					a = rotate16(a, 11);
					Uint16 o1 = table1[i1];
					int i2 = o1 & mask;
					Uint16 o2 = table2[i2];
					table1[i1] = o2;
					table2[i2] = o1 + a;
					a += o1 ^ o2;
					return a;
				}
				void genindF::walk_state(StateWalkingObject *walker) {
					walker->handle(i);
					i &= mask;
					walker->handle(a);
					for (int i = 0; i <= mask; i++) walker->handle(table1[i]);
					for (int i = 0; i <= mask; i++) walker->handle(table2[i]);
				}
				std::string genindF::get_name() const {
					std::stringstream buf;
					buf << "genindF(" << table_size_L2 << ")";
					return buf.str();
				}

			}//NotRecommended
		}//Polymorphic
	}//RNGs
}//PractRand
