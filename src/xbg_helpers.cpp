#include <string>
//#include <ostream>
//#include <sstream>
#include <vector>
//#include <list>
#include <set>
#include <map>
#include <cmath>
#include <memory>
#include <algorithm>
#include <cstdlib>

#include "PractRand.h"
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/test_helpers.h"
#include "PractRand/xbg_helpers.h"
//#include "PractRand/RNGs/efiix64x48.h"

#if defined _MSC_VER && _MSC_VER >= 1400
#include <intrin.h>
#endif


namespace PractRand {
	namespace Internals {
		namespace XBG_helpers {
			namespace BinarySignalOperations {
				void xor64x(long length, const Uint64 *k, Uint64 *value) {
					for (int i = 0; i < length; i++) value[i] ^= k[i];
				}
				void xor32x(long length, const Uint32 *k, Uint32 *value) {
					for (int i = 0; i < length; i++) value[i] ^= k[i];
				}
				void xor16x(long length, const Uint16 *k, Uint16 *value) {
					for (int i = 0; i < length; i++) value[i] ^= k[i];
				}

				void x2mod64x(long length, const Uint64 *mod, Uint64 *value) {
					// result is written to the same memory that the value was passed in
					// this basically leftshifts value by 1, and if the overflow bit was 1 it xors it by mod
					enum { WORDSIZE = 8 * sizeof(*mod) };
					int i;
					length -= 1;
					if (value[length] >> (WORDSIZE - 1)) {
						for (i = length; i; i--) value[i] = ((value[i - 1] >> (WORDSIZE - 1)) | (value[i] << 1)) ^ mod[i];
						value[0] = (value[0] << 1) ^ mod[i];
					}
					else {
						for (i = length; i; i--) value[i] = (value[i - 1] >> (WORDSIZE - 1)) | (value[i] << 1);
						value[0] <<= 1;
					}
				}
				void x2mod32x(long length, const Uint32 *mod, Uint32 *value) {
					enum { WORDSIZE = 8 * sizeof(*mod) };
					int i;
					length -= 1;
					if (value[length] >> (WORDSIZE - 1)) {
						for (i = length; i; i--) value[i] = ((value[i - 1] >> (WORDSIZE - 1)) | (value[i] << 1)) ^ mod[i];
						value[0] = (value[0] << 1) ^ mod[i];
					}
					else {
						for (i = length; i; i--) value[i] = (value[i - 1] >> (WORDSIZE - 1)) | (value[i] << 1);
						value[0] <<= 1;
					}
				}
				void x2mod16x(long length, const Uint16 *mod, Uint16 *value) {
					enum { WORDSIZE = 8 * sizeof(*mod) };
					int i;
					length -= 1;
					if (value[length] >> (WORDSIZE - 1)) {
						for (i = length; i; i--) value[i] = ((value[i - 1] >> (WORDSIZE - 1)) | (value[i] << 1)) ^ mod[i];
						value[0] = (value[0] << 1) ^ mod[i];
					}
					else {
						for (i = length; i; i--) value[i] = (value[i - 1] >> (WORDSIZE - 1)) | (value[i] << 1);
						value[0] <<= 1;
					}
				}

				void mulmod64x(long length, const Uint64* mod, const Uint64* a, const Uint64* b, Uint64* result) {
					typedef Uint64 Word;
					enum { WORDSIZE = 8 * sizeof(Word) };
					if (length <= 0) return;
					int highest_word = length - 1;
					while (!a[highest_word]) {
						highest_word--;
						if (highest_word < 0) return;
					}
					if (false) {
						std::vector<Word> _tmpvec1; _tmpvec1.resize(length * 2);
						Uint64* tmp1 = &_tmpvec1[0];
						Uint64* tmpresult = &_tmpvec1[length];
						for (int i = 0; i < length; i++) tmp1[i] = b[i];
						for (int i = 0; i < length; i++) tmpresult[i] = 0;
						for (int k = 0; k <= highest_word; k++) {
							Word cur = a[k];
							for (int j = WORDSIZE; j; j--) {
								if (cur & 1) xor64x(length, tmp1, tmpresult);
								cur >>= 1;
								x2mod64x(length, mod, tmp1);
							}
						}
						for (int i = 0; i < length; i++) result[i] = tmpresult[i];
					}
					else {
						std::vector<Word> _tmpvec; _tmpvec.resize(length);
						Uint64* tmp = &_tmpvec[0];
						for (int i = 0; i < length; i++) tmp[i] = b[i];
						int pos = PractRand::Internals::ilog2_64(a[highest_word]);// a[highest_word] is guaranteed to be non-zero
						Word testbit = 1ull << pos;
						testbit >>= 1;
						while (highest_word >= 0) {
							while (testbit) {
								x2mod64x(length, mod, tmp);
								if (a[highest_word] & testbit) xor64x(length, b, tmp);
								testbit >>= 1;
							}
							testbit = 1ull << (WORDSIZE - 1);
							highest_word -= 1;
						}
						for (int i = 0; i < length; i++) result[i] = tmp[i];
					}
				}
				void mulmod32x(long length, const Uint32* mod, const Uint32* a, const Uint32* b, Uint32* result) {
					typedef Uint32 Word;
					enum { WORDSIZE = 8 * sizeof(Word) };
					if (length <= 0) return;
					int highest_word = length - 1;
					while (!a[highest_word]) {
						highest_word--;
						if (highest_word < 0) return;
					}
					if (false) {
						std::vector<Word> _tmpvec1; _tmpvec1.resize(length * 2);
						Word* tmp1 = &_tmpvec1[0];
						Word* tmpresult = &_tmpvec1[length];
						for (int i = 0; i < length; i++) tmp1[i] = b[i];
						for (int i = 0; i < length; i++) tmpresult[i] = 0;
						for (int k = 0; k <= highest_word; k++) {
							Word cur = a[k];
							for (int j = WORDSIZE; j; j--) {
								if (cur & 1) xor32x(length, tmp1, tmpresult);
								cur >>= 1;
								x2mod32x(length, mod, tmp1);
							}
						}
						for (int i = 0; i < length; i++) result[i] = tmpresult[i];
					}
					else {
						std::vector<Word> _tmpvec; _tmpvec.resize(length);
						Word* tmp = &_tmpvec[0];
						for (int i = 0; i < length; i++) tmp[i] = b[i];
						int pos = PractRand::Internals::ilog2_32(a[highest_word]);// a[highest_word] is guaranteed to be non-zero
						Word testbit = 1ull << pos;
						testbit >>= 1;
						while (highest_word >= 0) {
							while (testbit) {
								x2mod32x(length, mod, tmp);
								if (a[highest_word] & testbit) xor32x(length, b, tmp);
								testbit >>= 1;
							}
							testbit = 1ull << (WORDSIZE - 1);
							highest_word -= 1;
						}
						for (int i = 0; i < length; i++) result[i] = tmp[i];
					}
				}
				void mulmod16x(long length, const Uint16* mod, const Uint16* a, const Uint16* b, Uint16* result) {
					typedef Uint16 Word;
					enum { WORDSIZE = 8 * sizeof(Word) };
					if (length <= 0) return;
					int highest_word = length - 1;
					while (!a[highest_word]) {
						highest_word--;
						if (highest_word < 0) return;
					}
					if (false) {
						std::vector<Word> _tmpvec1; _tmpvec1.resize(length * 2);
						Word* tmp1 = &_tmpvec1[0];
						Word* tmpresult = &_tmpvec1[length];
						for (int i = 0; i < length; i++) tmp1[i] = b[i];
						for (int i = 0; i < length; i++) tmpresult[i] = 0;
						for (int k = 0; k <= highest_word; k++) {
							Word cur = a[k];
							for (int j = WORDSIZE; j; j--) {
								if (cur & 1) xor16x(length, tmp1, tmpresult);
								cur >>= 1;
								x2mod16x(length, mod, tmp1);
							}
						}
						for (int i = 0; i < length; i++) result[i] = tmpresult[i];
					}
					else {
						std::vector<Word> _tmpvec; _tmpvec.resize(length);
						Word* tmp = &_tmpvec[0];
						for (int i = 0; i < length; i++) tmp[i] = b[i];
						int pos = PractRand::Internals::ilog2_32(a[highest_word]);// a[highest_word] is guaranteed to be non-zero
						Word testbit = 1ull << pos;
						testbit >>= 1;
						while (highest_word >= 0) {
							while (testbit) {
								x2mod16x(length, mod, tmp);
								if (a[highest_word] & testbit) xor16x(length, b, tmp);
								testbit >>= 1;
							}
							testbit = 1ull << (WORDSIZE - 1);
							highest_word -= 1;
						}
						for (int i = 0; i < length; i++) result[i] = tmp[i];
					}
				}

				void expmod64x(long length, const Uint64* mod, const Uint64* base, Uint64* result, Uint64 exponent_low, Uint64 exponent_high) {
					typedef Uint64 Word;
					if (true) {
						std::vector<Word> _scratchpad; _scratchpad.resize(length * 2, 0);
						Word* t = &_scratchpad[length * 0];
						Word* u = &_scratchpad[length * 1];
						t[0] = 1; for (int i = 1; i < length; i++) t[i] = 0;
						Word bit;
						Word p = exponent_high;
						if (!p) goto starting_low_64;
						if (true) {
							bit = 1ull << 63;
							while (bit && !(bit & p)) bit >>= 1;
						}
						else {
							bit = 1ull << 63;
							bit >>= PractRand::Internals::count_high_zeroes64(p);
						}
						while (true) {
							// if (bit & exponent) t *= base;
							if (bit & p) mulmod64x(length, mod, base, t, t);
							// u = t;
							for (int i = 0; i < length; i++) u[i] = t[i];
							// t *= u;
							mulmod64x(length, mod, u, t, t);
							if (!(bit >>= 1)) break;
						}
					starting_low_64:
						bit = 1ull << 63;
						p = exponent_low;
						if (!exponent_high) {
							while (bit && !(bit & p)) bit >>= 1;
						}
						while (true) {
							// if (bit & exponent) t *= base;
							if (bit & p) mulmod64x(length, mod, base, t, t);
							// u = t;
							for (int i = 0; i < length; i++) u[i] = t[i];
							// t *= u;
							if (!(bit >>= 1)) break;
							mulmod64x(length, mod, u, t, t);

						}
						for (int i = 0; i < length; i++) result[i] = t[i];
					}
					else {
						// this version seems slower, disabled for the moment
						std::vector<Word> _scratchpad; _scratchpad.resize(length * 3, 0);
						Word* wip = &_scratchpad[length * 0];
						Word* pure = &_scratchpad[length * 1];
						Word* tmp = &_scratchpad[length * 2];
						for (int i = 0; i < length; i++) pure[i] = base[i];
						wip[0] = 1; for (int i = 1; i < length; i++) wip[i] = 0;
						while (exponent_high) {
							if (exponent_low & 1) {
								// if (exp128 & 1) wip *= pure;
								mulmod64x(length, mod, wip, pure, tmp);
								Word* old = wip; wip = tmp; tmp = old;
							}
							// exp128 >>= 1;
							exponent_low >>= 1;
							exponent_low |= exponent_high << (8 * sizeof(*mod) - 1);
							exponent_high >>= 1;

							// pure *= pure;
							mulmod64x(length, mod, pure, pure, tmp);
							Word* old = pure; pure = tmp; tmp = old;
						}
						while (exponent_low) {
							if (exponent_low & 1) {
								// if (exp128 & 1) wip *= pure;
								mulmod64x(length, mod, wip, pure, tmp);
								Word* old = wip; wip = tmp; tmp = old;
							}
							// exp128 >>= 1;
							exponent_low >>= 1;

							// pure *= pure;
							mulmod64x(length, mod, pure, pure, tmp);
							Word* old = pure; pure = tmp; tmp = old;
						}
						for (int i = 0; i < length; i++) result[i] = wip[i];
					}
				}
				void expmod32x(long length, const Uint32* mod, const Uint32* base, Uint32* result, Uint64 exponent_low, Uint64 exponent_high) {
					typedef Uint32 Word;
					if (true) {
						std::vector<Word> _scratchpad; _scratchpad.resize(length * 2, 0);
						Word* t = &_scratchpad[length * 0];
						Word* u = &_scratchpad[length * 1];
						t[0] = 1; for (int i = 1; i < length; i++) t[i] = 0;
						Uint64 bit;
						Uint64 p = exponent_high;
						if (!p) goto starting_low_64;
						if (true) {
							bit = 1ull << 63;
							while (bit && !(bit & p)) bit >>= 1;
						}
						else {
							bit = 1ull << 63;
							bit >>= PractRand::Internals::count_high_zeroes64(p);
						}
						while (true) {
							// if (bit & exponent) t *= base;
							if (bit & p) mulmod32x(length, mod, base, t, t);
							// u = t;
							for (int i = 0; i < length; i++) u[i] = t[i];
							// t *= u;
							mulmod32x(length, mod, u, t, t);
							if (!(bit >>= 1)) break;
						}
					starting_low_64:
						bit = 1ull << 63;
						p = exponent_low;
						if (!exponent_high) {
							while (bit && !(bit & p)) bit >>= 1;
						}
						while (true) {
							// if (bit & exponent) t *= base;
							if (bit & p) mulmod32x(length, mod, base, t, t);
							// u = t;
							for (int i = 0; i < length; i++) u[i] = t[i];
							// t *= u;
							if (!(bit >>= 1)) break;
							mulmod32x(length, mod, u, t, t);

						}
						for (int i = 0; i < length; i++) result[i] = t[i];
					}
					else {
						// this version seems slower, disabled for the moment
						std::vector<Word> _scratchpad; _scratchpad.resize(length * 3, 0);
						Word* wip = &_scratchpad[length * 0];
						Word* pure = &_scratchpad[length * 1];
						Word* tmp = &_scratchpad[length * 2];
						for (int i = 0; i < length; i++) pure[i] = base[i];
						wip[0] = 1; for (int i = 1; i < length; i++) wip[i] = 0;
						while (exponent_high) {
							if (exponent_low & 1) {
								// if (exp128 & 1) wip *= pure;
								mulmod32x(length, mod, wip, pure, tmp);
								Word* old = wip; wip = tmp; tmp = old;
							}
							// exp128 >>= 1;
							exponent_low >>= 1;
							exponent_low |= exponent_high << (8 * sizeof(*mod) - 1);
							exponent_high >>= 1;

							// pure *= pure;
							mulmod32x(length, mod, pure, pure, tmp);
							Word* old = pure; pure = tmp; tmp = old;
						}
						while (exponent_low) {
							if (exponent_low & 1) {
								// if (exp128 & 1) wip *= pure;
								mulmod32x(length, mod, wip, pure, tmp);
								Word* old = wip; wip = tmp; tmp = old;
							}
							// exp128 >>= 1;
							exponent_low >>= 1;

							// pure *= pure;
							mulmod32x(length, mod, pure, pure, tmp);
							Word* old = pure; pure = tmp; tmp = old;
						}
						for (int i = 0; i < length; i++) result[i] = wip[i];
					}
				}
				void expmod16x(long length, const Uint16* mod, const Uint16* base, Uint16* result, Uint64 exponent_low, Uint64 exponent_high) {
					typedef Uint16 Word;
					if (true) {
						std::vector<Word> _scratchpad; _scratchpad.resize(length * 2, 0);
						Word* t = &_scratchpad[length * 0];
						Word* u = &_scratchpad[length * 1];
						t[0] = 1; for (int i = 1; i < length; i++) t[i] = 0;
						Uint64 bit;
						Uint64 p = exponent_high;
						if (!p) goto starting_low_64;
						if (true) {
							bit = 1ull << 63;
							while (bit && !(bit & p)) bit >>= 1;
						}
						else {
							bit = 1ull << 63;
							bit >>= PractRand::Internals::count_high_zeroes64(p);
						}
						while (true) {
							// if (bit & exponent) t *= base;
							if (bit & p) mulmod16x(length, mod, base, t, t);
							// u = t;
							for (int i = 0; i < length; i++) u[i] = t[i];
							// t *= u;
							mulmod16x(length, mod, u, t, t);
							if (!(bit >>= 1)) break;
						}
					starting_low_64:
						bit = 1ull << 63;
						p = exponent_low;
						if (!exponent_high) {
							while (bit && !(bit & p)) bit >>= 1;
						}
						while (true) {
							// if (bit & exponent) t *= base;
							if (bit & p) mulmod16x(length, mod, base, t, t);
							// u = t;
							for (int i = 0; i < length; i++) u[i] = t[i];
							// t *= u;
							if (!(bit >>= 1)) break;
							mulmod16x(length, mod, u, t, t);

						}
						for (int i = 0; i < length; i++) result[i] = t[i];
					}
					else {
						// this version seems slower, disabled for the moment
						std::vector<Word> _scratchpad; _scratchpad.resize(length * 3, 0);
						Word* wip = &_scratchpad[length * 0];
						Word* pure = &_scratchpad[length * 1];
						Word* tmp = &_scratchpad[length * 2];
						for (int i = 0; i < length; i++) pure[i] = base[i];
						wip[0] = 1; for (int i = 1; i < length; i++) wip[i] = 0;
						while (exponent_high) {
							if (exponent_low & 1) {
								// if (exp128 & 1) wip *= pure;
								mulmod16x(length, mod, wip, pure, tmp);
								Word* old = wip; wip = tmp; tmp = old;
							}
							// exp128 >>= 1;
							exponent_low >>= 1;
							exponent_low |= exponent_high << (8 * sizeof(*mod) - 1);
							exponent_high >>= 1;

							// pure *= pure;
							mulmod16x(length, mod, pure, pure, tmp);
							Word* old = pure; pure = tmp; tmp = old;
						}
						while (exponent_low) {
							if (exponent_low & 1) {
								// if (exp128 & 1) wip *= pure;
								mulmod16x(length, mod, wip, pure, tmp);
								Word* old = wip; wip = tmp; tmp = old;
							}
							// exp128 >>= 1;
							exponent_low >>= 1;

							// pure *= pure;
							mulmod16x(length, mod, pure, pure, tmp);
							Word* old = pure; pure = tmp; tmp = old;
						}
						for (int i = 0; i < length; i++) result[i] = wip[i];
					}
				}
			}//BinarySignalOperations

			bool generate_xbg_characteristic_polynomial64x(long xbg_length, void (*advance_state)(Uint64*), Uint64* cp) {
				const long bits = xbg_length * 64;
				BitMatrix matrix(bits * 2 + 64, bits + 1);

				//initializing matrix state
				matrix.clear_to_zeroes();
				std::unique_ptr<Uint64[]> tmp_ptr(new Uint64[xbg_length]);
				tmp_ptr[0] = 1;
				for (long xw = 1; xw < xbg_length; xw++) tmp_ptr[xw] = 0;
				for (long y = 0; y <= bits; y++) {//XBG output on left
					for (long xw = 0; xw < xbg_length; xw ++) {
						matrix.raw64_write(xw, y, tmp_ptr[xw]);
						//matrix._import64(x, y, *(tmp_ptr++));
					}
					advance_state(&tmp_ptr[0]);
				}
				for (long x = 0; x < bits; x++) matrix.bit_set(bits + x, x);//identity matrix on right

				//gaussian elimination
				// TODO : convert to using BitMatrix::normalize_and_rank() instead of having a separate implementation just for this
				long ranks_found = 0;
				for (long x = 0; x < bits; x++) {
					bool pivot_found = false;
					long y = ranks_found;
					for (; y < bits; y++) {
						if (matrix.bit_read(x, y)) {
							if (y == ranks_found) y++;
							else matrix.xor_rows(ranks_found, y);
							pivot_found = true;
							break;
						}
					}
					if (!pivot_found) continue;
					for (; y <= bits; y++) {
						if (matrix.bit_read(x, y)) {
							matrix.xor_rows(y, ranks_found);
						}
					}
					ranks_found++;
				}

				//output
				for (long w = 0; w < xbg_length; w++) cp[w] = matrix.raw64_read(xbg_length + w, bits);
				return ranks_found == bits;
			}
			void generate_xbg_distance_code64x(long xbg_length, const Uint64 *cp, Uint64 *destination, bool backwards, Uint64 distance_low, Uint64 distance_high) {
				if (!xbg_length) return;
				if (distance_low == 0 && distance_high == 0) {
					destination[0] = 1;
					for (int i = 1; i < xbg_length; i++) destination[i] = 0;
					return;
				}
				std::vector<Uint64> _tmp;
				_tmp.resize(xbg_length);
				Uint64 *tmp = &_tmp[0];
				if (!backwards) {
					tmp[0] = 2;
					for (int i = 1; i < xbg_length; i++) tmp[i] = 0;
				}
				else {
					for (int i = xbg_length - 1; i >= 0; i--) tmp[i] = cp[i] >> 1;
					for (int i = xbg_length - 2; i >= 0; i--) tmp[i] |= cp[i + 1] << (sizeof(*cp) * 8 - 1);
					tmp[xbg_length - 1] |= 1ull << (sizeof(*cp) * 8 - 1);
				}
				using namespace BinarySignalOperations;
				expmod64x(xbg_length, cp, tmp, destination, distance_low, distance_high);
			}
			void add_xbg_distance_code64x(long xbg_length, const Uint64* cp, Uint64* destination, const Uint64* dc1, const Uint64* dc2) {
				BinarySignalOperations::mulmod64x(xbg_length, cp, dc1, dc2, destination);
			}
			void multiply_xbg_distance_code64x(long xbg_length, const Uint64* cp, Uint64* destination, const Uint64* dc, Uint64 multiplier_low, Uint64 multiplier_high) {
				BinarySignalOperations::expmod64x(xbg_length, cp, dc, destination, multiplier_low, multiplier_high);
			}
			void fastforward_xbg64x(long xbg_length, Uint64* state, void (*advance_state)(Uint64*), const Uint64* distance_code) {
				std::vector<Uint64> _tmp; _tmp.resize(xbg_length, 0);
				Uint64* tmp = &_tmp[0];

				for (long x = 0; x < xbg_length; x++) {
					Uint64 word = distance_code[x];
					for (long k = 0; k < 8*sizeof(word); k++) {
						if (word & 1) BinarySignalOperations::xor64x(xbg_length, state, tmp);
						word >>= 1;
						advance_state(state);
					}
				}
				for (long x = 0; x < xbg_length; x++) state[x] = tmp[x];
			}

			/*XorshiftMatrix::XorshiftMatrix(int size_, bool identity) {
				size = size_;
				bits.resize(size * size, false);
				if (identity) for (int i = 0; i < size; i++) set(i, i, true);
			}
			void XorshiftMatrix::apply(const std::vector<bool> &input, std::vector<bool> &output) {
				if (input.size() != size) issue_error();
				output.resize(size);
				for (int i = 0; i < size; i++) {
					bool r = false;
					for (int j = 0; j < size; j++) r ^= input[j] & get(j, i);
					output[i] = r;
				}
			}
			bool XorshiftMatrix::operator==(const XorshiftMatrix &other) const {
				return size == other.size && bits == other.bits;
			}
			XorshiftMatrix XorshiftMatrix::operator*(const XorshiftMatrix &other) const {
				XorshiftMatrix rv(size, false);
				if (other.size != size) issue_error();
				for (int i = 0; i < size; i++) {
					for (int j = 0; j < size; j++) {
						bool r = false;
						for (int k = 0; k < size; k++) {
							r ^= get(j, k) & other.get(k, i);
						}
						rv.set(j, i, r);
					}
				}
				return rv;
			}
			XorshiftMatrix XorshiftMatrix::exponent(Uint64 exponent_value) const {
				XorshiftMatrix rv(size, true), tmp(*this);
				if (!exponent_value) return rv;
				while (true) {
					if (exponent_value & 1) rv = rv * tmp;
					exponent_value >>= 1;
					if (!exponent_value) break;
					tmp = tmp * tmp;
				}
				return rv;
			}
			XorshiftMatrix XorshiftMatrix::exponent2Xminus1(Uint64 x) const {
				//exponent(2**X-1), X may be >=64
				XorshiftMatrix rv(size, true);
				while (x) {
					rv = rv * rv * *this;
					x--;
				}
				return rv;
			}
			bool XorshiftMatrix::verify_period_factorization(const std::vector<Uint64> &factors) const {
				// we test that identity is reached if we use all factors
				// we test that identity is *not* reached if we skip any single factor
				// provided that all factors are prime, this is sufficient to guarantee that cycle length is equal to the product of all factors
				XorshiftMatrix identity(size, true);
				for (long i = -1; i < (long)factors.size(); i++) {
					XorshiftMatrix tmp = *this;
					for (long j = 0; j < (long)factors.size(); j++) {
						if (i != j) tmp = tmp.exponent(factors[j]);
					}
					bool is_identity = tmp == identity;
					if (i < 0) {
						if (!is_identity) return false;
					}
					else {
						if (is_identity) return false;
					}
				}
				return true;
			}//*/


			BitMatrix::BitMatrix(long w_, long h_) {
				w = w_;
				h = h_;
				ww = (w + WORD_BITS_MASK) >> WORD_BITS_L2;
				data.resize(ww*h, 0);
			}
			void BitMatrix::clear_to_zeroes() {
				for (auto it = data.begin(); it != data.end(); it++) *it = 0;
			}
			void BitMatrix::clear_to_identity() {
				clear_to_zeroes();
				int min = w > h ? h : w;
				for (int i = 0; i < min; i++) bit_set(i, i);
			}
			bool BitMatrix::operator==(const BitMatrix &other) const {
				return w == other.w && h == other.h && data == other.data;
			}
			Uint64 BitMatrix::raw64_read(long xw, long y) {
				Uint64 rv;
				if (sizeof(rv) == sizeof(Word)) {
					rv = data[y * ww + xw];
				}
				else if (sizeof(rv) > sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(rv)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << (L2DIFF > 0 ? L2DIFF : 0) // it's just 1<<L2DIFF, the extra stuff is there to silence compiler complaints when this codepath is inactive
					};
					xw <<= L2DIFF;
					rv = data[y * ww + xw + RATIO - 1];
					for (int i = (1 << L2DIFF) - 2; i >= 0; i--) rv = (rv << (8 * sizeof(Word))) | data[y * ww + xw + i];
				}
				else if (sizeof(rv) < sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(Word)>::value - Internals::StaticLog2<sizeof(rv)>::value,
						RATIO = 1 << L2DIFF
					};
					rv = data[y*ww + (xw >> L2DIFF)] >> ((xw & (RATIO - 1)) << (3 + Internals::StaticLog2<sizeof(rv)>::value));
				}
				return rv;
			}
			Uint32 BitMatrix::raw32_read(long xw, long y) {
				Uint32 rv;
				if (sizeof(rv) == sizeof(Word)) {
					rv = data[y * ww + xw];
				}
				else if (sizeof(rv) > sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(rv)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << (L2DIFF > 0 ? L2DIFF : 0) // it's just 1<<L2DIFF, the extra stuff is there to silence compiler complaints when this codepath is inactive
					};
					xw <<= L2DIFF;
					rv = data[y * ww + xw + RATIO - 1];
					for (int i = (1 << L2DIFF) - 2; i >= 0; i--) rv = (rv << (8 * sizeof(Word))) | data[y * ww + xw + i];
				}
				else if (sizeof(rv) < sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(Word)>::value - Internals::StaticLog2<sizeof(rv)>::value,
						RATIO = 1 << L2DIFF
					};
					rv = data[y*ww + (xw >> L2DIFF)] >> ((xw & (RATIO - 1)) << (3 + Internals::StaticLog2<sizeof(rv)>::value));
				}
				return rv;
			}
			Uint16 BitMatrix::raw16_read(long xw, long y) {
				Uint16 rv;
				if (sizeof(rv) == sizeof(Word)) {
					rv = data[y * ww + xw];
				}
				else if (sizeof(rv) > sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(rv)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << (L2DIFF > 0 ? L2DIFF : 0) // it's just 1<<L2DIFF, the extra stuff is there to silence compiler complaints when this codepath is inactive
					};
					xw <<= L2DIFF;
					rv = data[y * ww + xw + RATIO - 1];
					for (int i = (1 << L2DIFF) - 2; i >= 0; i--) rv = (rv << (8 * sizeof(Word))) | data[y * ww + xw + i];
				}
				else if (sizeof(rv) < sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(Word)>::value - Internals::StaticLog2<sizeof(rv)>::value,
						RATIO = 1 << L2DIFF
					};
					rv = data[y*ww + (xw >> L2DIFF)] >> ((xw & (RATIO - 1)) << (3 + Internals::StaticLog2<sizeof(rv)>::value));
				}
				return rv;
			}
			void BitMatrix::raw64_write(long xw, long y, Uint64 value) {
				if (sizeof(value) == sizeof(Word)) {
					data[y * ww + xw] = value;
				}
				else if (sizeof(value) > sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(value)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << (L2DIFF > 0 ? L2DIFF : 0) // it's just 1<<L2DIFF, the extra stuff is there to silence compiler complaints when this codepath is inactive
					};
					xw <<= L2DIFF;
					int shift = 0;
					data[y * ww + xw + 0] = Word(value >> shift);
					for (int i = 1; i < RATIO; i++) {
						shift += 8 * sizeof(Word);
						data[y * ww + xw + i] = Word(value >> shift);
					}
				}
				else if (sizeof(value) < sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(value)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << L2DIFF,
						SHIFT_SHIFT = 3 + Internals::StaticLog2<sizeof(value)>::value
					};
					const Word mask = (Word(1) << SHIFT_SHIFT) - 1;
					int shift = (xw & (RATIO - 1)) << SHIFT_SHIFT;
					data[y * ww + (xw >> L2DIFF)] &= ~(mask << shift);
					data[y * ww + (xw >> L2DIFF)] |= Word(value) << shift;
				}
			}
			void BitMatrix::raw32_write(long xw, long y, Uint32 value) {
				if (sizeof(value) == sizeof(Word)) {
					data[y * ww + xw] = value;
				}
				else if (sizeof(value) > sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(value)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << (L2DIFF > 0 ? L2DIFF : 0) // it's just 1<<L2DIFF, the extra stuff is there to silence compiler complaints when this codepath is inactive
					};
					xw <<= L2DIFF;
					int shift = 0;
					data[y * ww + xw + 0] = Word(value >> shift);
					for (int i = 1; i < RATIO; i++) {
						shift += 8 * sizeof(Word);
						data[y * ww + xw + i] = Word(value >> shift);
					}
				}
				else if (sizeof(value) < sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(Word)>::value - Internals::StaticLog2<sizeof(value)>::value,
						RATIO = 1 << L2DIFF,
						SHIFT_SHIFT = 3 + Internals::StaticLog2<sizeof(value)>::value
					};
					const Word mask = (Word(1) << SHIFT_SHIFT) - 1;
					int shift = (xw & (RATIO - 1)) << SHIFT_SHIFT;
					data[y * ww + (xw >> L2DIFF)] &= ~(mask << shift);
					data[y * ww + (xw >> L2DIFF)] |= Word(value) << shift;
				}
			}
			void BitMatrix::raw16_write(long xw, long y, Uint16 value) {
				if (sizeof(value) == sizeof(Word)) {
					data[y * ww + xw] = value;
				}
				else if (sizeof(value) > sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(value)>::value - Internals::StaticLog2<sizeof(Word)>::value,
						RATIO = 1 << (L2DIFF > 0 ? L2DIFF : 0) // it's just 1<<L2DIFF, the extra stuff is there to silence compiler complaints when this codepath is inactive
					};
					xw <<= L2DIFF;
					int shift = 0;
					data[y * ww + xw + 0] = Word(value >> shift);
					for (int i = 1; i < RATIO; i++) {
						shift += 8 * sizeof(Word);
						data[y * ww + xw + i] = Word(value >> shift);
					}
				}
				else if (sizeof(value) < sizeof(Word)) {
					enum {
						L2DIFF = Internals::StaticLog2<sizeof(Word)>::value - Internals::StaticLog2<sizeof(value)>::value,
						RATIO = 1 << L2DIFF,
						SHIFT_SHIFT = 3 + Internals::StaticLog2<sizeof(value)>::value
					};
					const Word mask = (Word(1) << SHIFT_SHIFT) - 1;
					int shift = (xw & (RATIO - 1)) << SHIFT_SHIFT;
					data[y * ww + (xw >> L2DIFF)] &= ~(mask << shift);
					data[y * ww + (xw >> L2DIFF)] |= Word(value) << shift;
				}
			}
			/*void BitMatrix::_import64(long x, long y, Uint64 packed_bits) {
				if (sizeof(packed_bits) == sizeof(Word)) {
					data[y * ww + (x >> WORD_BITS_L2)] = Word(packed_bits);
				}
				else if (sizeof(packed_bits) < sizeof(Word)) {
					enum { RATIO = sizeof(Word) / sizeof(packed_bits) };// must be a power of 2
					enum { NUM_PACKED_BITS_L2 = sizeof(packed_bits)* 8 == 128 ? 7 : (sizeof(packed_bits)* 8 == 64 ? 6 : (sizeof(packed_bits)* 8 == 32 ? 5 : (sizeof(packed_bits)* 8 == 16 ? 4 : (sizeof(packed_bits)* 8 == 8 ? 3 : -1)))) };
					int shift = (x & (RATIO - 1)) << NUM_PACKED_BITS_L2;
					data[y * ww + (x >> WORD_BITS_L2)] |= Word(packed_bits) << shift;
				}
				else {
					data[y * ww + (x >> WORD_BITS_L2)] = Word(packed_bits);
					packed_bits >>= 8 * sizeof(Word);
					for (int i = 1; i < sizeof(packed_bits) / sizeof(Word); i++) {
						data[y * ww + (x >> WORD_BITS_L2) + 1] = Word(packed_bits);
					}
				}
			}
			void BitMatrix::_import32(long x, long y, Uint32 packed_bits) {
				if (sizeof(packed_bits) == sizeof(Word)) {
					data[y * ww + (x >> WORD_BITS_L2)] = Word(packed_bits);
				}
				else if (sizeof(packed_bits) < sizeof(Word)) {
					enum { RATIO = sizeof(Word) / sizeof(packed_bits) };// must be a power of 2
					enum { NUM_PACKED_BITS_L2 = sizeof(packed_bits)* 8 == 128 ? 7 : (sizeof(packed_bits)* 8 == 64 ? 6 : (sizeof(packed_bits)* 8 == 32 ? 5 : (sizeof(packed_bits)* 8 == 16 ? 4 : (sizeof(packed_bits)* 8 == 8 ? 3 : -1)))) };
					int shift = (x & (RATIO - 1)) << NUM_PACKED_BITS_L2;
					data[y * ww + (x >> WORD_BITS_L2)] |= Word(packed_bits) << shift;
				}
				else {
					data[y * ww + (x >> WORD_BITS_L2)] = Word(packed_bits);
					packed_bits >>= 8 * sizeof(Word);
					for (int i = 1; i < sizeof(packed_bits) / sizeof(Word); i++) {
						data[y * ww + (x >> WORD_BITS_L2) + 1] = Word(packed_bits);
					}
				}
			}
			void BitMatrix::_import16(long x, long y, Uint16 packed_bits) {
				if (sizeof(packed_bits) == sizeof(Word)) {
					data[y * ww + (x >> WORD_BITS_L2)] = Word(packed_bits);
				}
				else if (sizeof(packed_bits) < sizeof(Word)) {
					enum { RATIO = sizeof(Word) / sizeof(packed_bits) };// must be a power of 2
					enum { NUM_PACKED_BITS_L2 = sizeof(packed_bits)* 8 == 128 ? 7 : (sizeof(packed_bits)* 8 == 64 ? 6 : (sizeof(packed_bits)* 8 == 32 ? 5 : (sizeof(packed_bits)* 8 == 16 ? 4 : (sizeof(packed_bits)* 8 == 8 ? 3 : -1)))) };
					int shift = (x & (RATIO - 1)) << NUM_PACKED_BITS_L2;
					data[y * ww + (x >> WORD_BITS_L2)] |= Word(packed_bits) << shift;
				}
				else {
					data[y * ww + (x >> WORD_BITS_L2)] = Word(packed_bits);
					packed_bits >>= 8 * sizeof(Word);
					for (int i = 1; i < sizeof(packed_bits) / sizeof(Word); i++) {
						data[y * ww + (x >> WORD_BITS_L2) + 1] = Word(packed_bits);
					}
				}
			}*/
			void BitMatrix::_import_words(long offset, const Word *input, long length) {
				for (long i = 0; i < length; i++) data[offset + i] = input[i];
			}
#if 0
			void BitMatrix::import_partial_row(long x, long y, Word *input, long bits, long bit_offset, bool zeroed) {
				//49 seconds
				//added zeroed
				//46 seconds
				//added shifts
				Word *dest = &data[y*ww + (x >> WORD_BITS_L2)];
				if (false) {
					/*	//clear partial words at begining & end of region
					long end = x + bits;
					Word start_mask = (Word(1) << (x & WORD_BITS_MASK)) - 1;
					Word end_mask = ~(~Word(0) >> (((end-1) & WORD_BITS_MASK) ^ WORD_BITS_MASK));
					int end_word = end >> WORD_BITS;
					if (!zeroed) {
					clear_rectangle(x, x+bits, y, y+1);
					//dest[0] &= start_mask;
					//dest[end_word] &= end_mask;
					}
					long start_word = x >> WORD_BITS_L2;
					long end_word = (end-1) >> WORD_BITS_L2;
					long out_offset = x & WORD_BITS_MASK;
					input -= start_word;
					if (bit_offset == out_offset) {
					if (start_word == end_word) data[start_word] |= input[start_word] & ~(start_mask & end_mask);
					else {
					data[start_word] |= input[start_word] & ~start_mask;
					for (long xw = start_word+1; xw < end_word; xw++) {
					data[xw] = input[xw];
					}
					data[end_word] |= input[end_word] & ~end_mask;
					}
					}
					else {
					long bit_delta = (bit_offset - out_offset) & WORD_BITS_MASK;
					long iww =
					for (long iw = 1; iw < ; xw++) {
					;
					if (start_word == end_word) {
					}
					else {
					data[start_word] |= (input[xw] >> bit_offset) & ~start_mask;
					for (long xw = start_word+1; xw < end_word; xw++) {
					data[xw] = (input[xw] >> bit_delta) | (input[xw+1] << (WORD_BITS-bit_delta));
					}
					}
					}*/
				}
				//else if (!bit_offset && 
				else {
					if (!zeroed) clear_rectangle(x, x + bits, y, y + 1);
					for (long i = 0; i < bits; i++) {
						long a = i + bit_offset;
						bool value = (input[a >> WORD_BITS_L2] >> (a & WORD_BITS_MASK)) & 1;
						long b = i + (x & WORD_BITS_MASK);
						dest[b >> WORD_BITS_L2] |= Word(value) << (b & WORD_BITS_MASK);
					}
				}
			}
#endif
			void BitMatrix::xor_rows(long destination, long source) {
				long base1 = ww * destination;
				long base2 = ww * source;
				for (long i = 0; i < ww; i++) data[base1 + i] ^= data[base2 + i];
			}
			void BitMatrix::xor_rows_skip_start(long destination, long source, long skip) {
				long base1 = ww * destination;
				long base2 = ww * source;
				for (long i = skip; i < ww; i++) data[base1 + i] ^= data[base2 + i];
			}
			void BitMatrix::clear_rectangle(long min_x, long max_x, long min_y, long max_y) {
				long start_word = min_x >> WORD_BITS_L2;
				long end_word = (max_x - 1) >> WORD_BITS_L2;
				Word start_mask = (Word(1) << (min_x & WORD_BITS_MASK)) - 1;
				Word end_mask = ~(~Word(0) >> (((max_x - 1) & WORD_BITS_MASK) ^ WORD_BITS_MASK));
				if (start_word == end_word) {
					Word mask = start_mask | end_mask;
					for (long y = min_y; y < max_y; y++) data[y*ww + start_word] &= mask;
					return;
				}
				for (long y = min_y; y < max_y; y++) data[y*ww + start_word] &= start_mask;
				for (long x = start_word + 1; x < end_word; x++) {
					for (long y = min_y; y < max_y; y++) data[y*ww + start_word] = 0;
				}
				for (long y = min_y; y < max_y; y++) data[y*ww + end_word] &= end_mask;
			}
			long BitMatrix::normalize_and_rank() {
				long ranks_found = 0;
				for (long x = 0; x < w; x++) {
					if (bit_read(x, ranks_found)) {//row already in correct position
						for (long y = ranks_found + 1; y < h; y++) if (bit_read(x, y)) xor_rows(y, ranks_found);
						ranks_found++;
					}
					else {
						long y = ranks_found + 1;
						while (y < h && !bit_read(x, y)) y++;
						if (y != h) {
							xor_rows(ranks_found, y);
							xor_rows(y, ranks_found);
							for (++y; y < h; y++) if (bit_read(x, y)) xor_rows(y, ranks_found);
							ranks_found++;
						}
						else;//skipped rank
					}
				}
				return ranks_found;
			}
			long BitMatrix::normalize_and_rank_and_report(std::vector<Sint32> &rank_positions) {
				long ranks_found = 0;
				for (long x = 0; x < w; x++) {
					if (bit_read(x, ranks_found)) {//row already in correct position
						for (long y = ranks_found + 1; y < h; y++) if (bit_read(x, y)) xor_rows(y, ranks_found);
						rank_positions[x] = ranks_found;
						ranks_found++;
					}
					else {
						long y = ranks_found + 1;
						while (y < h && !bit_read(x, y)) y++;
						//seek_from_position_counts[x] += y - ranks_found;
						if (y != h) {
							rank_positions[x] = y;
							xor_rows(ranks_found, y);
							xor_rows(y, ranks_found);
							for (++y; y < h; y++) if (bit_read(x, y)) xor_rows(y, ranks_found);
							ranks_found++;
						}
						else {
							rank_positions[x] = -1;
						}
					}
				}
				return ranks_found;
			}
			long BitMatrix::large_normalize_and_rank() {
				long ranks_found = 0;

				//trying to minimize sweeps through memory... but not minimize actual math
				//(minimizing math costs too much)
				enum { STEP = 24 };//rows at once
				//--: 306, 4: 263, 8: 249, 16: 246, 64: 239, 256: 237
				//with "shortened" added it's now... 1: 241, 4: 205, 8: 199, 16: 197, 24: 188, 256: 183
				long last_ranks_found = 0;
				std::vector<long> rank_x;//tells us the first bit of each deferred rank in the range [last_ranks_found,ranks_found)
				long virtual_h = STEP;//rows up to (but not including) virtual_h have already been partially normalized up to ranks_found
				if (virtual_h > h) virtual_h = h;

				//Word *base = &data[0];
				for (long x = 0; x < w; x++) {
					long shortened = last_ranks_found >> WORD_BITS_L2;
					if (ranks_found >= virtual_h) {//time for a sweep
						//because too many unused ranks built up
						for (long y = virtual_h; y < h; y++) {
							for (long rank = last_ranks_found; rank < ranks_found; rank++) {
								long x2 = rank_x[rank - last_ranks_found];
								//if (read_position(x2,y)) xor_rows(y, rank);
								//if (read_position(x2,y)) for (long i = shortened; i < ww; i++) base[y*ww+i] ^= base[rank*ww+i];
								if (bit_read(x2, y)) xor_rows_skip_start(y, rank, shortened);
							}
						}
						last_ranks_found = ranks_found;
						rank_x.clear();
						virtual_h = x + STEP;
						if (virtual_h > h) virtual_h = h;
					}
					if (bit_read(x, ranks_found)) {//row already in correct position
						for (long y = ranks_found + 1; y < virtual_h; y++) if (bit_read(x, y)) xor_rows_skip_start(y, ranks_found, shortened);
						rank_x.push_back(x);
						ranks_found++;
					}
					else {
						long y = ranks_found + 1;
						while (y < virtual_h && !bit_read(x, y)) y++;
						if (y == virtual_h) {
							if (last_ranks_found != ranks_found) {//time for a sweep
								//because we need to search past the end of the stuff that is usable atm
								for (long y = virtual_h; y < h; y++) {
									for (long rank = last_ranks_found; rank < ranks_found; rank++) {
										long x2 = rank_x[rank - last_ranks_found];
										if (bit_read(x2, y)) xor_rows_skip_start(y, rank, shortened);
									}
								}
								last_ranks_found = ranks_found;
								rank_x.clear();
								virtual_h = x + STEP;
								if (virtual_h > h) virtual_h = h;
							}
							while (y < h && !bit_read(x, y)) y++;//it's now safe to search past virtual_h, because the buffer is clear
						}
						if (y != h) {
							long rfww = ranks_found * ww;
							xor_rows_skip_start(ranks_found, y, shortened);
							//for (long i = shortened; i < ww; i++) base[rfww+i] ^= base[y*ww+i];
							xor_rows_skip_start(y, ranks_found, shortened);
							//for (long i = shortened; i < ww; i++) base[y*ww+i] ^= base[rfww+i];
							for (++y; y < virtual_h; y++) if (bit_read(x, y)) xor_rows_skip_start(y, ranks_found, shortened);
							//for (++y; y < virtual_h; y++) if (read_position(x,y)) for (int i = shortened; i < ww; i++) base[y*ww+i] ^= base[ranks_found*ww+i];
							rank_x.push_back(x);
							ranks_found++;
						}
						else;//skipped rank
					}
				}
				return ranks_found;
			}
			BitMatrix BitMatrix::operator*(const BitMatrix &other) const {
				BitMatrix rv(other.w, this->h);
				if (other.h != this->w) issue_error("BitMatrix::operator* - matrix sizes incompatible");
				rv.clear_to_zeroes();

				for (int y = 0; y < rv.h; y++) {
					for (int x = 0; x < rv.w; x++) {
						bool r = false;
						for (int k = 0; k < this->w; k++) {
							r ^= bit_read(x, k) & other.bit_read(k, y);
						}
						rv.bit_write(x, y, r);
					}
				}
				return rv;
			}
			BitMatrix BitMatrix::exponentiation(Uint64 pow) const {
				if (w != h) issue_error("BitMatrix::exponentiation only supported on square matrices");
				BitMatrix rv(w, h);
				rv.clear_to_identity();
				if (!pow) return rv;
				BitMatrix tmp = *this;
				while (true) {
					if (pow & 1) rv = rv * tmp;
					pow <<= 1;
					if (!pow) break;
					tmp = tmp * tmp;
				}
				return rv;
			}
			BitMatrix BitMatrix::exponentiation128(Uint64 pow_low, Uint64 pow_high) const {
				if (w != h) issue_error("BitMatrix::exponentiation128 only supported on square matrices");
				BitMatrix rv(w, h);
				rv.clear_to_identity();
				if (!pow_low && !pow_high) return rv;
				BitMatrix tmp = *this;
				while (true) {
					if (pow_low & 1) rv = rv * tmp;
					pow_low <<= 1;
					pow_low |= pow_high << 63;
					pow_high <<= 1;
					if (!pow_low && !pow_high) break;
					tmp = tmp * tmp;
				}
				return rv;
			}
			BitMatrix BitMatrix::exponentiation2totheXminus1(Uint64 X) const {
				//exponent(2**X-1), X may be >=64
				if (w != h) issue_error("BitMatrix::exponentiation2totheXminus1 only supported on square matrices");
				BitMatrix rv(w, h);
				rv.clear_to_identity();
				while (X) {
					rv = rv * rv * *this;
					X--;
				}
				return rv;
			}
			void BitMatrix::multiply_by_vector(unsigned int length_in_bits, Uint64 *in, Uint64 *out) const {
				if (length_in_bits > w || length_in_bits > h) issue_error();
				typedef Uint64 OutWord;
				enum {
					INTERFACE_WORD_BITS = sizeof(OutWord)* 8,
					INTERFACE_WORD_BITS_L2 = 6
				};
				typedef char _compiletime_assert_bitmatrix_mulvec[((1 << INTERFACE_WORD_BITS_L2) == INTERFACE_WORD_BITS && sizeof(OutWord) == sizeof(*in)) ? 1 : -1];
				unsigned int hw = (length_in_bits + (INTERFACE_WORD_BITS - 1)) >> INTERFACE_WORD_BITS_L2;
				for (int y = 0; y < hw; y++) out[y] = 0;
				for (int y = 0; y < length_in_bits; y++) {
					bool r = false;
					for (int x = 0; x < length_in_bits; x++) r ^= ((in[x >> INTERFACE_WORD_BITS_L2] >> (x & (INTERFACE_WORD_BITS - 1))) & 1) & bit_read(x, y);
					out[y >> INTERFACE_WORD_BITS_L2] |= OutWord(r ? 1 : 0) << (y & (INTERFACE_WORD_BITS - 1));
				}
			}
			void BitMatrix::multiply_by_vector(unsigned int length_in_bits, Uint32 *in, Uint32 *out) const {
				if (length_in_bits > w || length_in_bits > h) issue_error();
				typedef Uint32 OutWord;
				enum {
					INTERFACE_WORD_BITS = sizeof(OutWord)* 8,
					INTERFACE_WORD_BITS_L2 = 5
				};
				typedef char _compiletime_assert_bitmatrix_mulvec[((1 << INTERFACE_WORD_BITS_L2) == INTERFACE_WORD_BITS && sizeof(OutWord) == sizeof(*in)) ? 1 : -1];
				unsigned int hw = (length_in_bits + (INTERFACE_WORD_BITS - 1)) >> INTERFACE_WORD_BITS_L2;
				for (int y = 0; y < hw; y++) out[y] = 0;
				for (int y = 0; y < length_in_bits; y++) {
					bool r = false;
					for (int x = 0; x < length_in_bits; x++) r ^= ((in[x >> INTERFACE_WORD_BITS_L2] >> (x & (INTERFACE_WORD_BITS - 1))) & 1) & bit_read(x, y);
					out[y >> INTERFACE_WORD_BITS_L2] |= OutWord(r ? 1 : 0) << (y & (INTERFACE_WORD_BITS - 1));
				}
			}
			void BitMatrix::multiply_by_vector(unsigned int length_in_bits, Uint16 *in, Uint16 *out) const {
				if (length_in_bits > w || length_in_bits > h) issue_error();
				typedef Uint16 OutWord;
				enum {
					INTERFACE_WORD_BITS = sizeof(OutWord)* 8,
					INTERFACE_WORD_BITS_L2 = 4
				};
				typedef char _compiletime_assert_bitmatrix_mulvec[((1 << INTERFACE_WORD_BITS_L2) == INTERFACE_WORD_BITS && sizeof(OutWord) == sizeof(*in)) ? 1 : -1];
				unsigned int hw = (length_in_bits + (INTERFACE_WORD_BITS - 1)) >> INTERFACE_WORD_BITS_L2;
				for (int y = 0; y < hw; y++) out[y] = 0;
				for (int y = 0; y < length_in_bits; y++) {
					bool r = false;
					for (int x = 0; x < length_in_bits; x++) r ^= ((in[x >> INTERFACE_WORD_BITS_L2] >> (x & (INTERFACE_WORD_BITS - 1))) & 1) & bit_read(x, y);
					out[y >> INTERFACE_WORD_BITS_L2] |= OutWord(r ? 1 : 0) << (y & (INTERFACE_WORD_BITS - 1));
				}
			}
			bool BitMatrix::verify_period_factorization(const std::vector<Uint64> &factors) const {
				// we test that identity is reached if we use all factors
				// we test that identity is *not* reached if we skip any single factor
				// provided that all factors are prime, this is sufficient to guarantee that cycle length is equal to the product of all factors
				if (w != h) issue_error("BitMatrix::verify_period_factorization - only supported on square matrices");
				BitMatrix identity(w, h);
				identity.clear_to_identity();
				for (long i = -1; i < (long)factors.size(); i++) {
					BitMatrix tmp = *this;
					for (long j = 0; j < (long)factors.size(); j++) {
						if (i != j) tmp = tmp.exponentiation(factors[j]);
					}
					bool is_identity = tmp == identity;
					if (is_identity != (i < 0)) return false;
				}
				return true;
			}
			bool BitMatrix::verify_period_factorization128(const std::vector<Uint64> &factors) const {
				//factors are stored as pairs, the low 64 bits then the high 64 bits
				if (w != h) issue_error("BitMatrix::verify_period_factorization128 - only supported on square matrices");
				BitMatrix identity(w, h);
				identity.clear_to_identity();
				for (long i = -1; i < (long)factors.size(); i++) {
					BitMatrix tmp = *this;
					for (long j = 0; j < (long)factors.size(); j++) {
						if (i != j) tmp = tmp.exponentiation128(factors[j * 2 + 0], factors[j * 2 + 1]);
					}
					bool is_identity = tmp == identity;
					if (is_identity != (i < 0)) return false;
				}
				return true;
			}

		}//XBG_helpers
	}//Internals
}//PractRand