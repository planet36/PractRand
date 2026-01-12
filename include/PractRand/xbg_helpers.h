namespace PractRand {
	namespace Internals {
		namespace XBG_helpers {
			namespace BinarySignalOperations {
				void xor64x(long length, const Uint64 *k, Uint64 *value);
				void xor32x(long length, const Uint32 *k, Uint32 *value);
				void xor16x(long length, const Uint16 *k, Uint16 *value);

				void x2mod64x(long length, const Uint64 *mod, Uint64 *value);
				void x2mod32x(long length, const Uint32 *mod, Uint32 *value);
				void x2mod16x(long length, const Uint16 *mod, Uint16 *value);

				void mulmod16x(long length, const Uint16* mod, const Uint16* a, const Uint16* b, Uint16* result);
				void mulmod32x(long length, const Uint32* mod, const Uint32* a, const Uint32* b, Uint32* result);
				void mulmod64x(long length, const Uint64* mod, const Uint64* a, const Uint64* b, Uint64* result);

				//void expmod16x(long length, const Uint16* mod, const Uint16* base, Uint16* result, Uint64 exponent_low, Uint64 exponent_high);
				//void expmod32x(long length, const Uint32* mod, const Uint32* base, Uint32* result, Uint64 exponent_low, Uint64 exponent_high);
				void expmod64x(long length, const Uint64* mod, const Uint64* base, Uint64* result, Uint64 exponent_low, Uint64 exponent_high);
			}

			//returns true on success ; advance_state takes as input a pointer at an array of length xbg_length which represents an xor-based-generator (LFSR or equivalent), alters the state of that array to represent the next state in sequence
			bool generate_xbg_characteristic_polynomial64x(long xbg_length, void (*advance_state)(Uint64*), Uint64* cp);

			void generate_xbg_distance_code64x(long xbg_length, const Uint64 *cp, Uint64 *destination, bool backwards, Uint64 distance_low, Uint64 distance_high = 0);
			void add_xbg_distance_code64x(long xbg_length, const Uint64 *cp, Uint64 *destination, const Uint64 *dc1, const Uint64 *dc2);
			void multiply_xbg_distance_code64x(long xbg_length, const Uint64 *cp, Uint64 *destination, const Uint64 *dc, Uint64 multiplier_low, Uint64 multiplier_high=0);
			void fastforward_xbg64x(long xbg_length, Uint64 *state, void (*advance_state)(Uint64*), const Uint64 *distance_code);

			class BitMatrix {
				/*
					matrix math performed on binary math (aka modulo 2)
					Several uses:
						Binary Rank statistical tests
						state transition functions for Xor-Based Generators can be represented as a BitMatrix - this can be used with exponentiation to skip forward or backwards in a cycle
						deriving "characteristic polynomials" of Xor-Based Generators can be done with a BitMatrix - which can be used to skip forwards/backwards more efficiently
				*/
			public:
				typedef Uint64 Word;
			private:
				std::vector<Word> data;
				long w, h, ww;
			public:
				enum {
					WORD_BITS = sizeof(Word)* 8,
					WORD_BITS_MASK = WORD_BITS - 1,
					WORD_BITS_L2 = WORD_BITS == 64 ? 6 : (WORD_BITS == 32 ? 5 : (WORD_BITS == 16 ? 4 : (WORD_BITS == 8 ? 3 : -1)))
				};
				BitMatrix(long w_, long h_);
				void clear_to_zeroes();
				void clear_to_identity();//all zeroes except a diagonal of 1s - contrary to the name, this IS legal to use on non-square matrices
				bool operator==(const BitMatrix &other) const;
				bool operator!=(const BitMatrix &other) const { return !(*this == other); }
				BitMatrix operator*(const BitMatrix &other) const;

				bool bit_read(long x, long y) const { return (data[y*ww + (x >> WORD_BITS_L2)] >> (x & WORD_BITS_MASK)) & 1; }
				void bit_write(long x, long y, bool value) { Word &indexed = data[y*ww + (x >> WORD_BITS_L2)]; Word bit = 1; bit <<= x & WORD_BITS_MASK; Word old = indexed & ~bit; if (value) old |= bit; indexed = old; }
				void bit_clear(long x, long y) { Word &indexed = data[y*ww + (x >> WORD_BITS_L2)]; Word bit = 1; bit <<= x & WORD_BITS_MASK; indexed &= ~bit; }
				void bit_set(long x, long y) { Word &indexed = data[y*ww + (x >> WORD_BITS_L2)]; Word bit = 1; bit <<= x & WORD_BITS_MASK; indexed |= bit; }
				void bit_toggle(long x, long y) { Word &indexed = data[y*ww + (x >> WORD_BITS_L2)]; Word bit = 1; bit <<= x & WORD_BITS_MASK; indexed ^= bit; }

				Uint64 raw64_read(long xw, long y);
				Uint32 raw32_read(long xw, long y);
				Uint16 raw16_read(long xw, long y);
				void raw64_write(long xw, long y, Uint64 value);
				void raw32_write(long xw, long y, Uint32 value);
				void raw16_write(long xw, long y, Uint16 value);
				//void _import64(long x, long y, Uint64 packed_bits);// x must be a multiple of 64
				//void _import32(long x, long y, Uint32 packed_bits);// x must be a multiple of 32
				//void _import16(long x, long y, Uint16 packed_bits);// x must be a multiple of 16
				void _import_words(long offset, const Word *input, long length);
				//void import_partial_row(long x, long y, Word *input, long bits, long bit_offset, bool zeroed=false); // not worth the trouble
				void xor_rows(long destination, long source);
				void xor_rows_skip_start(long destination, long source, long skip_words);
				void clear_rectangle(long min_x, long max_x, long min_y, long max_y);// [min,max]

				long normalize_and_rank();
				long normalize_and_rank_and_report(std::vector<Sint32> &rank_positions);
				long large_normalize_and_rank();
				void multiply_by_vector(const std::vector<bool> &input, std::vector<bool> &output);
				void square_multiply_in_place(const BitMatrix &other);
				BitMatrix exponentiation(Uint64 pow) const;
				BitMatrix exponentiation128(Uint64 pow_low, Uint64 pow_high) const;
				BitMatrix exponentiation2totheXminus1(Uint64 X) const;
				void multiply_by_vector(unsigned int length_in_bits, Uint64 *in, Uint64 *out) const;
				void multiply_by_vector(unsigned int length_in_bits, Uint32 *in, Uint32 *out) const;
				void multiply_by_vector(unsigned int length_in_bits, Uint16 *in, Uint16 *out) const;

				bool verify_period_factorization(const std::vector<Uint64> &factors) const;
				bool verify_period_factorization128(const std::vector<Uint64> &factors) const;//factors are stored as pairs, the low 64 bits then the high 64 bits
			};
			/*class XorshiftMatrix {
				//matrix representing a state transition function for an RNG that uses only xors and fixed shifts
				//if the PRNG state is represented as a vector of bits, then multiplying it by its state transition matrix produces the next state
				std::vector<bool> bits;
				int size;
			public:
				XorshiftMatrix(int size_, bool identity);
				void apply(const std::vector<bool> &input, std::vector<bool> &output);
				XorshiftMatrix operator*(const XorshiftMatrix &other) const;
				bool operator==(const XorshiftMatrix &other) const;
				bool operator!=(const XorshiftMatrix &other) const { return !(*this == other); }
				XorshiftMatrix exponent(Uint64 exponent_value) const;
				XorshiftMatrix exponent2Xminus1(Uint64 X) const;
				bool verify_period_factorization(const std::vector<Uint64> &factors) const;
				bool verify_period_factorization128(const std::vector<Uint64> &factors) const;//factors are stored as pairs, the low 64 bits then the high 64 bits
				bool get(int in_index, int out_index) const { return bits[in_index + out_index * size]; }
				void set(int in_index, int out_index, bool value) { bits[in_index + out_index * size] = value; }
				void toggle(int in_index, int out_index) { bits[in_index + out_index * size] = !bits[in_index + out_index * size]; }
			};//*/
			/*class PackedBits {//used for LFSR Characteristic Polynomials only, so far
			typedef Uint64 Word;
			Word *packed_data;
			Uint64 length;
			enum {
			WORD_BITS = sizeof(Word)* 8,
			WORD_BITS_MASK = WORD_BITS - 1,
			WORD_BITS_L2 = WORD_BITS == 64 ? 6 : (WORD_BITS == 32 ? 5 : (WORD_BITS == 16 ? 4 : (WORD_BITS == 8 ? 3 : -1)))
			};
			public:
			PackedBits(Uint64 length_);//does not zero
			void zero_all();
			bool bit_read(long pos) { return (packed_data[pos >> WORD_BITS_L2] >> (pos & WORD_BITS_MASK)) & 1; }
			void bit_write(long pos, bool value) { Word &indexed = packed_data[pos >> WORD_BITS_L2]; Word bit = 1; bit <<= pos & WORD_BITS_MASK; Word old = indexed & ~bit; if (value) old |= bit; indexed = old; }
			void bit_clear(long pos) { Word &indexed = packed_data[pos >> WORD_BITS_L2]; Word bit = 1; bit <<= pos & WORD_BITS_MASK; indexed &= ~bit; }
			void bit_set(long pos) { Word &indexed = packed_data[pos >> WORD_BITS_L2]; Word bit = 1; bit <<= pos & WORD_BITS_MASK; indexed |= bit; }
			void mul_mod(PackedBits *multiplicand, PackedBits *modulus, PackedBits *result);
			void exp_mod(PackedBits *multiplicand, PackedBits *modulus, PackedBits *result, Uint64 pow_low64, Uint64 pow_high64);
			};*/

			/*
			//void xbg8_fastforward(long statelen, Uint8 *internal_state, RNGs::vRNG *rng, const Uint8 *charpoly, const Uint8 *how_far);
			BitMatrix *xbg_output_to_CP_matrix(RNGs::vRNG *rng);//uses raw64
			BitMatrix *xbg_state_to_CP_matrix(RNGs::vRNG *rng);//uses state-walking-objects
			void xbg_CP_matrix_to_CP64(const BitMatrix *matrix, Uint64 *charpoly);
			void xbg_CP_matrix_to_CP32(const BitMatrix *matrix, Uint32 *charpoly);
			void xbg_CP_matrix_to_CP16(const BitMatrix *matrix, Uint16 *charpoly);
			//void xbg_CP_matrix_to_CP8(const BitMatrix *matrix, Uint8 *charpoly);
			void xbg64_generate_distance_code(long statelen, const Uint64 *charpoly, Uint64 *distance_code);
			void xbg32_generate_distance_code(long statelen, const Uint32 *charpoly, Uint32 *distance_code);
			void xbg16_generate_distance_code(long statelen, const Uint16 *charpoly, Uint16 *distance_code);
			void xbg64_fastforward(long statelen, Uint64 *internal_state, RNGs::vRNG *rng, const Uint64 *charpoly, const Uint64 *distance_code);
			void xbg32_fastforward(long statelen, Uint32 *internal_state, RNGs::vRNG *rng, const Uint32 *charpoly, const Uint32 *distance_code);
			void xbg16_fastforward(long statelen, Uint16 *internal_state, RNGs::vRNG *rng, const Uint16 *charpoly, const Uint16 *distance_code);
			*/
		}//XBG_helpers
	}//Internals
}//PractRand