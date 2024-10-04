#pragma once

#include "PractRand/rng_basics.h"
#include <cmath>
#include <vector>

namespace PractRand::Tests {
		//categories = old # of entries in tables
		//return value = new # of entries in tables
		//combines adjacent entries
		//N should be on the order of the sum of all counts (but, constant)
		//low N will combine many probs, high N fewer combinations
		//if aggressive is true, it will treat N as a hard limit on how low probabilities can be
		//otherwise, it will treat it as a soft limit
		//linear combines only adjacent entries; non-linear is not yet implemented
		int simplify_prob_table ( unsigned long categories, double N, double *prob_table, Uint64 *counts, bool linear, bool aggressive );
		double chi_squared_test ( unsigned long categories, const double *prob_table, const Uint64 *counts );
		double rarity_test(unsigned long categories, const double *prob_table, const Uint64 *counts);
		class G_TEST {
			long double total, sum, prob_sum, minimum_prob, partial_count, partial_prob;
			long categories;
		public:
			G_TEST() { reset(); }
			void reset();
			void set_minimum_prob(double min) { minimum_prob = min; }
			void add_category(Uint64 count, long double probability);
			void finalize();
			double get_result() const;
			long get_DoF() const;
			long get_categories() const { return categories; }
		};
		double g_test(unsigned long categories, const double *prob_table, const Uint64 *counts);
		double g_test_flat(unsigned long categories, const Uint64 *counts);
		double g_test_flat_merge_normal(unsigned long categories, const Uint64 *counts, Uint64 total = Uint64(-1), double target_ratio = 32.0);//already converted to approximately normal distribution (mandatory since DoF is not returned)
		double my_test(unsigned long categories, const double *prob_table, const Uint64 *counts);//if events are independent, this should converge to a normal distribution (mean 0 variance 1) ; intended for extremely unequal probability distributions like {0.5,0.25,0.125,0.0625,..}
		double math_chisquared_to_pvalue ( double chisquared, double DoF );
		double math_chisquared_to_normal ( double chisquared, double DoF );
		double math_pvalue_to_chisquared ( double pvalue, double DoF );
		double math_normaldist_to_pvalue(double normal);
		double math_normaldist_to_suspicion(double normal);
		double math_pvalue_to_normaldist(double pvalue);
		double math_normaldist_pdf ( double normal );
		Uint64 math_nChooseR(int set_size, int num_choices);
		double math_factorial(double a);
		double math_factorial_log(Uint64 a);//log of a!
		class SampleSet;
		//long double gap_probs( int first, int last, long double baseprob = (255.0 / 256.0) );
		//double raw_test_edge_distribution( unsigned long categories, const double *prob_table, const Uint64 *counts );
		//double test_edge_distribution( unsigned long categories, const double *prob_table, const Uint64 *counts );
		double test_uniformity( const SampleSet &sorted_data );
		double test_table_uniformity( unsigned long categories, const double *prob_table, const Uint64 *counts );

		double calculate_center_bit_combination_chance(int num_bits_L2);
		void get_hamming_weight_chances(int num_bits, std::vector<double> &pdf, std::vector<double> &cdf);//vector size = 1+(num_bits/2))
		// switches between a variety of mathods based upon the magnitude of num_bits

		Uint8  reverse_bits8 (Uint8);
		Uint16 reverse_bits16(Uint16);
		Uint32 reverse_bits32(Uint32);
		Uint64 reverse_bits64(Uint64);


		class SampleSet {
		public:
			std::vector<double> rs;
			Uint32 duplicates{};
			double sum{};
			double sum_sqr{};
			void _count_duplicates();
			double _get_index ( double other_result ) const;//interpolates
		public:
			void _normalize();
			void _add(double result) {rs.push_back(result); sum += result; sum_sqr += result * result;}
			SampleSet() = default;
			void add(const SampleSet &other) {
				for (unsigned int i = 0; i < other.rs.size(); i++) _add(other.rs[i]);
				_normalize();
			}
			void add(const double *new_results, int n) {
				long s = rs.size();
				rs.resize(s + n);
				for (int i = 0; i < n; i++) {
					rs[s+i] = new_results[i];
					sum += new_results[i];
					sum_sqr += new_results[i] * new_results[i];
				}
				_normalize();
			}
			void add(double new_result) {
				_add(new_result);
				_normalize();
			}
			void reset() {duplicates = 0; sum = 0; sum_sqr = 0; rs.resize(0);}
			long size() const {return rs.size();}
			long num_duplicates() const {return duplicates;}
			long get_index ( double other_result ) const {return int(0.5+_get_index(other_result));}
			double get_mean() const {if (rs.empty()) return 0; else return sum / rs.size();}
			double get_stddev() const {if (rs.empty()) return 0; double avg = get_mean(), avg_sqr = sum_sqr / rs.size(); return std::sqrt(avg_sqr - avg * avg);}
			double get_percentile ( double other_result ) const;//0 to 1
			double get_result_by_index(int i) const {return rs[i];}
			double get_result_by_percentile(double d) const;
			long get_num_elements_less_than ( double other_result ) const;
			long get_num_elements_greater_than ( double other_result ) const;
			void get_num_elements_less_and_greater ( double other_result, int &num_lower, int &num_higher ) const;
		};

		template<typename LowIntType, int size>
		class FixedSizeCount {
			LowIntType low[size];
			Uint64 high[size];
		public:
			int get_size() {return size;}
			void reset_counts() {
				for (int i = 0; i < size; i++) low[i] = 0;
				for (int i = 0; i < size; i++) high[i] = 0;
			}
			FixedSizeCount() {reset_counts();}
			void increment(int index) {if (!++low[index]) high[index] += 1ull << (8*sizeof(LowIntType));}
			const Uint64 &operator[] (int index) {
				high[index] += low[index];
				low[index] = 0;
				return high[index];
			}
			void flush() {
				for (int i = 0; i < size; i++) {
					high[i] += low[i];
					low[i] = 0;
				}
			}
			const Uint64 *get_array() {flush(); return &high[0];}
		};
		template<typename LowIntType>
		class VariableSizeCount {
			LowIntType *low;
			Uint64 *high;
			int size;
		public:
			int get_size() {return size;}
			void reset_counts() {
				for (int i = 0; i < size; i++) low[i] = 0;
				for (int i = 0; i < size; i++) high[i] = 0;
			}
			void set_size(int size_) {
				size = size_;
				low = static_cast<LowIntType*>(std::realloc(low, sizeof(LowIntType) * size));
				high = static_cast<Uint64*>(std::realloc(high, sizeof(Uint64) * size));
				reset_counts();
			}
			VariableSizeCount() : low(nullptr), high(nullptr), size(0) {}
			VariableSizeCount(int size_) : low(nullptr), high(nullptr), size(0) {set_size(size_);}
			VariableSizeCount(const VariableSizeCount &other) = delete;//copy constructor disallowed
			void increment(int index) {if (!++low[index]) high[index] += 1ull << (8*sizeof(LowIntType));}
			const Uint64 &operator[] (int index) {
				high[index] += low[index];
				low[index] = 0;
				return high[index];
			}
			void flush() {
				for (int i = 0; i < size; i++) {
					high[i] += low[i];
					low[i] = 0;
				}
			}
			const Uint64 *get_array() {flush(); return &high[0];}
			void swap_array(VariableSizeCount<LowIntType> &other) {
				if (other.size != size) issue_error("VariableSizeCount::swap_array");
				LowIntType *tmp_low = low;
				Uint64 *tmp_high = high;
				low = other.low;
				high = other.high;
				other.low = tmp_low;
				other.high = tmp_high;
			}
			void force_count(int index, Uint64 value) { low[index] = 0; high[index] = value; }
		};

		class BitMatrix {
			typedef Uint32 Word;
			std::vector<Word> data;
			int w, h, ww;
		public:
			static constexpr int WORD_BITS = sizeof(Word)*8;
			static constexpr int WORD_BITS_MASK = WORD_BITS-1;
			static constexpr int WORD_BITS_L2 = WORD_BITS==64?6:(WORD_BITS==32?5:(WORD_BITS==16?4:(WORD_BITS==8?3:-1)));
			void init(int w_, int h_);
			void raw_import(int offset, Word *input, int length);
			void import_partial_row(int x, int y, Word *input, int bits, int bit_offset, bool zeroed=false);
			bool read_position(int x, int y) const;
			void xor_rows(int destination, int source);
			void xor_rows_skip_start(int destination, int source, int skip);//skip is measured in words?
			void clear_rectangle(int min_x, int max_x, int min_y, int max_y);
			int normalize_and_rank();
			int large_normalize_and_rank();
		};

		struct RawTestCalibrationData_117 {
			//for use on tests that produce (very) roughly a normal distribution
			//should be based upon at least 512 samples
			const char *name;  //e.g. "Gap-16:A"
			Uint64 blocks;     //e.g. 32 for 32 KB
			Uint64 num_samples;//e.g. 65536 for that many results of known good RNGs used to construct raw_table
			Uint64 num_duplicates;

			double table[117];
			double median{};//redundant
			double mean{};
			double stddev{};

			static constexpr double ref_p[117] = {
	0.00001, 0.00002, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005,
	0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
	0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20,
	0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30,
	0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40,
	0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50,
	0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60,
	0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70,
	0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80,
	0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90,
	0.91, 0.92, 0.93, 0.94, 0.95, 0.96,	0.97, 0.98, 0.99,
	0.995, 0.998, 0.999, 0.9995, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999
};

			double get_median_sample() const {return table[49 + 9];}

			//linear interpolation
			double sample_to_index(double sample) const;
			double index_to_sample(double index) const;
			static double pvalue_to_index(double pvalue);
			static double index_to_pvalue(double index);
			double pvalue_to_sample(double pvalue) const {return index_to_sample(pvalue_to_index(pvalue));}
			double sample_to_pvalue(double sample) const {return index_to_pvalue(sample_to_index(sample));}
		};
		struct RawTestCalibrationData_129 {// 117 wasn't quite enough, or rather in a few rare cases we can do better with a little more
			//for use on tests that produce (very) roughly a normal distribution (typically chi-squared tests on overlapping samples)
			//should be based upon at least 512 test results on known good RNGs, preferably a lot more
			const char *name;  //e.g. "Gap-16:A"
			Uint64 blocks;     //e.g. 32 for 32 KB
			Uint64 num_samples;//e.g. 65536 for that many results of known good RNGs used to construct raw_table
			Uint64 num_duplicates;

			double table[129];
			double median;//redundant
			double mean;
			double stddev;

			int limit;//0 if from 129, 6 if from 117, or should it reflect the raw number of samples?

			static RawTestCalibrationData_129 *convert117to129(const RawTestCalibrationData_117 *old);

			static constexpr double ref_p[129] = {
	0.0000001, 0.0000002, 0.0000005,
	0.000001, 0.000002, 0.000005, 0.00001, 0.00002, 0.00005,
	0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005,
	0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
	0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20,
	0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30,
	0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40,
	0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50,
	0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60,
	0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70,
	0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80,
	0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90,
	0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
	0.995, 0.998, 0.999, 0.9995, 0.9998, 0.9999,
	0.99995, 0.99998, 0.99999, 0.999995, 0.999998, 0.999999,
	0.9999995, 0.9999998, 0.9999999
};

			double get_median_sample() const { return table[64]; }

			//linear interpolation
			double sample_to_index(double sample) const;
			double index_to_sample(double index) const;
			static double pvalue_to_index(double pvalue, int limit);
			static double index_to_pvalue(double index, int limit);
			double pvalue_to_sample(double pvalue) const { return index_to_sample(pvalue_to_index(pvalue, limit)); }
			double sample_to_pvalue(double sample) const { return index_to_pvalue(sample_to_index(sample), limit); }
		};
}//PractRand
