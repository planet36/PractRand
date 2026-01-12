
namespace PractRand {
	namespace Tests {
		//categories = old # of entries in tables
		//return value = new # of entries in tables
		//combines adjacent entries
		//N should be on the order of the sum of all counts (but, constant)
		//low N will combine many probs, high N fewer combinations
		//if aggressive is true, it will treat N as a hard limit on how low probabilities can be
		//otherwise, it will treat it as a soft limit
		//linear combines only adjacent entries; non-linear is not yet implemented
		int simplify_prob_table ( unsigned long categories, double N, double *prob_table, Uint64 *counts, bool ordered, bool aggressive );
		double chi_squared_test ( unsigned long categories, const double *prob_table, const Uint64 *counts );
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
		double rarity_test(unsigned long categories, const double *prob_table, const Uint64 *counts, bool ordered, bool accept_zero);//if events are independent, this should converge to a normal distribution (mean 0 variance 1) ; intended for extremely unequal probability distributions like {0.5,0.25,0.125,0.0625,0.03125,..}
		//double rarity_test(unsigned long categories, const double *prob_table, const Uint64 *counts, const double *base_scores);//base_scores are normally the logs of the probabilities from prob_table - except if the probability is zero, which IS legal, in which case the base score should be a large (but finite) negative number
		double math_chisquared_to_pvalue(double chisquared, double DoF);
		double math_chisquared_to_normal ( double chisquared, double DoF );
		double math_pvalue_to_chisquared ( double pvalue, double DoF );
		double math_normaldist_to_pvalue(double normal);
		double math_normaldist_to_suspicion(double normal);
		double math_pvalue_to_normaldist(double pvalue);
		double math_normaldist_cdf(double normal);
		double math_normaldist_pdf(double normal);
		double math_normaldist_pdf_log(double normal);
		Uint64 math_nChooseR(int set_size, int num_choices);
		inline Uint64 math_nCk(long n, long k) { return math_nChooseR(n, k);}
		double math_factorial(double a);
		double math_factorial_log(Uint64 a);//log of a!
		double math_poisson_pmf(double lambda, int value);

		class SampleSet;
		//long double gap_probs( int first, int last, long double baseprob = (255.0 / 256.0) );
		//double raw_test_edge_distribution( unsigned long categories, const double *prob_table, const Uint64 *counts );
		//double test_edge_distribution( unsigned long categories, const double *prob_table, const Uint64 *counts );
		double uniformity_test( const SampleSet &sorted_data );
		double uniformity_test( unsigned long categories, const double *prob_table, const Uint64 *counts );
		double uniformity_test_with_brute_force(unsigned long categories, const double *prob_table, const Uint64 *counts, PractRand::RNGs::Polymorphic::vRNG *known_good);
		//double uniformity_test_with_brute_force2(unsigned long categories, const double *prob_table, const Uint64 *counts, PractRand::RNGs::Polymorphic::vRNG *known_good);

		double calculate_center_bit_combination_chance(int num_bits_L2);
		void get_hamming_weight_chances(long num_bits, std::vector<double> &pdf, std::vector<double> &cdf);//vector size = 1+(num_bits/2))
		// switches between a variety of methods based upon the magnitude of num_bits

		class SampleSet {
		public:
			std::vector<double> rs;
			Uint32 duplicates;
			double sum;
			double sum_sqr;
			void _count_duplicates();
			double _get_index ( double other_result ) const;//interpolates
		public:
			void _normalize();
			void _add(double result) {rs.push_back(result); sum += result; sum_sqr += result * result;}
			SampleSet() : duplicates(0), sum(0), sum_sqr(0) {}
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
			double get_variance() const { if (rs.empty()) return 0; double avg = get_mean(), avg_sqr = sum_sqr / rs.size(); return (avg_sqr - avg * avg); }
			double get_stddev() const { if (rs.empty()) return 0; double avg = get_mean(), avg_sqr = sum_sqr / rs.size(); return std::sqrt(avg_sqr - avg * avg); }
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
		private:
			VariableSizeCount(const VariableSizeCount &other);//copy constructor disallowed
		public:
			int get_size() {return size;}
			void reset_counts() {
				for (int i = 0; i < size; i++) low[i] = 0;
				for (int i = 0; i < size; i++) high[i] = 0;
			}
			void set_size(int size_) {
				size = size_;
				if (size_) {
					low = (LowIntType*)std::realloc(low, sizeof(LowIntType)* size);
					high = (Uint64*)std::realloc(high, sizeof(Uint64)* size);
					reset_counts();
				}
				else {
					std::free(low);
					std::free(high);
					low = NULL;
					high = NULL;
				}
			}
			VariableSizeCount() : size(0), low(NULL), high(NULL) {}
			VariableSizeCount(int size_) : size(0), low(NULL), high(NULL) { set_size(size_); }
			~VariableSizeCount() { std::free(low); std::free(high); low = NULL; high = NULL; }
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

		struct RawTestCalibrationData_117 {
			//for use on tests that produce (very) roughly a normal distribution
			//should be based upon at least 512 samples
			const char *name;  //e.g. "Gap-16:A"
			Uint64 blocks;     //e.g. 32 for 32 KB
			Uint64 num_samples;//e.g. 65536 for that many results of known good RNGs used to construct raw_table
			Uint64 num_duplicates;

			double table[117];
			double median;//redundant
			double mean;
			double stddev;

			static const double ref_p[117];

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
			//static void combined129s(const RawTestCalibrationData_129 *a, const RawTestCalibrationData_129 *b, RawTestCalibrationData_129 *destination);

			static const double ref_p[129];

			double get_median_sample() const { return table[64]; }

			//linear interpolation
			double sample_to_index(double sample) const;
			double index_to_sample(double index) const;
			static double pvalue_to_index(double pvalue, int limit);
			static double index_to_pvalue(double index, int limit);
			double pvalue_to_sample(double pvalue) const { return index_to_sample(pvalue_to_index(pvalue, limit)); }
			double sample_to_pvalue(double sample) const { return index_to_pvalue(sample_to_index(sample), limit); }
		};
		/*struct RawTestCalibrationData_101 {// rethinking this slightly
			//117 wasn't enough but it was plenty in some areas and insufficient in others ; 129 helped some extremes, but only in one specific way
			//101 attempts to balance everything a little better - less granularity but more digits in the center, much higher granularity at the edges

			//for use on tests that produce (very) roughly a normal distribution with mean near 0 (typically chi-squared tests on overlapping samples forced in to something normal-ish)
			//should be based upon at least 180 test results on known good RNGs, preferably a lot more
			const char *name;  //e.g. "Gap-16:A"
			Uint64 blocks;     //e.g. 32 for 32 KB
			Uint64 num_samples;//e.g. 65536 for that many results of known good RNGs used to construct raw_table
			Uint64 num_duplicates;

			double table[101];//30 lower region, 41 central region, 30 upper region
			double minimum;
			double median;//redundant
			double maximum;
			double mean;
			double stddev;

			int limit;//0 if from 129, 6 if from 117, or should it reflect the raw number of samples?

			static RawTestCalibrationData_129 *convert117to101(const RawTestCalibrationData_117 *old);
			static RawTestCalibrationData_129 *convert129to101(const RawTestCalibrationData_117 *old);
			//static void combined129s(const RawTestCalibrationData_129 *a, const RawTestCalibrationData_129 *b, RawTestCalibrationData_129 *destination);

			static const double ref_p[101];

			double get_median_sample() const { return table[64]; }

			//linear interpolation
			double sample_to_index(double sample) const;
			double index_to_sample(double index) const;
			static double pvalue_to_index(double pvalue, int limit);
			static double index_to_pvalue(double index, int limit);
			double pvalue_to_sample(double pvalue) const { return index_to_sample(pvalue_to_index(pvalue, limit)); }
			double sample_to_pvalue(double sample) const { return index_to_pvalue(sample_to_index(sample), limit); }
		};//*/
		class TestCalibrationData {
			//const PractRand::Tests::RawTestCalibrationData_117 * raw;
			const PractRand::Tests::RawTestCalibrationData_129 * raw;
			bool owns_raw_allocation;
			double p_threshold;
			double low_sample, high_sample;
			double low_a, high_a;
			double low_b, high_b;
			double mid_sample;
		public:
			void init();
			TestCalibrationData(const PractRand::Tests::RawTestCalibrationData_117 *raw_);
			TestCalibrationData(const PractRand::Tests::RawTestCalibrationData_129 *raw_);
			~TestCalibrationData();
			double sample_to_suspicion_no_extrapolation(double sample) const;
			double sample_to_suspicion(double sample) const;
			double suspicion_to_sample(double suspicion) const;
			double sample_to_pvalue(double sample) const;
			double pvalue_to_sample(double pvalue) const;
		};
		class TestCalibrationManager {
			std::map<std::string, std::map< Uint64, TestCalibrationData * > > data;
		public:
			TestCalibrationManager(const PractRand::Tests::RawTestCalibrationData_117 *raw_data_array117, const PractRand::Tests::RawTestCalibrationData_129 *raw_data_array129);
			~TestCalibrationManager();
			Uint64 get_minimum_length(std::string name);
			TestCalibrationData *get_calibration_data(std::string name, Uint64 blocks);
		};
		//	currently there's a single central listing of calibration data
		//	I'm thinking of having separate listings for each test instead, but that's not implemented
		extern TestCalibrationManager calibration_manager;
	}//Tests
}//PractRand
