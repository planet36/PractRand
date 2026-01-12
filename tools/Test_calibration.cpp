#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <sstream>
#include <algorithm>
//#include <map>

#if defined _MSC_VER && _MSC_VER >= 1800
#include <intrin.h>
#endif

//master header, includes everything in PractRand for both 
//  practical usage and research... 
//  EXCEPT it does not include specific algorithms
#include "PractRand_full.h"

//specific RNG algorithms, to produce (pseudo-)random numbers
#include "PractRand/RNGs/all.h"
#include "PractRand/RNGs/other/mt19937.h"
#include "PractRand/RNGs/other/transform.h"
#include "PractRand/RNGs/other/lcgish.h"
#include "PractRand/RNGs/other/oddball.h"
#include "PractRand/RNGs/other/simple.h"
#include "PractRand/RNGs/other/cbuf.h"
#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/other/special.h"

//specific testing algorithms, to detect bias in supposedly random numbers
#include "PractRand/Tests/BCFN.h"
#include "PractRand/Tests/Gap16.h"
#include "PractRand/Tests/DistC6.h"
#include "PractRand/Tests/transforms.h"
#include "PractRand/Tests/FPF.h"
#include "PractRand/Tests/FPMulti.h"
#include "PractRand/Tests/BRank.h"
#include "PractRand/Tests/CoupGap.h"
#include "PractRand/Tests/mod3.h"
#include "PractRand/Tests/NearSeq.h"

using namespace PractRand;
using namespace PractRand::Internals;
using namespace PractRand::Tests;

//some helpers for the sample programs:
#include "multithreading.h"
#include "TestManager.h"
#include "MultithreadedTestManager.h"





double ref_p117[117] = {
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
	0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
	0.995, 0.998, 0.999, 0.9995, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999
};
double ref_p129[129] = {
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
double ref_p129_with_formatting[] = {
	0.0000001, 0.0000002, 0.0000005, -1,
	0.000001, 0.000002, 0.000005, -1,
	0.00001, 0.00002, 0.00005, -1,
	0.0001, 0.0002, 0.0005, -1,
	0.001, 0.002, 0.005, -1,
	0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, -1,
	0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, -1,
	0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, -1,
	0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, -1,
	0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, -1,
	0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, -1,
	0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, -1,
	0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, -1,
	0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, -1,
	0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, -1,
	0.995, 0.998, 0.999, -1,
	0.9995, 0.9998, 0.9999, -1,
	0.99995, 0.99998, 0.99999, -1,
	0.999995, 0.999998, 0.999999, -1,
	0.9999995, 0.9999998, 0.9999999, -2
};

void print_ss(const SampleSet &ss, const std::string &name, Uint64 blocks) {
//	std::printf("{\"BCFN-%d/%d\",%7.0f,%5d, %d, {", tbits, 1<<stride_L2, double(Uint64(std::pow(2,length_L2) / 1024)), (int)ss.size(), (int)ss.num_duplicates());
//	for (int i = 0; i < 117; i++) std::printf("%s%+7.3f", i ? "," : "", ss.get_result_by_percentile(ref_p[i]));
//	std::printf("}, %+.4f, %+.4f, %.4f},\n", ss.get_result_by_percentile(0.5), ss.get_mean(), ss.get_stddev());
	std::printf("{\"%s\",%9.0f,%5ld,%4ld, {", name.c_str(), (double)blocks, (long)ss.size(), (long)ss.num_duplicates());
	for (int i = 0; i < 129; i++) {
		if (ref_p129[i] >= 0.01 && ref_p129[i] <= 0.99) std::printf("%s%+7.3f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
		else std::printf("%s%+10.5f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
	}
	std::printf("}, %+.4f, %+.4f, %.4f, %d},\n", ss.get_result_by_percentile(0.5), ss.get_mean(), ss.get_stddev(), 0);
}


double fake_bcfn(PractRand::RNGs::vRNG *known_good, int tbits, Uint64 n) {
	PractRand::RNGs::LightWeight::sfc32 rng(known_good);
	int size = 1 << tbits;
	int mask = size - 1;
	Uint32 cur = rng.raw32();
	std::vector<Uint64> table; table.resize(size, 0);
	std::vector<double> probs; probs.resize(size, 1.0/size);
	n = (n + 7) >> 4;
	while (n) {
		Uint32 n32 = Uint32(n);
		if (n32 != n) n32 = 1<<30;
		n -= n32;
		for (; n32 > 0; n32--) {
			table[(cur >>  0) & mask]++;
			table[(cur >>  1) & mask]++;
			table[(cur >>  2) & mask]++;
			table[(cur >>  3) & mask]++;
			table[(cur >>  4) & mask]++;
			table[(cur >>  5) & mask]++;
			table[(cur >>  6) & mask]++;
			table[(cur >>  7) & mask]++;
			table[(cur >>  8) & mask]++;
			table[(cur >>  9) & mask]++;
			table[(cur >> 10) & mask]++;
			table[(cur >> 11) & mask]++;
			table[(cur >> 12) & mask]++;
			table[(cur >> 13) & mask]++;
			table[(cur >> 14) & mask]++;
			table[(cur >> 15) & mask]++;
			cur >>= 16;
			cur |= rng.raw32() << 16;
		}
	}
	int reduced_size = size;
	//if (!level) reduced_size = simplify_prob_table(size, samples / 32, &probs[0], &tempcount[0], true, true);
	double rv = PractRand::Tests::g_test(reduced_size, &probs[0], &table[0]);
	double rn = PractRand::Tests::math_chisquared_to_normal(rv, reduced_size-1);
	return rn;
}
double fake_bcfn2(PractRand::RNGs::vRNG *known_good, int tbits, Uint64 n, double p) {
	if (p == 0.5) return fake_bcfn(known_good, tbits, n);
	if (p <= 0 || p >= 1) issue_error();
	PractRand::RNGs::LightWeight::sfc32 rng(known_good);
	int size = 1 << tbits;
	int mask = size - 1;
	Uint32 p_i = Uint32(std::floor(p * std::pow(2.0, 32)));
	Uint32 cur = rng.raw32();
	std::vector<Uint64> table; table.resize(size, 0);
	std::vector<double> probs; probs.resize(size);
	for (Uint32 i = 0; i < size; i++) {
		int ones = count_ones32(i);
		probs[i] = std::pow(p, ones) * std::pow(1-p, tbits-ones);
	}
	n = (n + 3) >> 3;
	while (n) {
		Uint32 n32 = Uint32(n);
		if (n32 != n) n32 = 1<<30;
		n -= n32;
		for (; n32 > 0; n32--) {
			table[(cur >>  0) & mask]++;
			table[(cur >>  1) & mask]++;
			table[(cur >>  2) & mask]++;
			table[(cur >>  3) & mask]++;
			table[(cur >>  4) & mask]++;
			table[(cur >>  5) & mask]++;
			table[(cur >>  6) & mask]++;
			table[(cur >>  7) & mask]++;
#define BIT (rng.raw32() < p_i ? 1 : 0)
			Uint32 next8 = (BIT);
			next8 |= (BIT) << 1;
			next8 |= (BIT) << 2;
			next8 |= (BIT) << 3;
			next8 |= (BIT) << 4;
			next8 |= (BIT) << 5;
			next8 |= (BIT) << 6;
			next8 |= (BIT) << 7;
#undef BIT
			cur >>= 8;
			cur |= next8 << 24;
		}
	}
	int reduced_size = size;
	//if (!level) reduced_size = simplify_prob_table(size, samples / 32, &probs[0], &tempcount[0], true, true);
	double rv = PractRand::Tests::g_test(reduced_size, &probs[0], &table[0]);
	double rn = PractRand::Tests::math_chisquared_to_normal(rv, reduced_size-1);
	return rn;
}
SampleSet fake_bcfn_dist(PractRand::RNGs::vRNG *known_good, int tbits, Uint64 n, Uint32 samples, double p) {
	SampleSet ss;
	if (p == 0.5) for (Uint32 i = 0; i < samples; i++) ss._add(fake_bcfn(known_good, tbits, n));
	else for (Uint32 i = 0; i < samples; i++) ss._add(fake_bcfn2(known_good, tbits, n, p));
	ss._normalize();
	return ss;
}
void print_fake_bcfn_dist(int tbits, int stride_L2, double length_L2, int samples, bool unbalanced) {
	PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);
	SampleSet ss;
	static double chance_skipped[15] = {
		0.0,              //1 bit
		0.5,              //2 bit
		0.375,            //4 bit
		0.2734375,        //8 bit
		0.196380615234375,//16 bit
		0.139949934091419,//32 bit
		0.0993467537479669,//64 bit
		0.0703860921700151,//128 bit
		0.0498191099361402,//256 bit
		0.0352446354858388,//512 bit
		0.0249278058726663,//1 Kbit
		0.0176287723815027,//2 Kbit
		0.0124661853439194,//4 Kbit
		0.0088151932052590,//8 Kbit
		0.0062333780055594//16 Kbit
	};
	int level = stride_L2 + 3;
	double even_chance = (level <= 15) ? chance_skipped[level] : (chance_skipped[14] * std::pow(0.5, 0.5 * (level-14)));
	double p = unbalanced ? (even_chance + 1)*0.5 : 0.5;
	double unskipped_chance = unbalanced ? 1 : 1 - even_chance;
	ss = fake_bcfn_dist(&known_good, tbits, std::pow(2, length_L2 + 3 - level) * unskipped_chance - tbits + 1, samples, p);
	std::printf("{\"%s-%d/%d\",%9.0f,%5d,%4d, {", unbalanced?"BCFNU":"BCFN", tbits, 1<<stride_L2, double(Uint64(std::pow(2,length_L2) / 1024)), (int)ss.size(), (int)ss.num_duplicates());
	for (int i = 0; i < 129; i++) std::printf("%s%+7.3f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
	std::printf("}, %+.4f, %+.4f, %.4f, %d},\n", ss.get_result_by_percentile(0.5), ss.get_mean(), ss.get_stddev(), 0);
}
void blah_bcfn() {
	for (int n = 1<<10; n <= 1<<24; n<<=2) {
		for (int stride = 3; stride <= 4; stride++) for (double len = 20+stride/2; len <= 20+stride/2; len++) {
			//print_fake_bcfn_dist(6, 2,len,n, false);
			//print_fake_bcfn_dist(6, 2,len+std::log(1.5)/std::log(2.0),n, false);
			//print_fake_bcfn_dist(9, stride,len,n, true);
			print_fake_bcfn_dist(9, stride,len+std::log(1.5)/std::log(2.0),n, true);
		}
	}
}
Uint64 generate_binomial_dist(PractRand::RNGs::vRNG *known_good, Uint64 sample_length) {
	//returns number of 0s in (sample_length) random bits
	if (sample_length > 1ull << 12) {
		double p = known_good->rand_float();
		double n = Tests::math_pvalue_to_normaldist(p);
		double mean = sample_length * 0.5;
		double dev = sqrt(sample_length * 0.5 * 0.5);
		//double delta = 1.0 / dev;
		double rv = mean + n * dev;
		return Uint64(rv);
	}
	Uint32 len = sample_length;
	Uint32 rv = 0;
	for (; len >= 32; len -= 32) rv +=count_ones32(known_good->raw32());
	for (; len >= 8; len -= 8) rv += count_ones8(known_good->raw8());
	for (; len >= 1; len -= 1) rv += known_good->raw32() & 1;	
	return rv;
}
void fake_fpf_raw(PractRand::RNGs::vRNG *known_good, int tbits, Uint64 sample_length, int trials, SampleSet &ss_raw) {
	std::vector<Uint64> counts;
	long size = 1 << tbits;
	counts.resize(size);
	long mask = size - 1;
	for (long i = 0; i < trials; i++) {
		for (long x = 0; x <= mask; x++) counts[x] = 0;
		Uint64 samples_left = sample_length;
		while (samples_left >= 1024) {
			for (int x = 0; x < 256; x++) {
				counts[known_good->raw32() & mask]++; counts[known_good->raw32() & mask]++; counts[known_good->raw32() & mask]++; counts[known_good->raw32() & mask]++;
			}
			samples_left -= 1024;
		}
		long samples_left2 = samples_left;
		for (long x = 0; x < samples_left2; x++) counts[known_good->raw32() & mask]++;//*/
		/*counts[0] = sample_length;
		for (int b = 0; b < tbits; b++) {
			for (int x = 0; x < (1<<b); x++) {
				Uint64 base = counts[x];
				Uint64 half = generate_binomial_dist(known_good, base);
				counts[x] = half;
				counts[x + (1<<b)] = base-half;
			}
		}//*/

		double raw = raw = Tests::g_test_flat(mask+1,&counts[0]);
		ss_raw._add(Tests::math_chisquared_to_normal(raw, mask));
		//ss_raw._add(Tests::math_chisquared_to_pvalue(raw, mask));
		//ss_raw._add(TestResult::pvalue_to_suspicion(Tests::math_chisquared_to_pvalue(raw, mask)));
		//ss_raw._add(raw);
	}
	ss_raw._normalize();
}
void print_fake_fpf_intra(int tbits, Uint64 sample_length, int trials) {
	PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);
	SampleSet ss_raw;
	fake_fpf_raw(&known_good, tbits, sample_length, trials, ss_raw);
	//double offset = sample_length * std::log(sample_length * std::pow(0.5, tbits));
	SampleSet ss;
	//for (int i = 0; i < ss_raw.size(); i++) ss._add(Tests::math_chisquared_to_normal((ss_raw.get_result_by_index(i)+offset)*2, (1<<tbits)-1));
	//ss._normalize();
	ss = ss_raw;


	std::printf("{\"FPF-%d+6/8+\",%9.0f,%5d,%4d, {", tbits, double(sample_length), (int)ss.size(), (int)ss.num_duplicates());
	for (int i = 0; i < 129; i++) {
		if (ref_p129[i] >= 0.01 && ref_p129[i] <= 0.99) std::printf("%s%+7.3f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
		else std::printf("%s%+10.5f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
	}
	std::printf("}, %+.4f, %+.4f, %.4f, %d},\n", ss.get_result_by_percentile(0.5), ss.get_mean(), ss.get_stddev(), 0);
}
void print_fake_fpf_others(int tbits, Uint64 sample_length, int trials) {
}
void blah_fpf_all2() {
	PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);
	for (int num_platters = 1; num_platters <= 50; num_platters++) {
		SampleSet ss;
		for (int trial = 0; trial < 1 << 26; trial++) {
			double sum = 0;
			for (int p = 0; p < num_platters; p++) {
				double susp = TestResult::pvalue_to_suspicion(known_good.rand_double());
				sum += susp * susp - 4.162737902123020;
			}
			sum /= std::sqrt(double(num_platters));
			sum /= 9.308158403091918;
			ss._add(sum);
		}
		ss._normalize();
		std::printf("{\"FPF:all2\",%d,%5d,%4d, {", num_platters, (int)ss.size(), (int)ss.num_duplicates());
		for (int i = 0; i < 129; i++) {
			if (ref_p129[i] >= 0.01 && ref_p129[i] <= 0.99) std::printf("%s%+7.3f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
			else std::printf("%s%+10.5f", i ? "," : "", ss.get_result_by_percentile(ref_p129[i]));
		}
		std::printf("}, %+.4f, %+.4f, %.4f, %d},\n", ss.get_result_by_percentile(0.5), ss.get_mean(), ss.get_stddev(), 0);
	}
}
void blah_fpf() {
	for (int trials = 1 << 0; trials <= 1 << 24; trials <<= 2) {
		printf("\n\n\n");
		for (int tbits = 2; tbits <= 3; tbits+=1) {
			for (Uint64 length = 1ull<<(tbits+9); length <= 1ull<<(tbits+9); length <<= 1) {
				print_fake_fpf_intra(tbits, length, trials);
				if (length > 1) print_fake_fpf_intra(tbits, length + (length >> 1), trials);
			}
		}
	}
}

void find_test_distributions() {
	std::time_t start_time = std::time(NULL);
	std::clock_t start_clock = std::clock();

	PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);

	PractRand::RNGs::Polymorphic::efiix64x48 rng(&known_good);

	//Tests::ListOfTests tests(new Tests::BCFN2(2,13));
	//Tests::ListOfTests tests(new Tests::Gap16());
	//Tests::ListOfTests tests(new Tests::BRank(40));
	//Tests::ListOfTests tests(new Tests::BCFN_FF(2, 13));
	//Tests::ListOfTests tests(new Tests::mod3_simple());
	//Tests::ListOfTests tests(new Tests::mod3n(1));
	//Tests::ListOfTests tests(new Tests::GapUniversal2());
	//Tests::ListOfTests tests(new Tests::NearSeq3());
	//Tests::ListOfTests tests = Tests::Batteries::get_core_tests();
	//Tests::ListOfTests tests = Tests::Batteries::get_expanded_core_tests();
	Tests::ListOfTests tests(new Tests::DistC6(9,0, 1,0,0));
	//Tests::ListOfTests tests(new Tests::DistC6(6,1, 1,0,0));
	//Tests::ListOfTests tests(new Tests::DistC6(5,2, 1,0,0));
	//Tests::ListOfTests tests(new Tests::DistC6(4,3, 0,0,1));
	//Tests::ListOfTests tests(new Tests::DistC6(5,3, 1,0,1));
	//Tests::ListOfTests tests(new Tests::CoupGap());
	//Tests::ListOfTests tests(new Tests::FPF(0,14,6));
	//Tests::ListOfTests tests(new Tests::FPF(3,14,6));
	//Tests::ListOfTests tests(new Tests::FPF(4,14,6));
	//Tests::ListOfTests tests(new Tests::FPF(5,14,6));
	//Tests::ListOfTests tests(new Tests::FPF(6,14,6));
	//Tests::ListOfTests tests(new Tests::FPMulti(4, 0));
	TestManager tman(&tests, &known_good);
	tman.reset(&rng);

	//Uint64 test_size = 1 << 16;
	//test_size *= Tests::TestBlock::SIZE;

	std::map<std::string, std::map<Uint64, SampleSet> > data;
	Uint64 next_checkpoint = 1;
	enum { LARGEST_SIZE_L2 = 36 };
	enum { CHUNKY_L2 = 49 - (LARGEST_SIZE_L2 >= 34) ? LARGEST_SIZE_L2 + 4 : LARGEST_SIZE_L2 / 2 + 21 };
	enum { CHUNKY = 1 << CHUNKY_L2 };
	for (Uint64 n = 0; n <= 1ull << 30; n++) {
		if (n == next_checkpoint) {
			if (next_checkpoint < CHUNKY) next_checkpoint <<= 1; else next_checkpoint += CHUNKY;
			std::printf("\n\n\n\n");
			std::printf("==================================================\n");
			if (n & 1023) std::printf("checkpoint @ %d, time: ", int(n));
			else std::printf("checkpoint @ %dK, time: ", int(n) >> 10);
			unsigned int seconds = static_cast<unsigned long>(std::clock() - start_clock) / CLOCKS_PER_SEC;
			unsigned int days = seconds / 86400; seconds -= days * 86400;
			unsigned int hours = seconds / 3600; seconds -= hours * 3600;
			unsigned int minutes = seconds / 60; seconds -= minutes * 60;
			if (days) std::printf("%d days, ", days);
			if (days || hours) std::printf("%d hours, ", hours);
			if (days || hours || minutes) std::printf("%d minutes, ", minutes);
			std::printf("%d seconds\n", minutes);
			/*if (test_size < 10ull << 20) std::printf("for length = %d KB\n", test_size >> 10);
			else if (test_size < 10ull << 30) std::printf("for length = %d MB\n", test_size >> 20);
			else if (test_size < 10ull << 40) std::printf("for length = %d GB\n", test_size >> 30);
			else std::printf("for length = %d TB\n", test_size >> 40);*/
			std::printf("==================================================\n");
			int test_name_index = 0;
			for (std::map<std::string,std::map<Uint64,SampleSet> >::iterator it = data.begin(); it != data.end(); it++, test_name_index++) {
				std::string name = it->first;
				for (std::map<Uint64,SampleSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
					Uint64 length = it2->first;
					//int length_L2 = it2->first;
					SampleSet &ss = it2->second;
					ss._normalize();
					//std::printf("//mean= %f; median= %f; stddev= %f;\n", ss.get_mean(), ss.get_result_by_percentile(0.50), ss.get_stddev());
					std::printf("  {\"%s\",%9.0f, %d,%4d, {", name.c_str(), double(length), (int)ss.rs.size(), (int)ss.num_duplicates());
					/*double p[] = { 
						0.00001, 0.00002, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, -1, 
						0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, -1, 
						0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, -1, 
						0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, -1, 
						0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, -1, 
						0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, -1, 
						0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, -1, 
						0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, -1, 
						0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, -1, 
						0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, -1, 
						0.91, 0.92, 0.93, 0.94, 0.95, 0.96,	0.97, 0.98, 0.99, -1, 
						0.995, 0.998, 0.999, 0.9995, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999, -2
					};*/
					double *p = ref_p129_with_formatting;
					for (int i = 0; p[i] != -2; i++) {
						if (p[i] == -1) continue;
						if (i) std::printf(",");
						if (p[i] >= 0.01 && p[i] <= 0.99) std::printf("%+7.3f", ss.get_result_by_percentile(p[i]));
						else std::printf("%+10.5f", ss.get_result_by_percentile(p[i]));
					}
					std::printf("},  %+8.4f, %+8.4f, %8.4f, %d},\n", ss.get_result_by_percentile(0.50), ss.get_mean(), ss.get_stddev(), 0);
					std::fflush(stdout);
					//std::printf("}},\n");

					/*for (int L2 = 0; (2ull << L2) <= ss.rs.size(); L2++) {
						std::printf("//bott.L2:%2d:%+8.3f      ", L2, ss.get_result_by_index( (1ull << L2)-1 ));
						std::printf("topL2:%2d : %+8.3f\n",      L2, ss.get_result_by_index(ss.rs.size() - (1ull << L2)));
					}
					for (int L10 = 1; std::pow(10.0, double(L10)) <= ss.rs.size(); L10++) {
						double p = std::pow(0.1, double(L10));
						std::printf("//%f: %+8.3f      ", p, ss.get_result_by_percentile( p ));
						std::printf("%f: %+8.3f\n",   1-p, ss.get_result_by_percentile( 1-p ));
					}*/
				}
			}
		}
		Uint64 blocks_so_far = 0;
		for (int length_L2 = 10; length_L2 <= 1ull << LARGEST_SIZE_L2; length_L2 += 1) {
			if (length_L2 >= 10+3 && length_L2 < 99) {
				Uint64 new_blocks = (5ull << (length_L2-3)) / Tests::TestBlock::SIZE;
				tman.test(new_blocks - blocks_so_far);
				blocks_so_far = new_blocks;
				for (int i = 0; i < tests.tests.size(); i++) {
					std::vector<PractRand::TestResult> results;
					tests.tests[i]->get_results(results);
					for (int j = 0; j < results.size(); j++) data[results[j].name][new_blocks]._add(results[j].get_raw());
				}
			}
			if (length_L2 >= 10+2 && length_L2 <= 99) {
				Uint64 new_blocks = (3ull << (length_L2-2)) / Tests::TestBlock::SIZE;
				tman.test(new_blocks - blocks_so_far);
				blocks_so_far = new_blocks;
				for (int i = 0; i < tests.tests.size(); i++) {
					std::vector<PractRand::TestResult> results;
					tests.tests[i]->get_results(results);
					for (int j = 0; j < results.size(); j++) data[results[j].name][new_blocks]._add(results[j].get_raw());
				}
			}
			if (length_L2 >= 10+3 && length_L2 < 99) {
				Uint64 new_blocks = (7ull << (length_L2-3)) / Tests::TestBlock::SIZE;
				tman.test(new_blocks - blocks_so_far);
				blocks_so_far = new_blocks;
				for (int i = 0; i < tests.tests.size(); i++) {
					std::vector<PractRand::TestResult> results;
					tests.tests[i]->get_results(results);
					for (int j = 0; j < results.size(); j++) data[results[j].name][new_blocks]._add(results[j].get_raw());
				}
			}//*/
			if (true) {
				Uint64 new_blocks = (1ull << length_L2) / Tests::TestBlock::SIZE;
				tman.test(new_blocks - blocks_so_far);
				blocks_so_far = new_blocks;
				for (int i = 0; i < tests.tests.size(); i++) {
					std::vector<PractRand::TestResult> results;
					tests.tests[i]->get_results(results);
					for (int j = 0; j < results.size(); j++) data[results[j].name][new_blocks]._add(results[j].get_raw());
				}
			}
		}
		/*for (int i = 0; i < tests.tests.size(); i++) {
			std::vector<PractRand::TestResult> results;
			tests.tests[i]->get_results(results);
			for (int j = 0; j < results.size(); j++) data[results[j].name]._add(results[j].value);
		}*/
		tman.reset();
	}
}

static void calibrate_set_uniformity(SampleSet *calib, int n, PractRand::RNGs::vRNG *known_good) {
	for (int i = 0; i < 1ull << 27; i++) {
		SampleSet tmp;
		for (int j = 0; j < n; j++) tmp._add(known_good->rand_double());
		tmp._normalize();
		calib->_add(Tests::uniformity_test(tmp));
	}
	calib->_normalize();
}
static void compare_set_uniformity_to_table_uniformity(int n, int m, PractRand::RNGs::vRNG *known_good) {
	SampleSet _calib1, _calib2;
	SampleSet *calib1 = &_calib1, *calib2 = &_calib2;
	std::vector<Uint64> counts; counts.resize(m);
	std::vector<double> probs; probs.resize(m, 1.0 / m);
	for (int i = 0; i < 1ull << 20; i++) {
		SampleSet tmp;
		for (int j = 0; j < m; j++) counts[j] = 0;
		for (int j = 0; j < n; j++) {
			double rf = known_good->rand_double();
			tmp._add(rf);
			counts[std::floor(rf * m)]++;
		}
		tmp._normalize();
		calib1->_add(Tests::uniformity_test(tmp));
		calib2->_add(Tests::uniformity_test_with_brute_force(m, &probs[0], &counts[0], known_good));
		//calib2->_add(Tests::uniformity_test(m, &probs[0], &counts[0]));
	}
	calib1->_normalize();
	calib2->_normalize();
	std::printf("compare_set_uniformity_to_table_uniformity(%6d,%4d)", n, m);
	std::printf("    {%+6.3f,%7.3f}  ->  {%+9.3f,%7.4f}\n", calib1->get_mean(), calib1->get_stddev(), calib2->get_mean(), calib2->get_stddev());
//	std::printf("calib1: p=.001:%+5.2f    p=0.01:%+5.2f    p=0.100:%+5.2f    p=0.500:%+5.2f    p=0.900:%+5.2f    p=0.990:%+5.2f    p=0.999:%+5.2f\n",
//		calib1->get_result_by_percentile(.001), calib1->get_result_by_percentile(.01), calib1->get_result_by_percentile(.1), calib1->get_result_by_percentile(.5), calib1->get_result_by_percentile(.9), calib1->get_result_by_percentile(.99), calib1->get_result_by_percentile(.999));
//	std::printf("calib2: p=.001:%+5.2f    p=0.01:%+5.2f    p=0.100:%+5.2f    p=0.500:%+5.2f    p=0.900:%+5.2f    p=0.990:%+5.2f    p=0.999:%+5.2f\n",
//		calib2->get_result_by_percentile(.001), calib2->get_result_by_percentile(.01), calib2->get_result_by_percentile(.1), calib2->get_result_by_percentile(.5), calib2->get_result_by_percentile(.9), calib2->get_result_by_percentile(.99), calib2->get_result_by_percentile(.999));
}
static double lowest_of_N_imp1(int n, PractRand::RNGs::Polymorphic::vRNG *known_good) {
	double lowest = known_good->rand_double();
	for (int i = 1; i < n; i++) {
		double cur = known_good->rand_double();
		if (cur < lowest) lowest = cur;
	}
	return lowest;
}
static double lowest_of_N_imp2(double n, PractRand::RNGs::Polymorphic::vRNG *known_good) {
	return 1 - std::pow(known_good->rand_double(), 1.0 / n);
	//return 1 - std::pow(known_good->randlf(), 0.99 / n);
	//return std::pow(known_good->randlf(), n);
}
static void validate_lowest_of_N(PractRand::RNGs::Polymorphic::vRNG *known_good) {
	for (int n = 2; n < 10; n++) {
		SampleSet set1, set2;
		long size = 1 << 24;
		for (long i = 0; i < size; i++) set1._add(lowest_of_N_imp1(n, known_good));
		for (long i = 0; i < size; i++) set2._add(lowest_of_N_imp2(n, known_good));
		set1._normalize();
		set2._normalize();
		double step = std::sqrt(1.0 / size);
		double highest_error = 0;
		double end = 1.0 - step * 0.5;
		for (double p = step * 0.5; p < end; p += step) {
			double v1 = set1.get_result_by_percentile(p);
			double v2 = set2.get_result_by_percentile(p);
			double dif = std::fabs(v1 - v2);
			if (dif > highest_error) highest_error = dif;
		}
		highest_error /= step * std::sqrt(double(n));
		std::printf("validate_lowest_of_N(): N=%d, highest_error= %.1f units    %s\n", n, highest_error, highest_error < 2.5 ? "pass" : highest_error < 6 ? "suspicious" : "FAIL");
	}
	std::printf("validate_lowest_of_N() finished\n");
}
static void simple_chisquare_test(PractRand::RNGs::vRNG *known_good) {
	enum {SIZE = 1<<4};
	Uint64 counts[SIZE];
	double probs[SIZE];
	SampleSet ssA, ssB;
	enum {N = 8};
	for (int x = 0; x < SIZE; x++) probs[x] = 1.0 / SIZE;
	for (int i = 0; i < N; i++) {
		for (int x = 0; x < SIZE; x++) counts[x] = 0;
		for (int x = 0; x < SIZE * 100; x++) counts[known_good->rand_i32(SIZE)]++;
		double r = PractRand::Tests::g_test(SIZE, &probs[0], &counts[0]);
		double p = Tests::math_chisquared_to_pvalue(r, SIZE-1);
		double r2 = Tests::math_chisquared_to_normal(r, SIZE-1);
		double p2 = Tests::math_normaldist_to_pvalue(r2);
		ssA._add(p);
		ssB._add(p2);
	}
	ssA._normalize();
	ssB._normalize();
	SampleSet calib;
	calibrate_set_uniformity(&calib, N, known_good);
	double unifA = Tests::uniformity_test(ssA);
	double unifB = Tests::uniformity_test(ssB);
	double pA = calib.get_percentile(unifA);
	double pB = calib.get_percentile(unifB);
	std::printf("simple_chisquare_test: uniformity%.3f->percentile%.4f  uniformity%.3f->percentile%.4f\n", unifA, pA, unifB, pB);
}
void verify_test_distributions() {
	//std::time_t start_time = std::time(NULL);
	std::clock_t start_clock = std::clock();

	PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);
	PractRand::RNGs::Polymorphic::efiix32x48 rng(PractRand::SEED_AUTO);

	simple_chisquare_test(&known_good);

	Tests::ListOfTests tests = Tests::Batteries::get_core_tests();
	//Tests::ListOfTests tests = Tests::Batteries::get_expanded_core_tests();
	//Tests::ListOfTests tests(new Tests::FPF(4, 14, 6));
	//Tests::ListOfTests tests(new Tests::BRank(18));
	//Tests::ListOfTests tests(new Tests::DistC6(9,0, 1,0,0));
	//Tests::ListOfTests tests(new Tests::DistC6(6,1, 1,0,0));
	//Tests::ListOfTests tests(new Tests::DistC6(5,2, 1,0,0));
	//Tests::ListOfTests tests(new Tests::DistC6(4,3, 0,0,1));
	//Tests::ListOfTests tests(new Tests::DistC6(5,3, 1,0,1));
	//Tests::ListOfTests tests(new Tests::Gap16());
	TestManager tman(&tests, &known_good);
	tman.reset(&rng);

	std::map<std::string,std::map<int,SampleSet> > data;
	Uint64 next_checkpoint = 1;
	for (Uint64 n = 0; n <= 1<<20; n++) {
		if (n == next_checkpoint) {
			SampleSet calib;
			std::printf("\n\n\n\n");
			calibrate_set_uniformity(&calib, n, &known_good);
			std::printf("calibration finished\n");
			enum {CHUNKY = 1 << 12};
			if (next_checkpoint < CHUNKY) next_checkpoint <<= 1; else next_checkpoint += CHUNKY;
			std::printf("==================================================\n");
			std::printf("checkpoint @ %d\n", int(n) );
			std::printf("calibration = {%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f}\n", calib.get_result_by_percentile(.001), calib.get_result_by_percentile(.01), calib.get_result_by_percentile(.1), calib.get_result_by_percentile(.5), calib.get_result_by_percentile(.9), calib.get_result_by_percentile(.99), calib.get_result_by_percentile(.999));
			/*if (test_size < 10ull << 20) std::printf("for length = %d KB\n", test_size >> 10);
			else if (test_size < 10ull << 30) std::printf("for length = %d MB\n", test_size >> 20);
			else if (test_size < 10ull << 40) std::printf("for length = %d GB\n", test_size >> 30);
			else std::printf("for length = %d TB\n", test_size >> 40);*/
			std::printf("==================================================\n");
			int test_name_index = 0;
			for (std::map<std::string,std::map<int,SampleSet> >::iterator it = data.begin(); it != data.end(); it++, test_name_index++) {
				std::string name = it->first;
				for (std::map<int,SampleSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
					int length_L2 = it2->first;
					SampleSet &ss = it2->second;
					ss._normalize();
					//if (ss.num_duplicates()) continue;
					std::printf("\n\n name=\"%s\"; length_L2=%d;\n", name.c_str(), length_L2);
					std::printf("total= %lld; duplicates= %ld;\n", Uint64(ss.rs.size()), ss.num_duplicates());
					if (!ss.num_duplicates()) {
						double sum = Tests::uniformity_test(ss);
						double p = calib.get_percentile(sum);
						std::printf("p = %9.7f      raw uniformity = %+6.3f\n", p, sum);
						if (fabs(p - .5) >= 0.49999) std::printf("!!!!!!!!!!!!!!!!!!");
						if (fabs(p - .5) >= 0.4999) std::printf("!!!!!!!!!!!!!!!!!!");
						if (fabs(p - .5) >= 0.499) std::printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
					}
					std::printf("\n");
				}
			}
		}
		Uint64 blocks_so_far = 0;
		for (int length_L2 = 21; length_L2 <= 30; length_L2 += 1) {
			//data["fred"][length_L2]._add(known_good.randf(0.95));continue;
			Uint64 new_blocks = (1ull << length_L2) / Tests::TestBlock::SIZE;
			tman.test(new_blocks - blocks_so_far);
			blocks_so_far = new_blocks;
			for (int i = 0; i < tests.tests.size(); i++) {
				std::vector<PractRand::TestResult> results;
				tests.tests[i]->get_results(results);
				for (int j = 0; j < results.size(); j++) data[results[j].name][length_L2]._add(results[j].get_pvalue());
			}
		}
		tman.reset();
	}
}

struct Data {
	Uint64 count;
	double vec[129];
	Data() : count(0) {
		for (int i = 0; i < 129; i++) vec[i] = 0;
	}
};
std::map<std::string,std::map<int, Data > > data;

void print_data() {
	for (std::map<std::string,std::map<int,Data > >::iterator it = data.begin(); it != data.end(); it++) {
		std::string name = it->first;
		for (std::map<int,Data >::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			int length_L2 = it2->first;
			const Data &d = it2->second;
			std::printf("{ \"%s\", 1ull << (%d - 10), %.0f, {", name.c_str(), length_L2, (double)d.count );
			for (int i = 0; i < 129; i++) {
				if (ref_p129[i] >= 0.01 && ref_p129[i] <= 0.99) std::printf("%+7.3f", d.vec[i]);
				else std::printf("%+10.5f", d.vec[i]);
				if (i != 128) std::printf(",");
			}
			std::printf("}},\n");
		}
	}
}

static int sum_of_bytes(Uint32 in) {
	const Uint32 mask1 = 0x00FF00FF;
	in = (in & mask1) + ((in>>8) & mask1);
	in += in >> 16;
	return (in & 1023);
}

//#define GENERATE_GAUSSIAN() math_pvalue_to_normaldist(rng.rand_double()) //133
#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_lookup3_linear(rng.raw64()) //277
//#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_high_quality(rng.raw64(), rng.raw64(), rng.raw64())
//#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_popcount_lookup_linear(rng.raw64(), rng.raw64()) //248
//#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_popcount_lookup_quadratic(rng.raw64(), rng.raw64()) //233
//#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_popcount_lookup_cubic(rng.raw64(), rng.raw64()) //222
//#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_popcount_sum(rng.raw64(), rng.raw64()) //228
//#define GENERATE_GAUSSIAN() PractRand::Internals::generate_gaussian_popcount_sum2(rng.raw64(), rng.raw64(), rng.raw64(), rng.raw64(), rng.raw64()) //114
//#define GENERATE_GAUSSIAN() generate_with_ziggurat(rng) //173

void test_normal_distribution_a() {
	//PractRand::RNGs::LightWeight::sfc64 rng( PractRand::SEED_AUTO );
	//PractRand::RNGs::Polymorphic::sfc64 rng( PractRand::SEED_AUTO );
	//PractRand::RNGs::Polymorphic::mrsf64 rng( PractRand::SEED_AUTO );
	PractRand::RNGs::LightWeight::mrsf64 rng(PractRand::SEED_AUTO);
	GENERATE_GAUSSIAN();
	std::clock_t bench_start, bench_end;
	bench_start = std::clock();
	while (bench_start == (bench_end = std::clock())) ;
	bench_start = bench_end;
	Uint32 count = 0;
	while (std::clock_t(CLOCKS_PER_SEC*1.25 + 1) > std::clock_t((bench_end = std::clock()) - bench_start)) {
		double a = GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN() + GENERATE_GAUSSIAN();
		count += a != -123456789;
	}
	double rate = count*8.0 / (double(std::clock_t(bench_end-bench_start))/CLOCKS_PER_SEC);
	std::printf("gaussian speed: %.3f M / second\n", rate / 1000000.0 );
	if (true) {
		SampleSet ss; ss.rs.reserve(1<<27);
		for (int n = 8; n <= 28; n+=4) {
			for (long i = (1 << n) - ss.size(); i > 0; i--) ss._add(GENERATE_GAUSSIAN());
			ss._normalize();
			printf("test_normal_distribution_a0 : 2^%02d: mean:%+.5f, stddev:%.5f", n, ss.get_mean(), ss.get_stddev());
			double cumulative_err = 0;
			double cumulative_err2 = 0;
			double highest_err = 0;
			double lowest = 1.0, highest = 0.0;
			for (long i = 0; i < ss.size(); i++) {
				double expected_p = (i + 0.5) / ss.size();
				double raw = ss.get_result_by_index(i);
				double observed_p = PractRand::Tests::math_normaldist_to_pvalue(raw);
				double err = std::fabs(expected_p - observed_p);
				if (observed_p < lowest) {
					lowest = observed_p;
				}
				if (observed_p > highest) highest = observed_p;
				if (err > highest_err) highest_err = err;
				cumulative_err += err;
				cumulative_err2 += err * err;
			}
			printf(" err: %.5f  err2: %.5f  e.high:%.5f  extr: %3.1f  ex.ratio: %.2f\n", cumulative_err / ss.size(), std::sqrt(cumulative_err2 / ss.size()), highest_err, (lowest + (1-highest)) * (ss.size()-1.0), lowest / (lowest + (1.0 - highest)));
		}
	}
	if (true) {
		SampleSet ss; ss.rs.reserve(1<<27);
		for (int n = 8; n <= 24; n+=4) {
			for (long i = (1 << n) - ss.size(); i > 0; i--) ss._add(Tests::math_normaldist_to_pvalue(GENERATE_GAUSSIAN()));
			ss._normalize();
			printf("test_normal_distribution_a1 : 2^%02d: uniformity %+f\n", n, Tests::uniformity_test(ss));
		}
	}
	if (true) {
		enum {TBITS=17};
		std::vector<Uint64> counts; counts.resize(1<<TBITS, 0);
		Uint64 total = 0;
		const double p_thresh_low = 0.0;
		const double p_thresh_high = 1.0;
		double n_thresh_low = -999999999999999.0; if (p_thresh_low < 1) n_thresh_low = Tests::math_pvalue_to_normaldist(p_thresh_low);
		double n_thresh_high = 999999999999999.0; if (p_thresh_high < 1) n_thresh_high = Tests::math_pvalue_to_normaldist(p_thresh_high);
		const double scale = std::pow(2.0, TBITS) / (p_thresh_high - p_thresh_low);
		for (int n = 16; n <= 40; n++) {
			while (!(total >> n)) {
				for (long i = 0; i < 1<<16; ) {
					double norm = GENERATE_GAUSSIAN();
					if (norm < n_thresh_low) continue;
					if (norm > n_thresh_high) continue;
					double p = Tests::math_normaldist_to_pvalue(norm);
					if (p > p_thresh_high || p < p_thresh_low) PractRand::issue_error();
					double a = p - p_thresh_low;
					double index_f = a * scale;
					long index = long(std::floor(index_f));
					if (index & ~((1L << TBITS) - 1)) PractRand::issue_error();
					counts[index]++;
					i++;
				}
				total += 1<<16;
			}
			printf("test_normal_distribution_a2 : 2^%02d: %+f\n", n, Tests::g_test_flat_merge_normal(1<<TBITS, &counts[0], 1ull<<n) );
		}
	}
}

Uint64 count_period(PractRand::RNGs::vRNG *rng) {
	enum {BYTES = 32};//must be a power of 2 greater than or equal to 8
	if (rng->get_native_output_size() == 8) {
		typedef Uint8 Word;
		enum { BUFSIZE = BYTES / sizeof(Word) };
		Word buff1[BUFSIZE];
		Word buff2[BUFSIZE];
		rng->autoseed();
		for (int i = 0; i < BUFSIZE; i++) buff1[i] = rng->raw8();
		for (Uint64 p = 0; p < (1ull << 48); p++) {
			if ((buff2[p & (BUFSIZE - 1)] = rng->raw8()) == buff1[BUFSIZE - 1]) {
				bool match = true;
				for (int i = 1; i < BUFSIZE; i++) if (buff2[(p - i) & (BUFSIZE - 1)] != buff1[BUFSIZE - 1 - i]) match = false;
				if (match) {
					Uint64 rv = p + 1;
					std::printf("cycle found in [%s] of length %.0f\n", rng->get_name().c_str(), double(rv));
					if (!(rv & ((1 << 20)-1))) std::printf("%.0fM exactly\n", double(rv >> 20));
					else if (!(rv & 1023)) std::printf("%.0fK exactly\n", double(rv >> 10));
					return rv;
				}
			}
		}
	}
	return -1;
}

void test_sfc16() {
	if (1) {//search for short cycles
		enum {
			BUFFER_SIZE = 16,//plenty to tell if the state matched (6 would probably be enough, but more doesn't hurt)
			SHORT_CYCLE_L2 = 40,//the log-based-2 of the cycle length we consider "too short"
		};
		Uint16 buffy[BUFFER_SIZE];
		PractRand::RNGs::Raw::sfc16 rng;
		for (Uint64 seed = 0; true; seed++) {
			rng.seed(seed);
			for (int i = 0; i < BUFFER_SIZE; i++) buffy[i] = rng.raw16();
			for (Uint64 i = 0; i < ((1ull << SHORT_CYCLE_L2) - BUFFER_SIZE); i++) rng.raw16();
			int match = true;
			for (int i = 0; i < BUFFER_SIZE; i++) if (buffy[i] != rng.raw16()) match = false;
			if (match || !(seed & 15)) {
				std::printf("seed 0x");
				if (seed >> 32) {
					std::printf("%X", Uint32(seed >> 32));
					std::printf("%08X", Uint32(seed >> 0));
				}
				else std::printf("%X", Uint32(seed));
				if (match) std::printf(" had a short cycle\n\n  !!!!!!!!!!!!!!!!!!!\n\n  !!!!!!!!!!!!!!!!!!!\n\n  !!!!!!!!!!!!!!!!!!!\n");
				else std::printf(" was good\n");
			}
		}

	}
	if (0) {
		long double p248 = std::pow(2, 48);
		long double base = std::pow(0.5, 48);
		long double base_inv = 1 - base;
		long double after2_16 = 1;
		for (int i = 1; i < 65536; i++) {
			after2_16 *= 1 - (1 / (p248 - i));
		}
		std::printf("chance of a seed leading to a cycle < 2**32: %g\n", double(1 / (1 - after2_16)));
		long double after2_24 = 1;
		for (int i = 1; i < 65536 * 256; i++) {
			after2_24 *= 1 - (1 / (p248 - i));
		}
		std::printf("chance of a seed leading to a cycle < 2**40: %g\n", double(1 / (1 - after2_24)));
		long double after2_32 = 1;
		for (unsigned long long i = 1; i < 65536 * 65536ull; i++) {
			after2_32 *= 1 - (1 / (p248 - i));
		}
		std::printf("chance of a seed leading to a cycle < 2**48: %f\n", double(1 / (1 - after2_32)));
	}
}


void compare_bernoulli_distribution_calculations() {
	long long num_bits_L2 = 38;
	long long num_bits = 1ull << num_bits_L2;
	long long half_bits = num_bits >> 1;
	long long region_boundary = half_bits - (1 << ((num_bits_L2 >> 1) - 1));
	long long exponent = -num_bits;
	long double p = 1.0;
	long double prev_cdf = 0;

	long double mean = num_bits / 2.0;
	long double dev = std::sqrt(num_bits * 0.5 * 0.5);
	long double delta = 1.0 / dev;

	long double prev_cdf2 = 0;
	long double prev_cdf3 = 0;
	long double pdfdiff[2] = { 0, 0 };
	long double cdfdiff[2] = { 0, 0 };

	for (long long i = 0; i <= half_bits; i++) {
		if (i) {
			p /= i; p *= num_bits + 1 - i;
		}
		int tmp_exp; p = std::frexpl(p, &tmp_exp); exponent += tmp_exp;
		long double tmp_p = 0; if (exponent > -1070) tmp_p = std::ldexpl(p, exponent);
		long double new_cdf = prev_cdf + tmp_p;

		double norm = (i - mean) * delta;
		long double cdf2 = math_normaldist_to_pvalue(norm + 0.5 * delta);
		long double pdf2 = math_normaldist_pdf(norm) * delta;
		//long double pdf3 = (math_normaldist_pdf(norm + 0.25 * delta) + math_normaldist_pdf(norm - 0.25 * delta) - 1 * math_normaldist_pdf(norm)) * delta / 1;
		long double cdf3 = prev_cdf3 + pdf2;
		//cdf[i] = math_normaldist_to_pvalue(norm + 0.5 * delta);
		//pdf[i] = math_normaldist_pdf(norm) * delta;

		long ri = i < region_boundary ? 0 : 1;
		pdfdiff[ri] += std::fabs(tmp_p - pdf2) * (i == half_bits ? 1 : 2);
		cdfdiff[ri] += std::fabs(new_cdf - cdf2) * (i == half_bits ? 1 : 2);
		if (i == half_bits) {
			std::printf("pre-final CDFs:\n\t0: %.15f\n\t1: %.15f\n\t2: %.15f\n\t~: %.15f\n", double(prev_cdf), double(prev_cdf2), double(prev_cdf3), 0.5 - 0.5 * calculate_center_bit_combination_chance(num_bits_L2));
		}
		prev_cdf = new_cdf;
		prev_cdf2 = cdf2;
		prev_cdf3 = cdf3;
	}
	std::printf("final CDFs:\n\t0: %.15f\n\t1: %.15f\n\t2: %.15f\n\t~: %.15f\n", double(prev_cdf), double(prev_cdf2), double(prev_cdf3), 0.5 + 0.5 * calculate_center_bit_combination_chance(num_bits_L2));
	std::printf("\n");
	std::printf("region 0: PDF delta: %g\n", double(pdfdiff[0]));
	std::printf("region 1: PDF delta: %g\n", double(pdfdiff[1]));
	std::printf("combined: PDF delta: %g\n", double(pdfdiff[0] + pdfdiff[1]));
	std::printf("\n");
	std::printf("region 0: CDF delta: %g\n", double(cdfdiff[0]));
	std::printf("region 1: CDF delta: %g\n", double(cdfdiff[1]));
	std::printf("combined: CDF delta: %g\n", double(cdfdiff[0] + cdfdiff[1]));
}

void count_hash_deviation() {
	Uint64 counts[65536];
	long lowest_uniques = 99999;
	long highest_uniques = 0;
	long sum_uniques = 0;
	for (long i = 0; i < 65536; i++) counts[i] = 0;
	for (long y = 0; y < 65536; y++) {
		Uint64 counts2[65536] = { 0 };
		long uniques = 0;
		for (long x = 0; x < 65536; x++) {
			Uint16 hashed = x + PractRand::Internals::rotate16(x+y*0, 8);
			counts[hashed]++;
			if (!counts2[hashed]++) uniques++;
		}
		if (uniques < lowest_uniques) lowest_uniques = uniques;
		if (uniques > highest_uniques) highest_uniques = uniques;
		sum_uniques += uniques;
	}
	long lowest = 999999;
	long highest = 0;
	Sint64 total_err = 0;
	for (long i = 0; i < 65536; i++) {
		Sint64 c = counts[i];
		if (c < lowest) lowest = c;
		if (c > highest) highest = c;
		total_err += std::abs(c - 65536);
	}
	std::printf("count_hash_deviation:\n\tlowest:  %ld\n\thighest: %ld\n\taverage error: %.5f\n\texpected: 65536.0\n", lowest, highest, total_err / 65536.0);
	std::printf("\n\tuniques: lowest:  %ld\n\taverage: %.5f\n\thighest: %ld\n", lowest_uniques, sum_uniques / 65536.0, highest_uniques);
}

#if 0
static const Uint8 count_high_zeroes_table[256] = {
	//	0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
	8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,//0
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,//1
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,//2
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,//3
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,//4
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,//5
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,//6
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,//7
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//8
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//9
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//10
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//11
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//12
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//13
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//14
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//15
};
unsigned long check_count_high_zeroes64(Uint64 value) {
	Uint8 count = 0;
	if (!(value >> 32)) { count += 32; value <<= 32; }
	if (!(value >> 48)) { count += 16; value <<= 16; }
	if (!(value >> 56)) { count += 8; value <<= 8; }
	return count_high_zeroes_table[value >> 56] + count;
}

void check_optimization() {
	Uint64 lcg_state[4] = { 0, 0, 0, 0 };
	Uint64 lcg_adder = 0x311A6E831C5ull;
	for (Uint64 i = 0; i < 100000; i++) {
		//just wanting to see the disassembly of this
		Uint8 carry = 0;
		check_count_high_zeroes64(i);
		//carry = _addcarry_u64(carry, lcg_state[0], lcg_adder, &lcg_state[0]);
		//carry = _addcarry_u64(carry, lcg_state[1], lcg_state[0], &lcg_state[1]);
		//carry = _addcarry_u64(carry, lcg_state[2], lcg_state[1], &lcg_state[2]);
		//carry = _addcarry_u64(carry, lcg_state[3], lcg_state[2], &lcg_state[3]);
	}
	if (lcg_state[0] == 0) std::printf("");
}
#endif

#pragma intrinsic(_umul128)
Uint64 _mulmod(Uint64 a, Uint64 b, Uint64 m) {
	if (m <= 1) issue_error();
	if (m >> 63) issue_error();//not confident this produces correct results when m is that close to the maximum
	if (a >= m) a %= m;
	if (b >= m) b %= m;
	Uint64 product = 0, tmp = a;
	while (b) {
		if (b & 1) {
			product += tmp;
			if (product >= m) product -= m;
		}
		b >>= 1;
		tmp <<= 1;
		if (tmp >= m) tmp -= m;
	}
	return product;
}
Uint32 _expmod32(Uint32 _x, Uint32 _e, Uint32 _m) {
	if (_m <= 1) issue_error();
	if (_e >= _m) issue_error();
	Uint64 x = _x, e = _e, m = _m;
	if (x >= m) x %= m;
	Uint64 result = 1, tmp = x;
	while (e) {
		if (e & 1) result = (result * tmp) % m;
		e >>= 1;
		tmp = (tmp * tmp) % m;
	}
	return result;
}
Uint64 _expmod64(Uint64 x, Uint64 e, Uint64 m) {
	if (m <= 1) issue_error();
	if (e >= m) issue_error();
	if (m >> 63) issue_error();//not confident this produces correct results when m is that close to the maximum
	if (x >= m) x %= m;
	Uint64 result = 1, tmp = x;
	while (e) {
		if (e & 1) result = _mulmod(result, tmp, m);
		e >>= 1;
		tmp = _mulmod(tmp, tmp, m);
	}
	return result;
}
void check_np2lcg_optimization(const Uint64 modulus, const Uint64 multiplier, const Uint64 *factors, const Uint64 num_factors) {
	if (double(multiplier) * double(modulus) >= 18446744073709551616.0) issue_error("overflow");
	enum { WARMUP = 200 };
	const Uint64 shift = std::ceil(std::log2(double(modulus)));
	const Uint64 post_shift_multiplier = (1ull << shift) % modulus;
	const Uint64 period = modulus - 1;
	const Uint64 mask = (1ull << shift) - 1;
	const Uint64 hard_max = mask + (multiplier - 1) * post_shift_multiplier;
	const Uint64 max_loops = 1ull << 42;
	Uint64 lowest_post_decrease = modulus, highest_post_increase = 1;
	Uint64 current_value = 1;
	Uint64 warmup = WARMUP;
	Uint64 loops = modulus + WARMUP;
	if (loops > max_loops) loops = max_loops;
	/*
		The classical modulus operation is a bit slow.  
		So instead, I use a modulus value just under a power of 2, and adjust based upon a shift and multiply (though I choose my numbers such that the multiplication can be done with an LEA in most cases).  
		This nominally works, but the result isn't strictly normalized.  That's fixable simply by loosening my definition of normalized... mostly.  
		Unfortunately, playing fast and loose with normalization causes things to skip around in ways that may cause problems.  
		But if I do the process twice, that's still faster than modulus, and the result seems to be perfect
		...although the old zero was at one end of the range, but is now in the middle (well, near-ish to the end, but not actually at the end) - that tiny hiccup doesn't seem to hurt much though.  
	*/
	for (Uint64 i = loops; i > 0; i--) {
		Uint64 next = current_value * multiplier;
		next = (next & mask) + (next >> shift) * post_shift_multiplier;
		next = (next & mask) + (next >> shift) * post_shift_multiplier;
		//next = (next & mask) + (next >> shift) * post_shift_multiplier;
		Uint64 prev = current_value;
		current_value = next;
		if (warmup) {warmup--;  continue;}
		/*if (next > mask) issue_error("next > mask");
		if (next <= post_shift_multiplier) {
			std::printf("prev_value: %lld   next_value: %lld\n", prev, next);
			issue_error("next < post_shift_multiplier");
		}*/
		//if (next < mask - modulus) issue_error("next < mask-modulus");
		if (next < prev) {//decrease
			if (next < lowest_post_decrease) {
				lowest_post_decrease = next;
			}
		}
		else if (next > prev) {//increase
			if (next > highest_post_increase) {
				highest_post_increase = next;
			}
		}
		else {//skip this
			issue_error("cycle length 1?");
		}
	}
	Uint64 spread = highest_post_increase - lowest_post_decrease + 1;
	std::printf("check_np2lcg_optimization result: spread = %lld, modulus = %lld, difference= %+lld, post-shift-multiplier: %lld, time: %.1f\n", spread, modulus, spread - modulus, post_shift_multiplier, std::clock() * (1.0 / CLOCKS_PER_SEC));
}
void check_np2lcg_multiplier() {
	// everyone wants non-power-of-2-modulus LCGs to have maximal cycle lengths
	// even I do, because without the benefits that grants (knowing the cycle length, knowing which cycle any given value is on, etc) I wouldn't bother even contemplating NP2LCGs
	// despite the fact that, in my experience, maximal cycle length multipliers actually yield lower quality output, on average, than other multipliers
	// and the fact that LCGs in general kind of suck, especially NP2LCGs

	// this code, when given a prime modulus M and the prime factorization of (M-1), will search for maximum length cycles

	//*
	Uint64 modulus = (1ull << 10) - 3;
	Uint64 factors[] = { 2, 2, 3, 5, 17, 0 };
	Uint64 chosen_multiplier = 177;//*/
	/*
	Uint64 modulus = (1ull << 11) - 9;
	Uint64 factors[] = { 2, 1019, 0 };
	Uint64 chosen_multiplier = 913;//*/
	/*
	Uint64 modulus = (1ull << 13) - 1;
	Uint64 factors[] = { 2, 3, 3, 5, 7, 13, 0 };
	Uint64 chosen_multiplier = 2185;//*/
	/*
	Uint64 modulus = (1ull << 17) - 1;
	Uint64 factors[] = { 2, 3, 5, 17, 257, 0 };
	Uint64 chosen_multiplier = 8485;//*/
	/*
	Uint64 modulus = (1ull << 19) - 1;
	Uint64 factors[] = { 2, 3, 3, 3, 7, 19, 73, 0 };
	Uint64 chosen_multiplier = 5601;//*/
	/*
	Uint64 modulus = (1ull << 22) - 3;
	Uint64 factors[] = { 2, 2, 3, 5, 5, 11, 31, 41, 0 };
	Uint64 chosen_multiplier = 916;//*/
	/*
	Uint64 modulus = (1ull << 24) - 3;
	Uint64 factors[] = { 2, 2, 3, 23, 89, 683, 0 };
	Uint64 chosen_multiplier = 17880;//*/
	/*
	Uint64 modulus = (1ull << 26) - 5;
	Uint64 factors[] = { 2, 479, 70051, 0 };
	Uint64 chosen_multiplier = 14791;//*/
	/*
	Uint64 modulus = (1ull << 29) - 3;
	Uint64 factors[] = { 2, 2, 7, 73, 262657, 0 };
	Uint64 chosen_multiplier = 78944;//*/
	/*
	Uint64 modulus = (1ull << 31) - 1;
	Uint64 factors[] = { 2, 3, 3, 7, 11, 31, 151, 331, 0 };
	Uint64 chosen_multiplier = 16807;//*/
	/*
	Uint64 modulus = (1ull << 33) - 9;
	Uint64 factors[] = { 2, 4294967291, 0 };
	Uint64 chosen_multiplier = 10349;//*/
	/*
	Uint64 modulus = (1ull << 36) - 5;
	Uint64 factors[] = { 2, 5, 6871947673, 0 };
	Uint64 chosen_multiplier = 2027339;//*/
	/*
	Uint64 modulus = (1ull << 39) - 7;
	Uint64 factors[] = { 2, 2, 2, 3, 3, 3, 5, 7, 13, 19, 37, 73, 109, 0 };
	Uint64 chosen_multiplier = 42017;//*/
	/*
	Uint64 modulus = 1099511627689ull;//2 to the 40th power minus 87
	Uint64 factors[] = { 2, 2, 2, 3, 3, 1487, 10269667, 0 };
	Uint64 chosen_multiplier = 234283;//*/
	/*
	Uint64 modulus = 17592186044399ull;//2 to the 44th power minus 17
	Uint64 factors[] = { 2, 7, 4583, 7993, 34303, 0 };
	Uint64 chosen_multiplier = 89164;//*/
	//*
	/*
	Uint64 modulus = 1125899906842079ull;//2 to the 50th power minus 545
	Uint64 factors[] = { 2, 2467, 228192117317ull, 0 };
	Uint64 chosen_multiplier = 11110;//*/
	//*

	long num_factors = 0;
	Uint64 product = 1;
	for (; factors[num_factors]; num_factors++) {
		Uint64 factor = factors[num_factors];
		product *= factor;
	}
	if (product != modulus - 1) issue_error("factors incorrect (1)");
	if (chosen_multiplier) {
		if (_expmod64(chosen_multiplier, product, modulus) != 1) issue_error("factors incorrect (2)");
		check_np2lcg_optimization(modulus, chosen_multiplier, factors, num_factors);
		return;//
	}
	Uint64 max_k = modulus;
	if (max_k > 99999) max_k = 99999;
	for (Uint64 k = 2; k < max_k; k++) {
		bool passing = true;
		if (modulus < (1ull << 32)) {
			if (_expmod32(k, product, modulus) != 1) continue;
			for (long fi = 0; fi < num_factors; fi++) {
				if (_expmod32(k, product / factors[fi], modulus) == 1) { passing = false; fi = num_factors; }
			}
		}
		else {
			if (_expmod64(k, product, modulus) != 1) continue;
			for (long fi = 0; fi < num_factors; fi++) {
				if (_expmod64(k, product / factors[fi], modulus) == 1) { passing = false; fi = num_factors; }
			}
		}
		if (!passing) continue;
		std::printf("%lld is full cycle for modulus %lld\n", k, modulus);
	}
}

int main(int argc, char **argv) {
	PractRand::initialize_PractRand();
	std::printf("initialized\n");
	PractRand::self_test_PractRand();
	std::printf("self-test completed\n");
	PractRand::RNGs::Polymorphic::sfc64 known_good(PractRand::SEED_AUTO);

	//nearseq2(); nearseq2(); nearseq2(); nearseq2(); nearseq2(); nearseq2(); nearseq2(); nearseq2();
	//test_sfc16();
	//blah_bcfn();
	//blah_fpf_all2();
	//blah_fpf();
	//find_test_distributions();
	//count_period(new PractRand::RNGs::Polymorphic::NotRecommended::xlcg8of64_varqual(28));
	//compare_bernoulli_distribution_calculations();
	//validate_lowest_of_N(&known_good);
	//for (int x = 0; x <= 22; x+=2) for (int y = 1; y <= 6 && y <= 1+(16/(x+4)); y++) compare_set_uniformity_to_table_uniformity(1 << x, 1 << y, &known_good);
	//verify_test_distributions();
	test_normal_distribution_a();
	//count_hash_deviation();
	//check_np2lcg_multiplier();
	//check_optimization();
	//print_data();

	/*	
	for (int i = 0; i < 117; i++) {
		std::printf(" %5f : %5f %5f %5f %5f %5f %5f\n", ref_p[i], Tests::math_pvalue_to_chisquared(ref_p[i], 1), 
			Tests::math_pvalue_to_chisquared(ref_p[i], 2), 
			Tests::math_pvalue_to_chisquared(ref_p[i], 3), 
			Tests::math_pvalue_to_chisquared(ref_p[i], 4),
			Tests::math_pvalue_to_chisquared(ref_p[i], 5),
			Tests::math_pvalue_to_chisquared(ref_p[i], 6)
		);
	}//*/


	return 0;
}




