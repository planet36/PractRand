#include <string>
//#include <ostream>
#include <sstream>
#include <vector>
#include <array>
#include <list>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/tests.h"
#include "PractRand/test_helpers.h"
#include "PractRand/test_batteries.h"

#include "PractRand/RNGs/arbee.h"

#include "PractRand/Tests/Gap16.h"
#include "PractRand/Tests/DistC6.h"
#include "PractRand/Tests/Pat5.h"
#include "PractRand/Tests/BCFN.h"
#include "PractRand/Tests/BCFN_MT.h"
#include "PractRand/Tests/FPF.h"
#include "PractRand/Tests/FPMulti.h"
#include "PractRand/Tests/Birthday.h"
#include "PractRand/Tests/CoupGap.h"
#include "PractRand/Tests/BRank.h"
#include "PractRand/Tests/mod3.h"
#include "PractRand/Tests/NearSeq.h"
#include "PractRand/Tests/coup16.h"
#include "PractRand/Tests/DistFreq4.h"
#include "PractRand/Tests/transforms.h"
#include "PractRand/Tests/permutation.h"

using namespace PractRand;
using namespace PractRand::Internals;

static long double gap_probs(int first, int last, long double chance_of_not_gap1) {
	return std::pow(chance_of_not_gap1, (long double)first) - std::pow(chance_of_not_gap1, (long double)last + 1);
}
static long double gap_probs64(Uint64 first, Uint64 last, long double chance_of_not_gap1) {
	return std::pow(chance_of_not_gap1, (long double)first) - std::pow(chance_of_not_gap1, (long double)last + 1);
}
static long double gap_expected(long double chance_of_gap1 = 1.0 / 65536) {
	long double e = 0;
	long double t = 0;
	long double end = 20.0 / chance_of_gap1 + 100;
	for (long double i = 1; i < end; i++) {
		long double p = chance_of_gap1 * std::pow(1 - chance_of_gap1, i - 1);
		e += i * p;
		t += p;
	}
	return e;
}
static long double gap_variance(long double chance_of_gap1 = 1.0 / 65536) {
	long double e = gap_expected(chance_of_gap1);
	long double var = 0;
	long double t = 0;
	long double end = 20.0 / chance_of_gap1 + 100;
	for (long double i = 1; i < end; i++) {
		long double p = chance_of_gap1 * std::pow(1 - chance_of_gap1, i - 1);
		long double sqr = i - e; sqr *= sqr;
		var += sqr * p;
		t += p;
	}
	return var;
}
static long double gap_log2_expected(long double chance_of_gap1 = 1.0 / 65536) {
	long double Le = 0;
	long double t = 0;
	long double end = 30.0 / chance_of_gap1 + 100;
	long double partial_Le = 0;
	long double partial_t = 0;
	for (long double i = 1; i < end; i++) {
		long double p = chance_of_gap1 * std::pow(1 - chance_of_gap1, i - 1);
		partial_Le += std::log2(i) * p;
		partial_t += p;
		if (partial_t >= 0.001) {
			Le += partial_Le;
			t += partial_t;
			partial_Le = 0;
			partial_t = 0;
		}
	}
	Le += partial_Le;
	t += partial_t;
	return Le;
}
static long double gap_log2_variance(long double chance_of_gap1 = 1.0 / 65536) {
	long double Le = gap_log2_expected(chance_of_gap1);
	long double var = 0;
	long double t = 0;
	long double end = 20.0 / chance_of_gap1 + 100;
	long double partial_var = 0;
	long double partial_t = 0;
	for (long double i = 1; i < end; i++) {
		long double p = chance_of_gap1 * std::pow(1 - chance_of_gap1, i - 1);
		long double sqr = std::log2(i) - Le; sqr *= sqr;
		partial_var += sqr * p;
		partial_t += p;
		if (partial_t >= 0.001) {
			var += partial_var;
			t += partial_t;
			partial_var = 0;
			partial_t = 0;
		}
	}
	var += partial_var;
	t += partial_t;
	return var;
}

static void truncate_table_bits(Uint64 *counts, double *probs, int old_bits, int new_bits) {
	int ns = 1 << new_bits;
	int os = 1 << old_bits;
	if (probs) for (int i = ns; i < os; i++) {
		int ni = i & (ns - 1);
		counts[ni] += counts[i];
		counts[i] = 0;
		probs[ni] += probs[i];
		probs[i] = 0;
	}
	else for (int i = ns; i < os; i++) {
		int ni = i & (ns - 1);
		counts[ni] += counts[i];
		counts[i] = 0;
	}
}
class TableIndexMasker {
	Uint32 lookup0[256];
	Uint32 lookup1[256];
	Uint32 lookup2[256];
	Uint32 lookup3[256];
public:
	void set_mask(int old_bits, Uint32 mask) {
		int bit_to_bit_lookup[32];
		int new_bit = 0;
		for (int old_bit = 0; old_bit < 32; old_bit++) {
			if (mask & (1 << old_bit)) bit_to_bit_lookup[old_bit] = new_bit++;
			else bit_to_bit_lookup[old_bit] = -1;
		}
		for (int i = 0; i < 256; i++) {
			Uint32 transformed = 0;
			for (int b = 0; b < 8; b++) if (i & (1 << b)) transformed |= Uint32(i) << bit_to_bit_lookup[b + 0];
			lookup0[i] = transformed;
		}
		for (int i = 0; i < 256; i++) {
			Uint32 transformed = 0;
			for (int b = 0; b < 8; b++) if (i & (1 << b)) transformed |= Uint32(i) << bit_to_bit_lookup[b + 8];
			lookup1[i] = transformed;
		}
		for (int i = 0; i < 256; i++) {
			Uint32 transformed = 0;
			for (int b = 0; b < 8; b++) if (i & (1 << b)) transformed |= Uint32(i) << bit_to_bit_lookup[b + 16];
			lookup2[i] = transformed;
		}
		for (int i = 0; i < 256; i++) {
			Uint32 transformed = 0;
			for (int b = 0; b < 8; b++) if (i & (1 << b)) transformed |= Uint32(i) << bit_to_bit_lookup[b + 24];
			lookup3[i] = transformed;
		}
	}
	Uint32 transform_index(Uint32 old_index) { return lookup0[Uint8(old_index >> 0)] | lookup1[Uint8(old_index >> 8)] | lookup2[Uint8(old_index >> 16)] | lookup3[Uint8(old_index >> 24)]; }
};




std::string PractRand::Tests::Gap16::get_name() const {
	return std::string("Gap-16");
}
void PractRand::Tests::Gap16::init( PractRand::RNGs::vRNG *known_good ) {
	int i;
	for (i = 0; i < 65536; i += 1) last[i] = 0;
	counts.reset_counts();
	autofail = false;
	blocks_tested = 0;
	warmup = 65536;
	extreme_lags.clear();
}
void PractRand::Tests::Gap16::increment_lag(Uint32 lag) {
	//if (lag < 1<<24) {index = (lag >> lookup_shift[lag>>16]) + lookup_offset[lag>>16];}
	//if (lag < (SIZE1 << SET1_SHIFT)) {//this should be true about 99.97% of the time
	//	counts.increment(lag >> SET1_SHIFT);
	//}
	//else {
		unsigned long tmp = lag - (SIZE1 << SET1_SHIFT);
		if (tmp < SIZE2 << SET2_SHIFT) {
			counts.increment(SIZE1 + (tmp >> SET2_SHIFT));
		}
		else {
			tmp = (tmp - (SIZE2 << SET2_SHIFT)) >> SET3_SHIFT;
			if (tmp < SIZE3) counts.increment(SIZE1 + SIZE2 + tmp);
			else extreme_lags.push_back(lag);
		}
	//}
}
void PractRand::Tests::Gap16::test_blocks(TestBlock *data, int numblocks) {
	if (autofail) return;
	while (warmup && numblocks > 64) {
		test_blocks(data, 64);
		data += 64;
		numblocks -= 64;
	}
	Uint32 ofs = Uint32(blocks_tested) * (TestBlock::SIZE / 2);
	unsigned long max = TestBlock::SIZE * numblocks / 2;
	Uint32 max2 = ofs + max;
	Uint16 *base = &data[0].as16[0];
	if (warmup) while (ofs != max2) {//warmup should last about 1.5 megabytes, but be highly variable
		Uint16 a = *(base++);
		Uint32 prior = last[a];
		last[a] = ++ofs;
		if (prior) {
			Uint32 lag = ofs - prior - 1;
			if (lag < (SIZE1 << SET1_SHIFT)) counts.increment(lag >> SET1_SHIFT);
			else increment_lag(lag);
		}
		else warmup--;
	}
	else while (ofs != max2) {
		Uint16 a;
		Uint32 prior, lag;

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		lag = ofs - prior - 1;
		if (lag < (SIZE1 << SET1_SHIFT)) counts.increment(lag >> SET1_SHIFT);
		else increment_lag(lag);

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		lag = Uint32(ofs - prior - 1);
		if (lag < (SIZE1 << SET1_SHIFT)) counts.increment(lag >> SET1_SHIFT);
		else increment_lag(lag);

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		lag = Uint32(ofs - prior - 1);
		if (lag < (SIZE1 << SET1_SHIFT)) counts.increment(lag >> SET1_SHIFT);
		else increment_lag(lag);

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		lag = Uint32(ofs - prior - 1);
		if (lag < (SIZE1 << SET1_SHIFT)) counts.increment(lag >> SET1_SHIFT);
		else increment_lag(lag);
	}
	Uint64 oblocks = blocks_tested;
	blocks_tested += numblocks;
	if ((oblocks>>19) != (blocks_tested>>19)) {//once every 512 megabytes or so... prevent overflow
		if (warmup) autofail = true;
		for (int i = 0; i < 65536; i++) {
			Uint32 n = last[i];
			if (Uint32(ofs - n) > 0x18000000) {
				autofail = true;
			}
		}
	}
}
//double PractRand::Tests::Gap16::get_result() {
void PractRand::Tests::Gap16::get_results( std::vector<TestResult> &results ) {
		//total weight: 1.000001
	if (blocks_tested < 3) return;
	double baseprob = 65535.0 / 65536.0;
	enum {TSIZE = SIZE1 + SIZE2 + SIZE3};
	std::vector<double> probs; probs.resize(TSIZE);
	//if (autofail) return 9876543210.;
	if (autofail) {
		results.push_back(TestResult(this->get_name() + ":!", autofail, autofail, TestResult::TYPE_PASSFAIL, 0.000001));
		return;
	}
	//correct probs for startup region:
	double lopped = 0;
	double inv_total_samples = 1.0 / (blocks_tested * (TestBlock::SIZE/2.));
	long highest_index = 0;
	for (int i = 0; i < TSIZE; i++) {
		int first, last;
		if (i < SIZE1) {
			first = i << SET1_SHIFT;
			last = first + (1 << SET1_SHIFT) - 1;
		}
		else if (i < SIZE1+SIZE2) {
			int tmp = i - SIZE1;
			first = (SIZE1 << SET1_SHIFT) + (tmp << SET2_SHIFT);
			last = first + (1 << SET2_SHIFT) - 1;
		}
		else if (i < TSIZE-1) {
			int tmp = i - SIZE1 - SIZE2;
			first = (SIZE1 << SET1_SHIFT) + (SIZE2 << SET2_SHIFT) + (tmp << SET3_SHIFT);
			last = first + (1 << SET3_SHIFT) - 1;
		}
		else {
			first = (SIZE1 << SET1_SHIFT) + (SIZE2 << SET2_SHIFT) + ((SIZE3-1) << SET3_SHIFT);
			last = 123456789;
		}
		double fraction = ((first+last)/2.+1.) * inv_total_samples;
		if (fraction > 1) fraction = 1;
		double p = gap_probs(first, last, baseprob);
		lopped += p * fraction;
		probs[i] = p * (1 - fraction);
		if (probs[i]) highest_index = i;
	}
	for (int i = 0; i < TSIZE; i++) {
		probs[i] /= 1 - lopped;
	}

	const Uint64 *count_ = counts.get_array();
	std::vector<Uint64> count; count.resize(TSIZE);
	Uint64 total_counts = 0;
	for (int i = 0; i < TSIZE; i++) {
		count[i] = count_[i];
		total_counts += count[i];
	}
	total_counts += extreme_lags.size();
	count[TSIZE - 1] += extreme_lags.size();
	if (!total_counts) {
		results.push_back(TestResult(get_name() + ":?", total_counts != 0, total_counts != 0, blocks_tested > 16 ? TestResult::TYPE_PASSFAIL : TestResult::TYPE_UNKNOWN, 0.0001));
		return;
	}

	double r_fine, r_coarse;
	double r_rarity = rarity_test(highest_index + 1, &probs[0], &counts[0], true, true);
	if (extreme_lags.size()) {
		double excess_gap = 0;
		for (int i = 0; i < extreme_lags.size(); i++) excess_gap += extreme_lags[i] - (1 << 17);
		//r_rarity += excess_gap * ??? / std::sqrt(double(total_counts));
	}
	//double r_uniformity = uniformity_test(highest_index + 1, &probs[0], &counts[0]);
	//double r_uniformity = uniformity_test_with_brute_force(highest_index + 1, &probs[0], &counts[0], &PractRand::RNGs::Polymorphic::arbee(13ull));
	int reduced_size = simplify_prob_table(
		highest_index + 1,
		(blocks_tested * TestBlock::SIZE/2.) / 40., 
		&probs[0], &count[0], true, true
	);
	r_fine = g_test(reduced_size, &probs[0], &count[0]);
//	double p1 = math_chisquared_to_pvalue(r1, reduced_size-1);
	r_fine = math_chisquared_to_normal(r_fine, reduced_size - 1);
	reduced_size = simplify_prob_table(
		reduced_size, 
		sqrt(blocks_tested * TestBlock::SIZE/2.), 
		&probs[0], &count[0], true, true
	);
	r_coarse = g_test(reduced_size, &probs[0], &count[0]);
//	double p2 = math_chisquared_to_pvalue(r2, reduced_size-1);
	r_coarse = math_chisquared_to_normal(r_coarse, reduced_size - 1);
	//double r;
	//if (r1 + 2.0 > fabs(r2)) r = r1 - 1;
	//else if (r1 - 2.0 < -fabs(r2)) r = r1 + 1;
	//else r = r2;
	//return r;
	TestCalibrationData *calib1 = calibration_manager.get_calibration_data("Gap-16:A", blocks_tested);
	TestCalibrationData *calib2 = calibration_manager.get_calibration_data("Gap-16:B", blocks_tested);
	double s1 = calib1->sample_to_suspicion(r_fine) * -1;
	double s2 = calib2->sample_to_suspicion(r_coarse) * -1;
	double cp1 = TestResult::suspicion_to_pvalue(s1);
	double cp2 = TestResult::suspicion_to_pvalue(s2);
	results.push_back(TestResult(get_name() + ":A", r_fine, s1, TestResult::TYPE_GOOD_S, 0.5));
	results.push_back(TestResult(get_name() + ":B", r_coarse, s2, TestResult::TYPE_GOOD_S, 0.5));
	results.push_back(TestResult(get_name() + ":C", r_rarity, r_rarity, TestResult::TYPE_BAD_S, 0.5));
	//results.push_back(TestResult(get_name() + ":D", r_uniformity, r_uniformity, TestResult::TYPE_BAD_S, 0.5));
}
/*double PractRand::Tests::Gap16::result_to_pvalue ( Uint64 blocks, double r ) {
	if (1) {//very crude aproximation:
		if (r < -9) return 1.0;
		if (r < -7) return 0.99999;
		if (r < -6) return 0.9999;
		if (r < -5) return 0.999;
		if (r < -4) return 0.99;
		if (r < -3) return 0.9;
		if (r > 12) return 0.0;
		if (r > 10) return 0.00001;
		if (r >  9) return 0.0001;
		if (r >  8) return 0.001;
		if (r >  7) return 0.01;
		if (r >  6) return 0.1;
		return 0.5;
	}
}*/




const double PractRand::Tests::GapUniversal1::expected_log_gap = 15.167378763679226;
const double PractRand::Tests::GapUniversal1::gap_log2_variance = 3.421308342472359;
std::string PractRand::Tests::GapUniversal1::get_name() const {
	return std::string("GapUni1");
}
void PractRand::Tests::GapUniversal1::init(PractRand::RNGs::vRNG *known_good) {
	TestBaseclass::init(known_good);
	int i;
	for (i = 0; i < 65536; i += 1) last[i] = 0;
	log_gap_sum = 0.0;
	num_gaps = 0;
	warmup = WARMUP;

	if (false) {
		long double expected = expected_log_gap;
		long double mean = 0;
		long double mean2 = 0;
		long double mean3 = 0;
		long double accum1 = 0;
		long double accum2 = 0;
		long double accum3 = 0;
		long double total_chance = 0;
		long biggest_gap = 65536 * 40;
		long double chance_modifier = 1 - std::pow((65535.0 / 65536), biggest_gap);
		for (int i = biggest_gap; i >= 1; i--) {
			long double lg = std::log2((long double)i);
			long double chance = std::pow(65535.0 / 65536, i) / 65535.0;
			chance /= chance_modifier;
			accum1 += lg * chance;
			if (accum1 > 0.001) { mean += accum1; accum1 = 0; }
			accum2 += std::pow(lg, 2.0) * chance;
			if (accum2 > 0.001) { mean2 += accum2; accum2 = 0; }
			accum3 += std::pow(lg - expected, 2.0) * chance;
			if (accum3 > 0.001) { mean3 += accum3; accum3 = 0; }
			// (A-B)^2 = A^2 - 2AB + B^2
			total_chance += chance;
		}
		mean += accum1; mean2 += accum2; mean3 += accum3;
		std::printf("\n mean = %.15f\n", double(mean));
		//std::printf("\n variance = %.15f\n", mean2 - mean * mean);
		//std::printf("\n variance = %.15f\n", mean3 + 2 * mean * expected - expected * expected - mean * mean);
		std::printf("\n variance = %.15f\n", double(mean3 - std::pow(mean - expected, 2.0)));
		//long double variance = mean3 - std::pow(mean - expected, 2.0);
		//long double cLK = 0.7 - 0.8 / 16 + (4 + 32.0 / 16) / 15.0 / std::pow(1024 * 1024 * 1024, 3.0 / 16);
		//std::printf("\n modified variance = %.15f\n", variance * cLK * cLK);
		std::printf("\n chance = %.15f\n", double(total_chance));
		std::printf("\n (1-chance_modifier) = %.15g\n", double(1 - chance_modifier));
		std::exit(0);
	}
}
void PractRand::Tests::GapUniversal1::handle_gap(Uint64 gap_value) {
	log_gap_sum += std::log2(double(gap_value)) - expected_log_gap;
	num_gaps += 1;
}
void PractRand::Tests::GapUniversal1::test_blocks(TestBlock *data, int numblocks) {
	Uint64 ofs = blocks_tested * (TestBlock::SIZE / 2);
	Uint64 max = TestBlock::SIZE * numblocks / 2;
	Uint64 max2 = ofs + max;
	Uint16 *base = &data[0].as16[0];
	while (warmup && ofs != max2) {
		Uint16 a = *(base++);
		last[a] = ++ofs;
		warmup -= 1;
	}
	while (ofs != max2) {
		Uint16 a = *(base++);
		Uint64 prior = last[a];
		last[a] = ++ofs;
		handle_gap(ofs - prior);
	}
	Uint64 oblocks = blocks_tested;
	blocks_tested += numblocks;
}
void PractRand::Tests::GapUniversal1::get_results(std::vector<TestResult> &results) {
	if (!num_gaps) return;
	long double cLK = 0.7 - 0.8 / 16 + (4 + 32.0 / 16) / 15.0 / std::pow(num_gaps, 3.0 / 16);
	double normalized_lg = log_gap_sum / (cLK * std::sqrt(double(num_gaps) * gap_log2_variance));
	results.push_back(TestResult(get_name(), normalized_lg, math_normaldist_to_pvalue(normalized_lg), TestResult::TYPE_GOOD_P, 1.0));
}



const double PractRand::Tests::GapUniversal2::expected_inverse = 1.0 / 36798.99546100125;
const double PractRand::Tests::GapUniversal2::gap_log2_variance = 3.421308342472359;
std::string PractRand::Tests::GapUniversal2::get_name() const {
	return std::string("GapUni2");
}
void PractRand::Tests::GapUniversal2::init(PractRand::RNGs::vRNG *known_good) {
	TestBaseclass::init(known_good);
	int i;
	for (i = 0; i < 65536; i += 1) last[i] = DUMMY_VALUE32;
	for (i = 0; i < 65536; i += 1) first[i] = 0;
	for (i = 0; i < 65536; i += 1) state[i] = 0;
	gap_product = 1.0;
	gap_product_extracted_L2 = 0;
	autofail = false;
	warmup = 65536;
}
void PractRand::Tests::GapUniversal2::normalize() {
	int n;
	gap_product = std::frexp(gap_product, &n);
	gap_product_extracted_L2 += n;
}
void PractRand::Tests::GapUniversal2::handle_gap(Uint32 gap_value) {
	gap_product *= gap_value * expected_inverse;
}
void PractRand::Tests::GapUniversal2::test_blocks(TestBlock *data, int numblocks) {
	if (autofail) return;
	Uint32 ofs = Uint32(blocks_tested) * (TestBlock::SIZE / 2);
	unsigned long max = TestBlock::SIZE * numblocks / 2;
	Uint32 max2 = ofs + max;
	Uint16 *base = &data[0].as16[0];
	if (warmup) while (ofs != max2) {
		Uint16 a;
		Uint32 prior, lag;

		a = *(base++);
		ofs += 1;
		prior = last[a];
		if (prior == DUMMY_VALUE32) {
			if (state[a] == 0) {
				first[a] = ofs;
				state[a] = 1;
			}
			else if (state[a] == 1) {
				last[a] = ofs;
				handle_gap( Uint32(ofs - first[a]));
				state[a] = 2;
				warmup -= 1;
			}
			else {
				issue_error("GapUniversal2: during warmup, prior == DUMMY_VALUE32 but state > 1");
			}
		}
		else {
			last[a] = ofs;
			handle_gap(Uint32(ofs - prior));
		}
		if (!(ofs & 63)) normalize();
	}
	else while (ofs != max2) {
		Uint16 a;
		Uint32 prior;

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		if (!(ofs & 63)) normalize();
	}
	Uint64 oblocks = blocks_tested;
	blocks_tested += numblocks;
	if ((oblocks >> 19) != (blocks_tested >> 19)) {//once every 512 megabytes or so... prevent overflow
		for (int i = 0; i < 65536; i++) {
			Uint32 n = last[i];
			if (Uint32(ofs - n) > 0x18000000) {
				autofail = true;
			}
			if (first[i] == DUMMY_VALUE32) autofail = true;
		}
	}
}
void PractRand::Tests::GapUniversal2::get_results(std::vector<TestResult> &results) {
	if (autofail) {
		results.push_back(TestResult(get_name() + "!", 1, 1, TestResult::TYPE_PASSFAIL, 1.0));
		return;
	}
	if (!blocks_tested) return;
	normalize();
	Uint64 base_length = blocks_tested * TestBlock::SIZE / sizeof(Uint16);
	Uint64 running_gap = 0;
	Uint64 first_extra_gap = 0;
	double extra_lg = 0;
	Uint64 extra_hits = 0;
	for (int i = 0; i < 65536; i++) {
		if (state[i] == 0) {
			running_gap += base_length;
		}
		else {
			if (!first_extra_gap) {
				first_extra_gap = running_gap + first[i];
				if (state[i] == 1) running_gap = Uint32(base_length - first[i]);
				else if (state[i] == 2) running_gap = Uint32(base_length - last[i]);
			}
			else if (state[i] == 1) {
				Uint64 g = running_gap + first[i];
				extra_lg += std::log2(g * expected_inverse);
				extra_hits++;
				running_gap = Uint32(base_length - first[i]);
			}
			else if (state[i] == 2) {
				Uint64 g = running_gap + first[i];
				extra_lg += std::log2(g * expected_inverse);
				extra_hits++;
				running_gap = Uint32(base_length - last[i]);
			}
		}
	}
	if (true) {
		Uint64 g = running_gap + first_extra_gap;
		extra_lg += std::log2(g * expected_inverse);
		extra_hits++;
	}

	if (true) {
		std::ostringstream o;
		o << get_name() << "";
		double combined_lg = gap_product_extracted_L2 + std::log2(gap_product) + extra_lg;
		//long double cLK = 0.7 - 0.8 / 16 + (4 + 32.0 / 16) / 15.0 / std::pow(base_length, 3.0 / 16);
		combined_lg -= 0.7;
		//double normalized_lg = combined_lg / (std::sqrt(double(base_length) * gap_log2_variance));
		double normalized_lg = combined_lg / (std::sqrt(double(base_length) * gap_log2_variance * 0.3916));
		TestCalibrationData *calib = calibration_manager.get_calibration_data("GapUni2", blocks_tested);
		double s = calib->sample_to_suspicion(normalized_lg);
		results.push_back(TestResult(o.str(), normalized_lg, s, TestResult::TYPE_GOOD_S, 1.0));
	}

	//double total_lg = std::log2(gap_product) + gap_product_extracted_L2 + extra_lg;
	//double normalized_lg = total_lg / std::sqrt(double(base_length - 1) * gap_log2_variance);

	//results.push_back(TestResult(this->get_name() + ":LB:F", ))

	//results.push_back(TestResult(this->get_name(), normalized_lg, math_normaldist_to_pvalue(normalized_lg), TestResult::TYPE_GOOD_P, 1.0));

	/*const int N = 1 << 16;
	const double x = 1.0 / N;
	double guess = std::log2(N) - 1/6;
	long double sum = 0;
	long double sum2 = 0;
	for (int i = 1; i <= 1<<24; i++) {
	double chance = x * std::pow(1 - x, i - 1);
	double L = std::log2(double(i)) - guess;
	sum += chance * L;
	sum2 += std::pow(L, 2.0) * chance;
	double mean_log2 = sum + guess;
	double mean_square = sum2 - guess * guess + 2 * guess * mean_log2;
	double variance1 = mean_square - mean_log2*mean_log2;
	double variance2 = sum2 - sum * sum;
	if (!(i & (i - 1))) std::printf("\n %6d:  %.12f %.12f        dev: %.7f", i, mean_log2, mean_log2 * std::log(2.0), variance2);
	}
	//*/
}






//const double PractRand::Tests::GapUniversal3::expected_inverse = 1.0 / 36798.99546100125;
//const double PractRand::Tests::GapUniversal3::gap_log2_variance = 3.421308342472359;
std::string PractRand::Tests::GapUniversal3::get_name() const {
	return std::string("GapUni3");
}
void PractRand::Tests::GapUniversal3::init(PractRand::RNGs::vRNG *known_good) {
	int i;
	for (i = 0; i < 65536; i += 1) last[i] = DUMMY_VALUE32;
	for (i = 0; i < 65536; i += 1) first[i] = 0;
	for (i = 0; i < 65536; i += 1) state[i] = 0;
	count_pow2 = 0;
	if (true) {//mandatory, atm - theoretically these are all constants that could be hardwired in, but they're not atm, have to be calculated here
		//plus, almost all of this depends upon the value of LOW_BITS
		overall_mean = 0;
		overall_variance = 0;
		const long double expected_mean = 15.16737;
		long double overall_total_chance = 0;
		for (i = 0; i < 1 << LOW_BITS; i++) {
			low_bit_indexed_data[i].gap_product = 1.0;
			low_bit_indexed_data[i].gap_product_extracted_L2 = 0;
			low_bit_indexed_data[i].hits = 0;
			low_bit_indexed_data[i].expected_gap_inverse = 1.0 / (37288.5 + i);

			const long double root = 65535.0 / 65536;
			long double mean = 0;
			long double mean2 = 0;
			long double accum = 0;
			long double accum2 = 0;
			long double chance_accum = 0;
			long double total_chance = 0;
			long gap_value = i + 1;
			while (gap_value < 65536 * 40) {
				long double gl = std::log2(double(gap_value));
				long double chance = std::pow(root, gap_value) / 65535.0;
				accum += gl * chance;
				gl -= expected_mean;
				accum2 += (gl * gl) * chance;
				chance_accum += chance;
				if (accum > 0.001 / (1 << LOW_BITS)) { mean += accum; accum = 0; }
				if (accum2 > 0.001 / (1 << LOW_BITS)) { mean2 += accum2; accum2 = 0; }
				if (chance_accum > 0.001 / (1 << LOW_BITS)) { total_chance += chance_accum; chance_accum = 0; }
				gap_value += (1 << LOW_BITS);
			}
			mean += accum; mean2 += accum2; total_chance += chance_accum;
			overall_mean += mean;
			overall_variance += mean2;
			overall_total_chance += total_chance;
			mean /= total_chance;
			mean2 /= total_chance;
			long double variance = mean2 - std::pow(mean - expected_mean, 2.0);
			low_bit_indexed_data[i].expected_gap_inverse = std::pow(0.5, mean);
			mean_lg[i] = mean;
			freq_probs[i] = total_chance;
			variances[i] = variance;
		}
		overall_variance -= std::pow(overall_mean - expected_mean, 2.0);
	}

	autofail = false;
	blocks_tested = 0;
	warmup = 65536;

}
void PractRand::Tests::GapUniversal3::normalize(LBI &lbi) {
	if (!std::isnormal(lbi.gap_product)) autofail = true;
	int n;
	lbi.gap_product = std::frexp(lbi.gap_product, &n);
	lbi.gap_product_extracted_L2 += n;
}
void PractRand::Tests::GapUniversal3::handle_gap(long gap_value) {
	LBI &lbi = low_bit_indexed_data[(gap_value - 1) & ((1 << LOW_BITS)-1)];
	lbi.hits++;
	lbi.gap_product *= gap_value * lbi.expected_gap_inverse;
	count_pow2 += (gap_value & (gap_value - 1)) ? 0 : 1;
	if (!(lbi.hits & 63)) normalize(lbi);
}
void PractRand::Tests::GapUniversal3::test_blocks(TestBlock *data, int numblocks) {
	if (autofail) return;
	Uint32 ofs = Uint32(blocks_tested) * (TestBlock::SIZE / 2);
	unsigned long max = TestBlock::SIZE * numblocks / 2;
	Uint32 max2 = ofs + max;
	Uint16 *base = &data[0].as16[0];
	if (warmup) while (ofs != max2) {
		Uint16 a;
		Uint32 prior, lag;

		a = *(base++);
		ofs += 1;
		prior = last[a];
		if (prior == DUMMY_VALUE32) {
			if (state[a] == 0) {
				first[a] = ofs;
				state[a] = 1;
			}
			else if (state[a] == 1) {
				last[a] = ofs;
				handle_gap(Uint32(ofs - first[a]));
				state[a] = 2;
				warmup -= 1;
			}
			else {
				issue_error("GapUniversal3: during warmup, prior == DUMMY_VALUE32 but state > 1");
			}
		}
		else {
			last[a] = ofs;
			handle_gap(Uint32(ofs - prior));
		}
	}
	else while (ofs != max2) {
		Uint16 a;
		Uint32 prior;

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));

		a = *(base++);
		prior = last[a];
		last[a] = ++ofs;
		handle_gap(Uint32(ofs - prior));
	}
	Uint64 oblocks = blocks_tested;
	blocks_tested += numblocks;
	if ((oblocks >> 19) != (blocks_tested >> 19)) {//once every 512 megabytes or so... prevent overflow
		for (int i = 0; i < 65536; i++) {
			Uint32 n = last[i];
			if (Uint32(ofs - n) > 0x18000000) {
				autofail = true;
			}
			if (first[i] == DUMMY_VALUE32) autofail = true;
		}
	}
}
void PractRand::Tests::GapUniversal3::get_results(std::vector<TestResult> &results) {
	if (autofail) {
		results.push_back(TestResult(get_name() + "!", 1, 1, TestResult::TYPE_PASSFAIL, 1.0));
		return;
	}
	if (!blocks_tested) return;
	for (int i = 0; i < (1 << LOW_BITS); i++) {
		normalize(low_bit_indexed_data[i]);
	}
	Uint64 base_length = blocks_tested * TestBlock::SIZE / sizeof(Uint16);
	Sint64 _base_length = 0; for (int i = 0; i < (1 << LOW_BITS); i++) _base_length += low_bit_indexed_data[i].hits; 
	LBI rollback_data[1 << LOW_BITS];
	for (int i = 0; i < (1 << LOW_BITS); i++) rollback_data[i] = low_bit_indexed_data[i];
	Uint64 running_gap = 0;
	Uint64 first_extra_gap = 0;
	//for (int i = 0; i < (1 << LOW_BITS); i++) {extra_lg[i] = 0; extra_hits[i] = 0;}
	for (int i = 0; i < 65536; i++) {
		if (state[i] == 0) {
			running_gap += base_length;
		}
		else {
			if (!first_extra_gap) {
				first_extra_gap = running_gap + first[i];
				if (state[i] == 1) running_gap = Uint32(base_length - first[i]);
				else if (state[i] == 2) running_gap = Uint32(base_length - last[i]);
			}
			else if (state[i] == 1) {
				Uint64 g = running_gap + first[i];
				handle_gap(g);
				running_gap = Uint32(base_length - first[i]);
			}
			else if (state[i] == 2) {
				Uint64 g = running_gap + first[i];
				handle_gap(g);
				running_gap = Uint32(base_length - last[i]);
			}
		}
	}
	if (true) {
		Uint64 g = running_gap + first_extra_gap;
		handle_gap(g);
	}

	int target_lowbits = (std::log2(base_length) - 10) / 3;
	if (target_lowbits < 0) target_lowbits = 0;
	if (target_lowbits > LOW_BITS) target_lowbits = LOW_BITS;
	std::vector<double> reduced_lg; reduced_lg.resize(1 << target_lowbits, 0);
	std::vector<double> reduced_variance; reduced_variance.resize(1 << target_lowbits, 0);
	std::vector<Uint64> reduced_hits; reduced_hits.resize(1 << target_lowbits, 0);
	//std::vector<double> reduced_lg_bias; reduced_lg_bias.resize(1 << target_lowbits, 0);

	for (int i = 0; i < (1 << LOW_BITS); i++) {
		int i2 = i & ((1 << target_lowbits) - 1);
		reduced_hits[i2] += low_bit_indexed_data[i].hits;
		reduced_lg[i2] += low_bit_indexed_data[i].gap_product_extracted_L2 + std::log2(low_bit_indexed_data[i].gap_product);
		reduced_variance[i2] += variances[i] * low_bit_indexed_data[i].hits;
	}

	double single_lg = 0;
	double single_variance = 0;
	Uint64 single_hits = 0;
	double single_lg_bias = 0;
	for (int i = 0; i < (1 << LOW_BITS); i++) {
		single_lg += low_bit_indexed_data[i].gap_product_extracted_L2 + std::log2(low_bit_indexed_data[i].gap_product);
		single_variance += low_bit_indexed_data[i].hits * variances[i];
		single_hits += low_bit_indexed_data[i].hits;
		single_lg_bias += low_bit_indexed_data[i].hits * (mean_lg[i] - overall_mean);
	}

	//=================================================================================
	for (int i = 0; i < (1 << LOW_BITS); i++) low_bit_indexed_data[i] = rollback_data[i];
	//=================================================================================

	double lowest_lbr, highest_lbr;
	double lowest_lbp = 2.0;
	double highest_lbp = -1.0;
	//std::vector<double> reduced_normalized_lg; reduced_normalized_lg.resize(1 << target_lowbits, 0);
	for (int i = 0; i < (1 << target_lowbits); i++) {
		//single_lg += reduced_lg[i];
		//single_hits += reduced_hits[i];
		//reduced_normalized_lg[i] 
		if (!reduced_variance[i]) autofail = true;
		double r = reduced_lg[i] / std::sqrt(reduced_variance[i]);
		double p = math_normaldist_to_pvalue(r);
		if (p < lowest_lbp) {
			lowest_lbp = p;
			lowest_lbr = r;
		}
		if (p > highest_lbp) {
			highest_lbp = p;
			highest_lbr = r;
		}
	}
	if (autofail) {
		results.push_back(TestResult(get_name() + "!", 1, 1, TestResult::TYPE_PASSFAIL, 1.0));
		return;
	}
	if (target_lowbits > 1) {
		std::ostringstream o;
		o << get_name() << ":" << target_lowbits << ":high";
		results.push_back(TestResult(o.str(), highest_lbr, 1 - std::pow(1 - highest_lbp, 1.0 / (1 << target_lowbits)), TestResult::TYPE_BAD_P, 1.0));
	}
	if (target_lowbits > 1) {
		std::ostringstream o;
		o << get_name() << ":" << target_lowbits << ":low";
		results.push_back(TestResult(o.str(), lowest_lbr, std::pow(lowest_lbp, 1.0 / (1 << target_lowbits)), TestResult::TYPE_BAD_P, 1.0));
	}
	if (true) {
		std::ostringstream o;
		double single_normalized_lg = single_lg / std::sqrt(single_variance);
		results.push_back(TestResult(this->get_name() + ":0:sin", single_normalized_lg, math_normaldist_to_pvalue(single_normalized_lg), TestResult::TYPE_BAD_P, 1.0));
	}
}







std::string PractRand::Tests::LimitedBigGapPrototype::get_name() const {
	return std::string("LBGP");
}
void PractRand::Tests::LimitedBigGapPrototype::init(PractRand::RNGs::vRNG *known_good) {
	if (!last) last = new PractRand::Tests::LimitedBigGapPrototype::Entry[1ull << INDEX_BITS];
	if (!last) issue_error("LimitedBigGapPrototype - failed to allocated memory block 1");
	if (!first) first = new PractRand::Tests::LimitedBigGapPrototype::Entry[1ull << INDEX_BITS];
	if (!first) issue_error("LimitedBigGapPrototype - failed to allocated memory block 2");
	for (long i = 0; i < SIZE; i++) {
		last[i].position = Uint64(Sint64(-1));
		last[i].value = Uint64(Sint64(-1));
		first[i].position = Uint64(Sint64(-1));
		first[i].value = Uint64(Sint64(-1));
	}
	for (long i = 0; i < NUM_DISTANCE_BUCKETS * MATCH_LEVELS; i++) counts[i] = 0;
	warmup1 = warmup2 = 1 << INDEX_BITS;
	TestBaseclass::init(known_good);
}
void PractRand::Tests::LimitedBigGapPrototype::deinit() {
	delete[] last;
	delete[] first;
}
long PractRand::Tests::LimitedBigGapPrototype::distance_to_bucket(Uint64 gap_value) {
	long rv = ilog2_64(gap_value);
	return rv;
}
void PractRand::Tests::LimitedBigGapPrototype::bucket_to_distance_range(long distance_bucket, Uint64 &min_distance, Uint64 &max_distance) {
	//inclusive on both sides

	min_distance = 1ull << distance_bucket;
	max_distance = (min_distance << 1) - 1;
}
void PractRand::Tests::LimitedBigGapPrototype::test_blocks(TestBlock *data, int numblocks) {
	while (numblocks >= 65536) {
		test_blocks(data, 65536);
		data += 65536;
		numblocks -= 65536;
	}
	if (!numblocks) return;
	Uint64 base_position = blocks_tested * (TestBlock::SIZE / sizeof(Uint64));
	long max = numblocks * (TestBlock::SIZE / sizeof(Uint64));
	for (long pos = 0; pos < max; pos++) {
		Uint64 cur = data->as64[pos];
		if (cur & ((1 << ZERO_BITS) - 1)) continue;
		//if (TOTAL_BITS < 64) cur &= (1ull << TOTAL_BITS) - 1;
		long index = long(cur >> ZERO_BITS) & (SIZE - 1);
		Entry &e = last[index];
		Uint64 full_position = base_position + pos;
		Uint64 trimmed_value = cur >> (ZERO_BITS + INDEX_BITS);
		Uint64 prev_pos, prev_value;
		if (e.position == Uint64(Sint64(-1))) {
			Entry &e2 = first[index];
			if (e2.position == Uint64(Sint64(-1))) {
				// first occurance
				e2.position = full_position;
				e2.value = trimmed_value;
				warmup1--;
				continue;
			}
			else {
				// second occurance
				prev_pos = e2.position;
				prev_value = e2.value;
				warmup2--;
			}
		}
		else {
			// third or later occurance
			prev_pos = e.position;
			prev_value = e.value;
		}
		e.position = full_position;
		e.value = trimmed_value;

		Uint64 gap_distance = full_position - prev_pos;
		long distance_bucket = distance_to_bucket(gap_distance);
		Uint64 mismatch = prev_value ^ trimmed_value;
		long match_bits = count_low_zeroes64(mismatch);
		long match_level = match_bits >> CHECK_COUNT_SHIFT;
		if (match_level > MAX_MATCH_LEVEL) match_level = MAX_MATCH_LEVEL;
		long full_bucket = distance_bucket << MATCH_LEVELS_L2;
		full_bucket += match_level;
		counts[full_bucket]++;
	}
	blocks_tested += numblocks;
}
void PractRand::Tests::LimitedBigGapPrototype::get_results(std::vector<TestResult> &results) {
	Uint64 position = blocks_tested * (TestBlock::SIZE / sizeof(Uint64));
	enum {
		GAP_BITS = ZERO_BITS + INDEX_BITS, 
	};
	// any check across distances will be distorted due to the curvature of gap distances
	// things to check:
	//		0. basic sanity check - did we get about the right amount of samples?
	//		1. overall distribution of match levels
	//		2. gap distribution at individual match levels (including hits from higher levels?)
	//		3. test everything at once

	/*
		possible changes now:
			first occurance is now recorded for each index, including both position and value
			thus, with wrap-around, it is possible to get one more sample per index, though that must be done at get_results() time
			also, warmup now knows, for each index, whether it has 0 occurances, 1 occurance, or 2+ occurances
			this *could* be used in the same way GapUniversal2 does it, but I'd rather not
			instead I'd like to do a chi-squared test on the distribution of 0/1/2+ occurances, and merge the result with the g[0+] test.  
	*/

	Uint64 total_count = 0;
	for (long i = 0; i < NUM_DISTANCE_BUCKETS * MATCH_LEVELS; i++) total_count += counts[i];

	double expected_count = double(position) / (1ull << ZERO_BITS);
	if (1) {// 0. basic sanity check - did we get about the right amount of samples?
		double ratio = expected_count / (1ull << INDEX_BITS);
		expected_count -= (1ull << INDEX_BITS) * (1 - std::exp(-ratio));
		double difference = expected_count - total_count;
		double variance = expected_count;
		double norm = difference / std::sqrt(variance);

		if (expected_count < 100) return;
		results.push_back(TestResult("LBGP:sanity", norm, norm, TestResult::TYPE_RAW_NORMAL, 1.0));
		//std::printf("LBGP: tcount: %lld   warmup: %lld   log2(pos): %f\n", total_count, warmup, std::log2(double(position)));
	}

	Uint64 match_level_counts[MATCH_LEVELS];
	for (long i = 0; i < MATCH_LEVELS; i++) match_level_counts[i] = 0;
	for (long i = 0; i < NUM_DISTANCE_BUCKETS * MATCH_LEVELS; i++) match_level_counts[i & (MATCH_LEVELS - 1)] += counts[i];
	double match_level_prob_ratio = 1ull << (1 << CHECK_COUNT_SHIFT);
	double match_level_prob_ratio_inv = 1.0 / match_level_prob_ratio;
	if (1) {// 1. overall distribution of match levels
		double prob_table[MATCH_LEVELS];
		double remaining = 1.0;
		for (long i = 0; i < MATCH_LEVELS; i++) {
			prob_table[i] = remaining * (1 - match_level_prob_ratio_inv);
			remaining *= match_level_prob_ratio_inv;
		}
		prob_table[MATCH_LEVELS-1] += remaining;
		long cat = MATCH_LEVELS;
		Uint64 tmp_counts[MATCH_LEVELS];
		for (long i = 0; i < MATCH_LEVELS; i++) tmp_counts[i] = match_level_counts[i];
		cat = simplify_prob_table(cat, expected_count / 1000.0, prob_table, tmp_counts, true, true);
		double chis = g_test(cat, prob_table, tmp_counts);
		double norm = math_chisquared_to_normal(chis, cat - 1);
		results.push_back(TestResult("LBGP:cross", norm, norm, TestResult::TYPE_RAW_NORMAL, 1.0));
	}

	if (1) {// 2. gap distribution at individual match levels (including hits from higher levels?)
		double prob_table[NUM_DISTANCE_BUCKETS];
		double base_gap_prob = ((1ull << GAP_BITS) - 1) / double(1ull << GAP_BITS);//the probability that a random gap distance is greater than 1
		int max_cat = NUM_DISTANCE_BUCKETS;
		for (long i = 0; i < NUM_DISTANCE_BUCKETS; i++) {
			Uint64 min_dist, max_dist;
			bucket_to_distance_range(i, min_dist, max_dist);
			prob_table[i] = gap_probs64(min_dist - 1, max_dist - 1, base_gap_prob);//distances are counted weirdly atm, maybe because Gap16 stored them efficiently or something
			double mid_dist_ratio = position / ((double(min_dist) + double(max_dist)) * 0.5);
			double scale;
			scale = 1 - std::exp(1 - mid_dist_ratio);
			if (scale < 0) scale = 0;
			if (max_dist > position) scale = 0;
			prob_table[i] *= scale;
			if (!prob_table[i] && i < max_cat) max_cat = i;
		}
		double total_prob = 0;
		for (long i = 0; i < max_cat; i++) total_prob += prob_table[i];
		for (long i = 0; i < max_cat; i++) prob_table[i] /= total_prob;

		enum { HIGHER_LEVELS_STILL_COUNT_TOWARDS_LOWER = 1 };
		for (long match_level = 0; match_level < MATCH_LEVELS; match_level++) {
			double expected_at_match_level = expected_count * std::pow(match_level_prob_ratio_inv, match_level);
			if (HIGHER_LEVELS_STILL_COUNT_TOWARDS_LOWER && match_level != MATCH_LEVELS - 1) expected_at_match_level *= (match_level_prob_ratio - 1) / match_level_prob_ratio;
			if (expected_at_match_level < 1000) continue;
			double tmp_probs[NUM_DISTANCE_BUCKETS];
			for (long i = 0; i < NUM_DISTANCE_BUCKETS; i++) tmp_probs[i] = prob_table[i];
			Uint64 tmp_counts[NUM_DISTANCE_BUCKETS];
			for (long i = 0; i < NUM_DISTANCE_BUCKETS; i++) {
				tmp_counts[i] = counts[(i << MATCH_LEVELS_L2) + match_level];
				// also adding in all higher match levels:
				if (HIGHER_LEVELS_STILL_COUNT_TOWARDS_LOWER)
					for (long match_level2 = match_level + 1; match_level2 < MATCH_LEVELS; match_level2++) tmp_counts[i] += counts[(i << MATCH_LEVELS_L2) + match_level2];
			}
			double norm_u = uniformity_test(max_cat, tmp_probs, tmp_counts);
			double norm_r = rarity_test(max_cat, tmp_probs, tmp_counts, true, true);
			long cat = max_cat;
			cat = simplify_prob_table(cat, std::pow(expected_at_match_level, 0.5) / 40.0, tmp_probs, tmp_counts, true, true);
			double chis = g_test(cat, tmp_probs, tmp_counts);
			double norm_g = math_chisquared_to_normal(chis, cat - 1);
			results.push_back(TestResult(std::string("LBGP:gaps[") + std::to_string(match_level) + std::string("+]:g"), norm_g, norm_g, TestResult::TYPE_RAW_NORMAL, 1.0));
			//results.push_back(TestResult(std::string("LBGP:gaps[") + std::to_string(match_level) + std::string("+]:r"), norm_r, norm_r, TestResult::TYPE_RAW_NORMAL, 1.0));
			//results.push_back(TestResult(std::string("LBGP:gaps[") + std::to_string(match_level) + std::string("+]:u"), norm_u, norm_u, TestResult::TYPE_RAW_NORMAL, 1.0));
		}
	}

#if 0
	if (position < 1ull << 24) return;



	//results.push_back(TestResult("LG64P:thits", total_hits, 0, TestResult::TYPE_UNKNOWN, 1.0));
	//results.push_back(TestResult("LG64P:tmisses", total_misses, 0, TestResult::TYPE_UNKNOWN, 1.0));
	std::printf("LCG64P:thits: %lld   misses: %lld   warmup: %lld\n", total_hits, total_misses, warmup);

	//if (position < 1ull << 28) return;//not trying to deal with short tests yet
	if (total_hits < 10) return;

	double hit_prob = 1.0 / double(1ull << (CHECK_BITS));
	double base_gap_prob = ((1ull << GAP_BITS) - 1) / double(1ull << GAP_BITS);
	double prob_table[64];
	double total_prob;
	Uint64 tmp_counts[64];
	for (long i = 0; i < 64; i++) {
		//sanity checking distance buckets
		Uint64 firstdist = 1ull << i;
		Uint64 lastdist = (firstdist << 1) - 1;
		if (distance_to_bucket(firstdist) != i) issue_error();
		if (distance_to_bucket(lastdist) != i) issue_error();
		if (firstdist > 1 && distance_to_bucket(firstdist - 1) == i) issue_error();
		if (lastdist + 1 != 0 && distance_to_bucket(lastdist + 1) == i) issue_error();
	}
	for (long i = 0; i < 64; i++) {
		Uint64 firstdist = 1ull << i;
		Uint64 lastdist = (firstdist << 1) - 1;
		prob_table[i] = gap_probs64(firstdist, lastdist, base_gap_prob);
		tmp_counts[i] = hit_counts[i];
		/*if (lastdist >= position) {
			//tmp_counts[i - 1] += tmp_counts[i];
			tmp_counts[i] = 0;
			prob_table[i] = 0;
		}//*/
	}
	total_prob = 0;
	for (long i = 0; i < 64; i++) total_prob += prob_table[i];
	for (long i = 0; i < 64; i++) prob_table[i] /= total_prob;
	int cat;
	double chis, norm;
	cat = 64;
	cat = simplify_prob_table(64, (position >> (GAP_BITS + CHECK_BITS)) / 1000, prob_table, tmp_counts, true, true);
	chis = g_test(cat, prob_table, tmp_counts);
	norm = math_chisquared_to_normal(chis, cat - 1);
	results.push_back(TestResult("LG64P:hitdist", norm, norm, TestResult::TYPE_RAW_NORMAL, 1.0));

	for (long i = 0; i < 64; i++) {
		Uint64 firstdist = 1ull << i;
		Uint64 lastdist = (firstdist << 1) - 1;
		prob_table[i] = gap_probs64(firstdist, lastdist, base_gap_prob);
		tmp_counts[i] = miss_counts[i];
	}
	cat = 64;
	cat = simplify_prob_table(64, (position >> (ZERO_BITS + CHECK_BITS)) / 1000, prob_table, tmp_counts, true, true);
	chis = g_test(cat, prob_table, tmp_counts);
	norm = math_chisquared_to_normal(chis, cat - 1);
	//results.push_back(TestResult("LG64P:misdist", norm, norm, TestResult::TYPE_RAW_NORMAL, 1.0));

	if (CHECK_BITS < 1) return;
	double prob_table2[64];
	for (long i = 0; i < 64; i++) {
		Uint64 firstdist = 1ull << i;
		Uint64 lastdist = (firstdist << 1) - 1;
		prob_table2[i] = (miss_counts[i] + hit_counts[i]) / double(total_hits + total_misses);
		tmp_counts[i] = hit_counts[i];
	}
	cat = 64;
	cat = simplify_prob_table(64, (position >> (ZERO_BITS + CHECK_BITS)) / 1000, prob_table2, tmp_counts, true, true);
	chis = g_test(cat, prob_table2, tmp_counts);
	norm = math_chisquared_to_normal(chis, cat - 1);
	results.push_back(TestResult("LG64P:hitdist2", norm, norm, TestResult::TYPE_RAW_NORMAL, 1.0));
#endif
}








PractRand::Tests::DistC7::DistC7(int length_, int unitsL_, int bits_clipped_0_, int bits_clipped_1_, int bits_clipped_2_)
:
DistC6(length_, unitsL_, bits_clipped_0_, bits_clipped_1_, bits_clipped_2_)
{
}
void PractRand::Tests::DistC7::init(PractRand::RNGs::vRNG *known_good) {
	DistC6::init(known_good);
	odd_counts.set_size(counts.get_size());
	odd_counts.reset_counts();
	odd = false;
}
std::string PractRand::Tests::DistC7::get_name() const {
	std::ostringstream tmp;
	tmp << "DC7-" << length << "x" << (1 << unitsL) << "Bytes-" <<
		(bits_clipped_0 + 10 * bits_clipped_1 + 100 * bits_clipped_2);
	return tmp.str();
}
void PractRand::Tests::DistC7::test_blocks(TestBlock *data, int numblocks) {
	int max = numblocks * (TestBlock::SIZE >> unitsL);
	int i = 0;
	while (warmup) {
		//int max2 = 1 * (TestBlock::SIZE >> unitsL);
		int max2 = ((warmup + (1 << 4) - 1) >> 4) << 4;//round up to a multiple of 16
		if (max2 > max) max2 = max;
		if (!ENABLE_8_BIT_BYPASS || unitsL) while (max2 > i) {
			int bits;
			switch (unitsL) {
				case 0: bits = count_ones8(data->as8[i]); break;
				case 1: bits = count_ones16(data->as16[i]); break;
				case 2: bits = count_ones32(data->as32[i]); break;
				case 3: bits = count_ones64(data->as64[i]); break;
				default: {
					issue_error();
					bits = 0;//just to make the compiler happy
				} break;
			}
			advance_index(bits);
			i++;
			if (warmup) warmup--;
			else {
				if (odd) odd_counts.increment(last_index);
				else counts.increment(last_index);
			}
			odd = !odd;
		}
		else while (max2 > i) {
			last_index = _advance_index(last_index, lookup_table[data->as8[i++]]);
			if (warmup) warmup--;
			else {
				if (odd) odd_counts.increment(last_index);
				else counts.increment(last_index);
			}
			odd = !odd;
		}
	}
	if (odd) issue_error("DC7 - odd should be false post-warmup, right?");
	Uint32 index = last_index;
	switch (unitsL) {
		case 0: {//8bit
				if (ENABLE_8_BIT_BYPASS) for (; i < max;) {
					index = _advance_index(index, lookup_table[data->as8[i++]]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[data->as8[i++]]);
					odd_counts.increment(index);
					index = _advance_index(index, lookup_table[data->as8[i++]]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[data->as8[i++]]);
					odd_counts.increment(index);
				}
				else for (; i < max;) {
					index = _advance_index(index, lookup_table[count_ones8(data->as8[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones8(data->as8[i++])]);
					odd_counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones8(data->as8[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones8(data->as8[i++])]);
					odd_counts.increment(index);
				}
		}
		break;
		case 1: {//16bit
				for (; i < max;) {
					index = _advance_index(index, lookup_table[count_ones16(data->as16[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones16(data->as16[i++])]);
					odd_counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones16(data->as16[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones16(data->as16[i++])]);
					odd_counts.increment(index);
				}
		}
		break;
		case 2: {//32bit
				for (; i < max;) {
					index = _advance_index(index, lookup_table[count_ones32(data->as32[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones32(data->as32[i++])]);
					odd_counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones32(data->as32[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones32(data->as32[i++])]);
					odd_counts.increment(index);
				}
		}
		break;
		case 3: {//64bit
				for (; i < max;) {
					index = _advance_index(index, lookup_table[count_ones64(data->as64[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones64(data->as64[i++])]);
					odd_counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones64(data->as64[i++])]);
					counts.increment(index);
					index = _advance_index(index, lookup_table[count_ones64(data->as64[i++])]);
					odd_counts.increment(index);
				}
		}
		break;
		default: {
				 issue_error();
				 return;
		}
		break;
	}
	if (i != max) issue_error("DC7 went past end?");
	last_index = index;
	blocks_tested += numblocks;
}
void PractRand::Tests::DistC7::get_results(std::vector<TestResult> &results) {
	long initial_results_size = results.size();
	long old_results_size;
	Uint64 tmp = blocks_tested;
	blocks_tested >>= 1;
	if (blocks_tested) {
		old_results_size = results.size();
		DistC6::get_results(results);
		if (results.size() == old_results_size + 1) results.back().name += ":even";
		old_results_size = results.size();
		counts.swap_array(odd_counts); DistC6::get_results(results); counts.swap_array(odd_counts);
		if (results.size() == old_results_size + 1) results.back().name += ":odd";
	}
	blocks_tested = tmp;
	const Uint64 *evens = counts.get_array();
	const Uint64 *odds = odd_counts.get_array();
	VariableSizeCount<Uint8> tmp_counts; tmp_counts.set_size(size);
	for (int i = 0; i < size; i++) tmp_counts.force_count(i, evens[i] + odds[i]);
	old_results_size = results.size();
	counts.swap_array(tmp_counts); DistC6::get_results(results); counts.swap_array(tmp_counts);
	if (results.size() == old_results_size + 1) results.back().name += ":both";
	if (results.size() != initial_results_size + 3) return;
	double raw = -2 * (std::log(results[initial_results_size].get_pvalue()) + std::log(results[initial_results_size + 1].get_pvalue()));
	double n = math_chisquared_to_normal(raw, 4);
	results.push_back(TestResult(get_name() + ":indep", n, 1 - math_chisquared_to_pvalue(raw, 4), TestResult::TYPE_GOOD_P, 0.5));
	return;
}




#if 0
PractRand::Tests::BCFN_MT::BCFN_MT( int unitsL2_, int tbits_ ) {
	unitsL2 = unitsL2_;
	tbits = tbits_;
}
static std::vector<int> BCFN_MT_calculate_thresholds(int max_thresholds, Uint64 shift, int word_bits_L2, double target_fraction = 1.0/3) {
	Uint64 word_bits = 1ull << word_bits_L2;
	std::vector<int> rv;
	rv.push_back(1);
	if (word_bits <= 16384) {
		std::vector<double> pdf, cdf;
		Tests::get_hamming_weight_chances(1 << word_bits_L2, pdf, cdf);
		Uint64 n = word_bits/2;
		double target = cdf[n-1] * target_fraction;
		int max = n >> shift;
		for (Uint64 i = 1; i <= max; i++) {
			double cur = cdf[n-(i << shift)];
			if (cur >= target) continue;
			rv.push_back(i<<shift);
			if (rv.size() >= max_thresholds) break;
			target = cur * target_fraction;
		}
	}
	else {
		Uint64 n = word_bits/2;
		double target = (0.5 - 0.5 * Tests::calculate_center_bit_combination_chance(word_bits_L2)) * target_fraction;
		int max = n >> shift;
		double mean = word_bits/2.0;
		double dev = sqrt(word_bits * 0.5 * 0.5);
		double delta = 1.0 / dev;
		for (Uint64 i = 1; i < max; i++) {
			int ti = n+1-i;
			double norm = (ti - mean) / dev;
			double cur = Tests::math_normaldist_to_pvalue(norm + 0.5 * delta);
			if (cur >= target) continue;
			rv.push_back(i<<shift);
			if (rv.size() >= max_thresholds) break;
			target = cur * target_fraction;
		}
	}
	return rv;
}
void PractRand::Tests::BCFN_MT::init( PractRand::RNGs::vRNG *known_good ) {
	int tsize = 1 << tbits;
	mask = tsize - 1;
	for (int level = 0; level < LEVELS; level++) {
		int elevel = unitsL2 + level;
		if (true) {
			int shift = (elevel <= 4) ? 0 : (elevel/2 - 2);
			bitcount_shift[level] = shift;
			std::vector<int> thresholds = BCFN_MT_calculate_thresholds(THRESHOLDS, shift, elevel+3, 0.333);

			for (int threshold_index = 0; threshold_index < INDEX_SIZE; threshold_index++) threshold_lookup[threshold_index + level * THRESHOLDS] = 255;
			int threshold_n = 0;
			int next_threshold = thresholds[threshold_n]>>shift;//BCFNMT_find_next_threshold(1, (1 << shift), (1 << bitcount_shift[level]), 8 << elevel) >> shift;
			for (int threshold_index = 0; threshold_index < INDEX_SIZE; threshold_index++) {
				if (threshold_index == next_threshold) {
					if (threshold_n < THRESHOLDS-1) threshold_n++;
					next_threshold = threshold_n < thresholds.size() ? thresholds[threshold_n]>>shift : 0;
					//next_threshold = BCFNMT_find_next_threshold(next_threshold << shift, next_threshold << shift, 1 << shift, 8 << elevel) >> shift;
				}
				threshold_lookup[threshold_index + level * THRESHOLDS] = threshold_n;
			}
		}
		for (int threshold = 0; threshold < THRESHOLDS; threshold++) {
			int lti = level * THRESHOLDS + threshold;
			warmup[lti] = tbits - 1;
			cur[lti] = 0;
			//total[lti] = 0;
			counts[lti].set_size(tsize);
			counts[lti].reset_counts();
		}
		odd[level] = false;
		leftovers[level] = 0;
	}
	blocks_tested = 0;
}
void PractRand::Tests::BCFN_MT::deinit() {
	for (int i = 0; i < LEVELS; i++) {
		counts[i].reset_counts();
	}
}
std::string PractRand::Tests::BCFN_MT::get_name() const {
	std::ostringstream f;
	f << "BCFN_MT(" << unitsL2 << "," << "," << tbits << ")";
	return f.str();
	//return make_string("BCFN-%d", tbits);
}
void PractRand::Tests::BCFN_MT::test_blocks(TestBlock *data, int numblocks) {
	//while (warmup[4]) {
	while (true) {
#define GET_BITS8(pos)  (count_ones8 (data[0].as8 [i+(pos)]) - 4)
#define GET_BITS16(pos) (count_ones16(data[0].as16[i+(pos)]) - 8)
#define GET_BITS32(pos) (count_ones32(data[0].as32[i+(pos)]) - 16)
#define GET_BITS64(pos) (count_ones64(data[0].as64[i+(pos)]) - 32)
//#define HANDLE_BITS(level,var) if (var){tmp=var>>31;cur[level]=((cur[level]<<1)-tmp)&mask[level];if (warmup[level]) warmup[level]--; else counts[level].increment(cur[level]);}
		switch (unitsL2) {
			case 0: {
				for (int i = 0; i < TestBlock::SIZE / 1; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS8(0);
					//HANDLE_BITS(0,bits0);
					handle_high_levels(0, bits0);
				}
			}
			break;
			case 1: {
				for (int i = 0; i < TestBlock::SIZE / 2; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS16(0);
					//HANDLE_BITS(0,bits0);
					handle_high_levels(0, bits0);
				}
			}
			break;
			case 2: {
				for (int i = 0; i < TestBlock::SIZE / 4; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS32(0);
					//HANDLE_BITS(0,bits0);
					handle_high_levels(0, bits0);
				}
			}
			break;
			case 3: {
				for (int i = 0; i < TestBlock::SIZE / 8; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS64(0);
					//HANDLE_BITS(0,bits0);
					handle_high_levels(0, bits0);
				}
			}
			break;
			default:
				issue_error();
		}
#undef HANDLE_BITS
		blocks_tested += 1;
		data += 1;
		numblocks -= 1;
		if (!numblocks) return;
	}
}
void PractRand::Tests::BCFN_MT::get_results(std::vector<TestResult> &results) {
		//total weight: varies, generally 0.5 on short samples to 1.5 on long samples
	return;
	/*if (!blocks_tested) return;
	//results.push_back(TestResult(this->get_name() + ":!", autofail ? 1 : 0, autofail ? 1 : 0, TestResult::TYPE_PASSFAIL, 0.000001));
	static const double chance_skipped[15] = {
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
	const double ref_chance = 1 - chance_skipped[5];
	std::vector<Uint64> tempcount; tempcount.resize(1<<tbits);
	std::vector<double> probs; probs.resize(1<<tbits);
	for (int level = 0; level < LEVELS; level++) {
		counts[level].flush();
		for (int i = 0; i < tempcount.size(); i++) tempcount[i] = counts[level][i];
		double p = 1.0 / probs.size();
		for (int i = 0; i < probs.size(); i++) probs[i] = p;

		int nlevel = level + unitsL2 + 3;
		double chance = chance_skipped[(nlevel < 15) ? nlevel : 14];
		if (nlevel > 14) chance *= std::pow(0.5, (nlevel-14) * 0.5);
		chance = 1 - chance;
		double samples = blocks_tested * TestBlock::SIZE * pow(0.5, level+unitsL2) * chance - tbits + 1;
		if (samples < 0) samples = 0;
		int effective_bits;
		//	when tbits is 1, the result is a close aproximation of a normal distribution
		//	but at higher tbits, the samples are coorelated and the deviation from normality becomes large
		//	once the number of samples becomes large though, it stabilizes to a consistent pattern
		//	if the number of samples is not large enough to reach a stable pattern then I transform the data to make it act like tbits was smaller
		if (false) ;
		else if (samples >= 256 * 1024 * 1024) effective_bits = 15;
		else if (samples >= 128 * 1024 * 1024) effective_bits = 14;
		else if (samples >= 64 * 1024 * 1024) effective_bits = 13;
		else if (samples >= 16 * 1024 * 1024) effective_bits = 12;
		else if (samples >= 8 * 1024 * 1024) effective_bits = 11;
		else if (samples >= 2 * 1024 * 1024) effective_bits = 10;
		else if (samples >= 1 * 1024 * 1024) effective_bits = 9;
		else if (samples >= 256 * 1024) effective_bits = 8;
		else if (samples >= 64 * 1024) effective_bits = 7;
		else if (samples >= 32 * 1024) effective_bits = 6;
		else if (samples >= 8 * 1024) effective_bits = 5;
		else if (samples >= 1536) effective_bits = 4;
		else if (level) continue; else effective_bits = 4;
		if (samples < 16 * 1024 && level) continue;
		if (effective_bits > tbits) effective_bits = tbits;

		TestCalibrationData *calib = NULL;
		if (effective_bits == 4) {
			if (samples > 1<<13) calib = calibration_manager.get_calibration_data("BCFN-4/4", samples / ref_chance / (1024/4));
			else calib = calibration_manager.get_calibration_data("_BCFN-4/4", samples / ref_chance / (1024/4));
		}
		else if (effective_bits > 4) {
			std::ostringstream internal_name;
			internal_name << "_BCFN-" << effective_bits << "/4";
			calib = calibration_manager.get_calibration_data(internal_name.str(), samples / ref_chance / (1024/4));
		}
		
		std::ostringstream name;
		name << "BCFN(";
		//double w = 0.5 / (1.0 + level * 0.25 + level * level * 0.05);
		double w;
		if (unitsL2 + level < 5) w = (unitsL2 + level + 2) / 8.0;
		else w = std::pow(0.5, double((unitsL2 + level + 2)/2));
		if (w > 0.375) w = 0.375;
		if (!level) w += 0.125;

		if (calib) {
			name << unitsL2 << "+" << level << "," << tbits << "-" << (tbits - effective_bits) << ")";
			truncate_table_bits(&tempcount[0], &probs[0], tbits, effective_bits);
			double rv = g_test(1 << effective_bits, &probs[0], &tempcount[0]);
			double rn = math_chisquared_to_normal(rv, (1<<effective_bits)-1);
			double rs = calib->sample_to_suspicion(rn) * -1;
			results.push_back(TestResult(name.str(), rn, rs, TestResult::TYPE_GOOD_S, w ) );
		}
		else {
			name << unitsL2 << "+" << level << "," << tbits << ")";
			double rv = g_test(1 << tbits, &probs[0], &tempcount[0]);
			double rn = math_chisquared_to_normal(rv, (1<<tbits)-1);
			results.push_back(TestResult(name.str(), rn, std::fabs(rn) > 25, TestResult::TYPE_PASSFAIL, w ) );
		}
	}*/
}
//*/
#endif















PractRand::Tests::BCFN_FF::BCFN_FF( int unitsL2_, int tbits_, bool unbalanced_ ) {
	unitsL2 = unitsL2_;
	tbits = tbits_;
	unbalanced = unbalanced_;
}
void PractRand::Tests::BCFN_FF::init( PractRand::RNGs::vRNG *known_good ) {
	int tsize = 1 << tbits;
	mask = tsize - 1;
	for (int level = 0; level < LEVELS; level++) {
		leftovers[level] = 0;
		warmup[level] = tbits - 1;
		cur[level] = 0;
		even[level] = false;
		counts[level].set_size(tsize);
		counts[level].reset_counts();

		int tmp = level + unitsL2 + 3 - 9;//2048-4096 bits should be shift=1 for COUNTS2_SIZE==256
		shifts[level] = tmp > 0 ? tmp/2 : 0;
		for (int x = 0; x < COUNTS2_SIZE; x++) counts2[level][x] = 0;
		extreme_counts2[level].clear();
	}
	blocks_tested = 0;
}
void PractRand::Tests::BCFN_FF::deinit() {
	for (int i = 0; i < LEVELS; i++) {
		counts[i].reset_counts();
	}
}
std::string PractRand::Tests::BCFN_FF::get_name() const {
	std::ostringstream f;
	f << "BCFN_FF(" << unitsL2 << "," << tbits << ")";
	return f.str();
	//return make_string("BCFN-%d", tbits);
}
void PractRand::Tests::BCFN_FF::get_results(std::vector<TestResult> &results) {
		//total weight: varies, generally 0.5 on short samples to 1.5 on long samples
	if (!blocks_tested) return;
	//results.push_back(TestResult(this->get_name() + ":!", autofail ? 1 : 0, autofail ? 1 : 0, TestResult::TYPE_PASSFAIL, 0.000001));
	static const double chance_skipped[15] = {
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
	const double ref_chance = 1 - chance_skipped[5];

	const double ref_chance_unbalanced = 1 - calculate_center_bit_combination_chance(5);
	std::vector<Uint64> tempcount; tempcount.resize(1 << tbits);
	std::vector<double> probs; probs.resize(1 << tbits);
	double overall_raw = 0;
	double overall_bins = 0;
	for (int level = 0; level < LEVELS; level++) {
		int nlevel = level + unitsL2 + 3;
		double chance_balanced = calculate_center_bit_combination_chance(nlevel);
		double samples = blocks_tested * TestBlock::SIZE * pow(0.5, level + unitsL2);
		if (!unbalanced) samples *= 1 - chance_balanced;
		samples -= tbits + 1;
		if (samples < 0) samples = 0;
		int effective_bits;
		/*
		when tbits is 1, the result is a close aproximation of a normal distribution
		but at higher tbits, the samples are coorelated and the deviation from normality becomes large
		once the number of samples becomes large though, it stabilizes to a consistent pattern
		if the number of samples is not large enough to reach a stable pattern then I transform the data to make it act like tbits was smaller
		*/
		double safety_factor = 1;
		static const float unbalanced_safety_table[6] = { 8, 6, 4.5, 3, 2, 1.5 };
		if (unbalanced) safety_factor = unbalanced_safety_table[nlevel >= 9 ? 5 : nlevel - 3];
		double adjusted_samples = samples / safety_factor;
		if (false);
		else if (adjusted_samples >= 256 * 1024 * 1024) effective_bits = 15;
		else if (adjusted_samples >= 128 * 1024 * 1024) effective_bits = 14;
		else if (adjusted_samples >= 64 * 1024 * 1024) effective_bits = 13;
		else if (adjusted_samples >= 16 * 1024 * 1024) effective_bits = 12;
		else if (adjusted_samples >= 8 * 1024 * 1024) effective_bits = 11;
		else if (adjusted_samples >= 2 * 1024 * 1024) effective_bits = 10;
		else if (adjusted_samples >= 1 * 1024 * 1024) effective_bits = 9;
		else if (adjusted_samples >= 256 * 1024) effective_bits = 8;
		else if (adjusted_samples >= 64 * 1024) effective_bits = 7;
		else if (adjusted_samples >= 32 * 1024) effective_bits = 6;
		else if (adjusted_samples >= 8 * 1024) effective_bits = 5;
		else if (adjusted_samples >= 1536) effective_bits = 4;
		else if (level || samples < 768) continue; else effective_bits = 4;
		if (effective_bits > tbits) effective_bits = tbits;

		TestCalibrationData *calib = NULL;
		if (effective_bits == 4) {
			if (adjusted_samples > 1 << 13) calib = calibration_manager.get_calibration_data("BCFN-4/4", adjusted_samples / ref_chance_unbalanced / (1024 / 4));
			else calib = calibration_manager.get_calibration_data("_BCFN-4/4", adjusted_samples / ref_chance_unbalanced / (1024 / 4));
		}
		else if (effective_bits > 4) {
			std::ostringstream internal_name;
			internal_name << "_BCFN-" << effective_bits << "/4";
			calib = calibration_manager.get_calibration_data(internal_name.str(), adjusted_samples / ref_chance_unbalanced / (1024 / 4));
		}

		std::ostringstream name;
		name << "BCFN_FF(";
		//double w = 0.5 / (1.0 + level * 0.25 + level * level * 0.05);
		double w;
		if (unitsL2 + level < 5) w = (unitsL2 + level + 2) / 8.0;
		else w = std::pow(0.5, double((unitsL2 + level + 2) / 2));
		if (w > 0.375) w = 0.375;
		if (!level) w += 0.125;
		/*
		unitsL2+level		Base Weight
		0					3/8
		1					3/8
		2*					3/8
		3					3/8
		4					2/8
		5					1/8
		6					1/16
		7					1/16
		8					1/32
		9					1/32
		10					1/64
		11					1/64
		12					1/128
		13					1/128
		...					...
		*/
		counts[level].flush();
		for (int i = 0; i < tempcount.size(); i++) tempcount[i] = counts[level][i];
		if (unbalanced) {
			double p1 = 0.5 - 0.5 * chance_balanced;
			double p1L = std::log(p1), p0L = std::log(1 - p1);
			for (int i = 0; i < probs.size(); i++) {
				int b = count_ones32(i);
				probs[i] = std::exp(p1L * b + p0L * (tbits - b));
			}
		}
		else {
			double p = 1.0 / probs.size();
			for (int i = 0; i < probs.size(); i++) probs[i] = p;
		}

		if (calib) {
			name << unitsL2 << "+" << level << "," << tbits << "-" << (tbits - effective_bits) << (unbalanced ? ",T" : ",F") << ")";
			truncate_table_bits(&tempcount[0], &probs[0], tbits, effective_bits);
			double rv = g_test(1 << effective_bits, &probs[0], &tempcount[0]);
			//for (int i = 0; i < 1<<effective_bits; i++) if (tempcount[i]) overall_raw += tempcount[i] * std::log(tempcount[i] / (probs[i] * samples));
			//overall_bins += (1 << effective_bits) - 1;
			double rn = math_chisquared_to_normal(rv, (1 << effective_bits) - 1);
			double rs = calib->sample_to_suspicion(rn) * -1;
			results.push_back(TestResult(name.str(), rn, rs, TestResult::TYPE_GOOD_S, w));
		}
		else {
			name << unitsL2 << "+" << level << "," << tbits << ")";
			double rv = g_test(1 << tbits, &probs[0], &tempcount[0]);
			double rn = math_chisquared_to_normal(rv, (1 << tbits) - 1);
			results.push_back(TestResult(name.str(), rn, std::fabs(rn) > 25 + (samples > 12345 ? 0 : 5), TestResult::TYPE_PASSFAIL, w));
		}
	}
	if (overall_bins) {
		double overall_norm = math_chisquared_to_normal(overall_raw * 2, overall_bins);
		results.push_back(TestResult(get_name() + ":all", overall_norm, overall_norm, TestResult::TYPE_RAW_NORMAL, 0.1));
	}

	for (int level = 0; level < LEVELS; level++) {
		//now the counts2 stuff:
		Uint64 total = 0;
		for (int i = 0; i < COUNTS2_SIZE; i++) total += counts2[level][i];
		total += extreme_counts2[level].size();
		if (!total) continue;
		int num_bits_L2 = unitsL2 + level + 3;
		int num_bits = 1 << num_bits_L2;
		int n = num_bits >> 1;
		if (num_bits_L2 <= 16) {
			if (total > 64) {
				std::vector<double> pdf, cdf;
				get_hamming_weight_chances(1 << num_bits_L2, pdf, cdf);
				double probs[COUNTS2_SIZE+2];
				Uint64 counts2_dup[COUNTS2_SIZE+2];
				for (int i = 0; i < COUNTS2_SIZE; i++) {
					counts2_dup[i+1] = counts2[level][i];
					int min = (i - COUNTS2_SIZE/2) << shifts[level];
					int max = min + (1 << shifts[level]) - 1;
					if (max < -n || min > n) probs[i+1] = 0;
					else {
						double p = 0;
						for (int j = min; j <= max; j++) p += abs(j) <= n ? pdf[n-abs(j)] : 0;
						probs[i+1] = p;
					}
				}
				/*if (COUNTS2_SIZE << shifts[level] >= num_bits) {
					probs[0] = 0;
					counts2_dup[0] = 0;
					probs[COUNTS2_SIZE+1] = 0;
					counts2_dup[COUNTS2_SIZE+1] = 0;
				}
				else {
					probs[0] = pdf[((-COUNTS2_SIZE/2) << shifts[level]) + n - 1];
					probs[COUNTS2_SIZE+1] = pdf[((-COUNTS2_SIZE/2) << shifts[level]) + n - 1];
					Uint64 lc = 0, hc = 0;
					for (
				}*/
				double samples = blocks_tested * TestBlock::SIZE * pow(0.5, level+unitsL2) - tbits + 1;
				int reduced_size = simplify_prob_table(COUNTS2_SIZE, samples / 40., &probs[1], &counts2_dup[1], true, false);
				double tr = g_test(reduced_size, &probs[1], &counts2_dup[1]);
				double n = math_chisquared_to_normal(tr, reduced_size-1);
				double p = math_chisquared_to_pvalue(tr, reduced_size - 1);
				double p2 = math_normaldist_to_pvalue(-n / 1.6);
				std::ostringstream name;
				name << "BCFN_FF(" << unitsL2 << "+" << level << "):freq";
				double w;
				if (unitsL2 + level < 5) w = (unitsL2 + level + 2) / 8.0;
				else w = std::pow(0.5, double((unitsL2 + level + 2)/2));
				if (w > 0.375) w = 0.375;
				if (!level) w += 0.125;
				//results.push_back(TestResult(name.str(), n, n, TestResult::TYPE_RAW_NORMAL, w*0.1 ) );
				results.push_back(TestResult(name.str(), n, p2, TestResult::TYPE_BAD_P, w*0.05 ) );
				//results.push_back(TestResult(name.str(), n, p, TestResult::TYPE_GOOD_P, w*0.1 ) );
			}
		}
	}
}
void PractRand::Tests::BCFN_FF::handle_high_levels ( int level, int bits ) {
	if (level >= LEVELS) return;
	if (!even[level]) {
		even[level] = true;
		leftovers[level] = bits;
		return;
	}
	else {
		even[level] = false;
		bits += leftovers[level];
	}

	if (true) {
		unsigned int c2i = (bits >> shifts[level]) + COUNTS2_SIZE/2;
		if (c2i < COUNTS2_SIZE) counts2[level][c2i]++;
		else extreme_counts2[level].push_back(bits);
	}
	if (true) {
		int tmp = bits >> 31;
		int index = ((cur[level] << 1) - tmp) & mask;
		cur[level] = index;
		if (warmup[level]) warmup[level] --;
		else {
			counts[level].increment(index);
		}
	}
	handle_high_levels(level+1, bits);
}
void PractRand::Tests::BCFN_FF::test_blocks(TestBlock *data, int numblocks) {
	while (blocks_tested < 16 && numblocks > 1) {
		test_blocks(data, 1);
		data += 1;
		numblocks -= 1;
	}
	while (warmup[4]) {
#define GET_BITS8(pos)  (count_ones8 (data[0].as8 [i+(pos)]) - 4)
#define GET_BITS16(pos) (count_ones16(data[0].as16[i+(pos)]) - 8)
#define GET_BITS32(pos) (count_ones32(data[0].as32[i+(pos)]) - 16)
#define GET_BITS64(pos) (count_ones64(data[0].as64[i+(pos)]) - 32)
#define HANDLE_BITS(level,var) {counts2[level][var+COUNTS2_SIZE/2]++; if (true){tmp=var>>31;cur[level]=((cur[level]<<1)-tmp)&mask;if (warmup[level]) warmup[level]--; else counts[level].increment(cur[level]);}}
		switch (unitsL2) {
			case 0: {
				for (int i = 0; i < TestBlock::SIZE / 1; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS8(0);
					HANDLE_BITS(0,bits0);
					handle_high_levels(1, bits0);
				}
			}
			break;
			case 1: {
				for (int i = 0; i < TestBlock::SIZE / 2; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS16(0);
					HANDLE_BITS(0,bits0);
					handle_high_levels(1, bits0);
				}
			}
			break;
			case 2: {
				for (int i = 0; i < TestBlock::SIZE / 4; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS32(0);
					HANDLE_BITS(0,bits0);
					handle_high_levels(1, bits0);
				}
			}
			break;
			case 3: {
				for (int i = 0; i < TestBlock::SIZE / 8; i+=1) {
					int bits0;
					int tmp;
					//0
					bits0 = GET_BITS64(0);
					HANDLE_BITS(0,bits0);
					handle_high_levels(1, bits0);
				}
			}
			break;
			default:
				issue_error();
		}
#undef HANDLE_BITS
		blocks_tested += 1;
		data += 1;
		numblocks -= 1;
		if (!numblocks) return;
	}
	if (true) {
		unsigned long max = numblocks * TestBlock::SIZE >> unitsL2;
		switch (unitsL2) {
#define GET_BITS(a) GET_BITS8(a)
			case 0: {
				for (unsigned long i = 0; i < max; i+=8) {
					long bits0, bits1, bits2, bits3;
					long tmp;
#define HANDLE_BITS(level,var) {counts2[level][var+COUNTS2_SIZE/2]++; if (true){tmp=var>>31;cur[level]=((cur[level]<<1)-tmp)&mask;if (warmup[level]) warmup[level]--; else counts[level].increment(cur[level]);}}
					//0
					bits1 = bits0 = GET_BITS(0);
					HANDLE_BITS(0,bits0);
					//1
					bits2 = bits1 += bits0 = GET_BITS(1);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//2
					bits1 = bits0 = GET_BITS(2);
					HANDLE_BITS(0,bits0);
					//3
					bits3 = bits2 += bits1 += bits0 = GET_BITS(3);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					//4
					bits1 = bits0 = GET_BITS(4);
					HANDLE_BITS(0,bits0);
					//5
					bits2 = bits1 += bits0 = GET_BITS(5);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//6
					bits1 = bits0 = GET_BITS(6);
					HANDLE_BITS(0,bits0);
					//7
					bits3 += bits2 += bits1 += bits0 = GET_BITS(7);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					HANDLE_BITS(3,bits3);
					handle_high_levels ( 4, bits3 );
				}
			}
			break;
#undef GET_BITS
#define GET_BITS(a) GET_BITS16(a)
			case 1: {
				for (unsigned long i = 0; i < max; i+=8) {
					long bits0, bits1, bits2, bits3;
					long tmp;
					//0
					bits1 = bits0 = GET_BITS(0);
					HANDLE_BITS(0,bits0);
					//1
					bits2 = bits1 += bits0 = GET_BITS(1);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//2
					bits1 = bits0 = GET_BITS(2);
					HANDLE_BITS(0,bits0);
					//3
					bits3 = bits2 += bits1 += bits0 = GET_BITS(3);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					//4
					bits1 = bits0 = GET_BITS(4);
					HANDLE_BITS(0,bits0);
					//5
					bits2 = bits1 += bits0 = GET_BITS(5);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//6
					bits1 = bits0 = GET_BITS(6);
					HANDLE_BITS(0,bits0);
					//7
					bits3 += bits2 += bits1 += bits0 = GET_BITS(7);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					HANDLE_BITS(3,bits3);
					handle_high_levels ( 4, bits3 );
				}
			}
			break;
#undef GET_BITS
#define GET_BITS(a) GET_BITS32(a)
			case 2: {
				for (unsigned long i = 0; i < max; i+=8) {
					long bits0, bits1, bits2, bits3;
					long tmp;
					//0
					bits1 = bits0 = GET_BITS(0);
					HANDLE_BITS(0,bits0);
					//1
					bits2 = bits1 += bits0 = GET_BITS(1);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//2
					bits1 = bits0 = GET_BITS(2);
					HANDLE_BITS(0,bits0);
					//3
					bits3 = bits2 += bits1 += bits0 = GET_BITS(3);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					//4
					bits1 = bits0 = GET_BITS(4);
					HANDLE_BITS(0,bits0);
					//5
					bits2 = bits1 += bits0 = GET_BITS(5);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//6
					bits1 = bits0 = GET_BITS(6);
					HANDLE_BITS(0,bits0);
					//7
					bits3 += bits2 += bits1 += bits0 = GET_BITS(7);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					HANDLE_BITS(3,bits3);
					handle_high_levels ( 4, bits3 );
				}
			}
			break;
#undef GET_BITS
#define GET_BITS(a) GET_BITS64(a)
			case 3: {
				for (unsigned long i = 0; i < max; i+=8) {
					int bits0, bits1, bits2, bits3;
					int tmp;
					//0
					bits1 = bits0 = GET_BITS(0);
					HANDLE_BITS(0,bits0);
					//1
					bits2 = bits1 += bits0 = GET_BITS(1);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//2
					bits1 = bits0 = GET_BITS(2);
					HANDLE_BITS(0,bits0);
					//3
					bits3 = bits2 += bits1 += bits0 = GET_BITS(3);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					handle_high_levels ( 3, bits2 );//bits2 is already 256 bits, any more might overflow the array index on counts2[3]
					//4
					/*bits1 = bits0 = GET_BITS(4);
					HANDLE_BITS(0,bits0);
					//5
					bits2 = bits1 += bits0 = GET_BITS(5);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					//6
					bits1 = bits0 = GET_BITS(6);
					HANDLE_BITS(0,bits0);
					//7
					bits3 += bits2 += bits1 += bits0 = GET_BITS(7);
					HANDLE_BITS(0,bits0);
					HANDLE_BITS(1,bits1);
					HANDLE_BITS(2,bits2);
					HANDLE_BITS(3,bits3);
					handle_high_levels ( 4, bits3 );*/
				}
			}
			break;
#undef GET_BITS
#undef GET_BITS8
#undef GET_BITS16
#undef GET_BITS32
#undef GET_BITS64
#undef HANDLE_BITS
			default:
				issue_error();
		}
	}
	blocks_tested += numblocks;
}








PractRand::Tests::FPMulti::FPMulti() //(int stride_bits_L2_, int skip_platters_)
	//:
	//stride_bits_L2(stride_bits_L2_)
{
	static const double gap_L2_expected_lookup[30] = {
		// I printed lots of precision, but the actual calculations weren't that accurate.  I'm pretty sure it's all at least as good as single-precision though, probably better.  
		0.732649482117484,
		1.537438290932739,
		2.401606811975667,
		3.311224720400715,
		4.253426596472763,
		5.217705249861148,
		6.196250654101414,
		7.183665553491496,
		8.176424757912608,
		9.172324308194467,
		10.170032291922604,
		11.168764874403314,
		12.168070314222019,
		13.167692567126194,
		14.167488448594183,
		15.167378763675536,
		16.167320107527143,
		17.167288872374503,
		18.167272301195908,
		19.167263538788724,
		20.167258919171243,
		21.167256490152202,
		22.167255216034651,
		23.167254549165989,
		24.167254200815417,
		25.167254019159184,
		26.167253924545253,
		27.167253875254346,
		28.167253848788377,
		29.167253834281478,
	};
	for (int i = 0; i <= MAX_EXP; i++) {
		int L = i + BASE_SIG_BITS + 1 - (i == MAX_EXP ? 1 : 0);
		double e;
		enum { THRESHOLD = 30 };
		//if (L <= THRESHOLD) e = gap_log2_expected(std::pow(0.5, L));
		//else e = gap_log2_expected(std::pow(0.5, 20)) + L - THRESHOLD;
		if (L < THRESHOLD) e = gap_L2_expected_lookup[L - 1];
		else e = gap_L2_expected_lookup[THRESHOLD - 1] + L - THRESHOLD;
		//std::printf("%d %.15f\n", i, e);
		platter[i].gap_expected_inverse = std::pow(0.5, e);
		//platter[i].gap_expected_inverse = std::pow(0.5, L-1);
	}
}
void PractRand::Tests::FPMulti::Platter::reset(PractRand::RNGs::vRNG *known_good, unsigned long e) {
	total_count = 0;

	//gap test stuff:
	for (int i = 0; i < (1 << GAP_SIG_BITS); i++) {
		//gap_global_history[i] = 1ull << 63;
		gap_global_history[i] = 0;
		//gap_global_history[i] = 0 - known_good->randli(1ull << (e + GAP_SIG_BITS + 1));
		//gap_global_history[i] = 0 - 1ull << (e + GAP_SIG_BITS);
		gap_first_hit[i] = 0;
	}
	gap_product = 0.5;
	gap_product_extracted_L2 = 1;
	gap_hits = 0;
	gap_unique_sigs = 0;
	//gap_warmed_up = 0;//
	
	//for (int i = 0; i < COUP_MASK_SIZE; i++) coup_mask[i] = 0;
	//freq_count.reset_counts();
	//coup_count.reset_counts();
	//last_coup = 0;
}
void PractRand::Tests::FPMulti::process(Uint64 position, unsigned long e, unsigned long sig) {
	Platter &p = platter[e];
	p.total_count++;

	if (1) {
		unsigned long gap_sig = sig >> (BASE_SIG_BITS - GAP_SIG_BITS);
		Uint64 old_pos = p.gap_global_history[gap_sig];
		p.gap_global_history[gap_sig] = position;// & 0x7FffFFffFFffFFffull;
		//enum { WARMUP_SETS = 2 };//
		//if (!p.gap_warmed_up) {
		//	Uint64 warmup_distance = WARMUP_SETS << (e + GAP_SIG_BITS);
		//	if (position >= warmup_distance) p.gap_warmed_up = true;
		//}
		//if (p.gap_warmed_up) {
		if (true) {
			if (!old_pos) {
				/*
					ideas for how to handle samples with no preceding value:
						1. don't count them
								NO: throws off end result too much, or increases data requirements
						2. count them as having a preceding value at position zero
								NO: still throws off end results too much, or increases data requirements
						3. multiply their position by a fudge factor, treat result as a gap value
								NO: though it works a bit better
						4. add the mean (or median?) gap value to their position, treat that as a gap value
								MAYBE: better theoretical basis, needs more testing
						5. generate a random value on a distribution similar to gap value, add it to their position
								MAYBE: sounds really sleazy, generates unwanted noise, but probably better than many alternatives
						6. generate random values on distriubtion like gap values, discard those less than position
								MAYBE: should be equivalent to #5, I think?
						7. generate many values spread out on the distribution like gap values, replace anything over the position with the position, take the logs, normalize, measure the variance of the resulting curve, add it to the expected variance of all samples, add the log of the position to the cumulative product L2
								EHHHHH: a high chance I screwed up my logic in there somewhere and wouldn't figure out where until spending way too long on it, if ever
						8. add the distance *after* the last occurance to the distance before the first occurance... the sum should follow the correct distribution, right?  but what if there is zero occurances... use another method as a fallback path?
								YEAH: I'm feeling good about this, so far, though it needs fleshing out - ideally checking the tail and the start against the same distribution would be better
								no... this doesn't help when only short gaps are possible due to insufficient data relative to platter level... does it?  
						9. figured out exact distribution based upon the idea that distance before first occurance had the same distribution as distance after last occurance, and the sum of the two was equal to the distribution of distances between two consecutive occurances, therefore the desired distribution convolved with itself must equal the distribution between consecutive occurances
								K = (chance of the target symbol occuring), x = (1 - K), M = sqrt(K)
								base gap distribution: {K * x^0, K * x^1, K * x^2, K * x^3, K * x^3, ...}
								half-gap distribution: {M * x^0 * n[0], M * x^1 * n[1], M * x^2 * n[2], M * x^3 * n[3], ...}
								n[0] = 1.0 ; n[i] = n[i-1] * (i * 2 - 1) / (i * 2)
								Anyway, given the exact distribution, I can make an ordered normalized scoring system for them, create a normally distributed overall score for all preceding semi-gaps and postceding semi-gaps in a given platter, convert that to a p-value
								Though, is that really any better than combining preceding and postceding semi-gaps to a single gap?
								INTERESTING: but... I think that's only true for large sample sizes ; so... it doesn't help anything that truly needs help, only makes the easier cases a tiny bit better quantified
								instead, need a way to handle short sequences systematically and well, regardless of number of occurances within a sequence
						10. maybe go across significand values?  like... the tail of sig0 to the pre of sig1 to make a gap distributed sample - mostly worthless, no better than #8 above ; but... if there are *zero* samples at a given sig, that length could be added to the preceding sigs tail and following sig pre, and if they have zero too then you can just keep going
								the results should be very close to gap distributed, and work even in cases where the expected per-sig occurrance rate is very low, so long as the expected whole-platter occurances is at least... 20 or so - beyond that, just return a '!' subtest result for the chance of there being zero occurances
								...maybe I should assume demand a minimum of 3 actual occurances?  no, higher... 10? that requires a much higher number of expected occurances though... actually... I think 120 expected should be enough even for a simple autofail on less than 10 actual, 40 might be enough with more sophisticated reporting
								still, with 9 bit significands, 120 expected occurances means 0.23 expected per significand value, which is a big improvement over the regular way
								actually, no, just key the whole thing off the actual count not the expected ; tests for deviation between the two should be handled seperately as part of the cross-platter-frequency-test instead of anything difficult
E	A=0		A=1		A=2		A=3		A=4		...		A=10
10	5e-5	5e-4	3e-3	8e-3	2e-2			1e-1
20	2e-9	4e-8	5e-7	3e-6	1e-5			6e-3
30	9e-14	3e-12	4e-11	4e-10	3e-9			2e-5
40	4e-18	2e-16	3e-15	5e-14	5e-13			1e-8
50	2e-22	1e-20	2e-19	4e-18	5e-17			5e-12
60	9e-27	5e-25	2e-23	3e-22	5e-21			2e-15

E=40 A=10: 1e-8
E=60 A=10: 2e-15
E=90 A=10: 8e-27

.999^10k
+ 10k * .999^(10k-1) * 0.001
+ 10k * (10k-1) / 2 * .999^(10k-2) * 0.001^2
+ 10k * (10k-1) * (10k-2) / 6 * 0.999^(10k-3) * 0.001^3





				*/
				p.gap_first_hit[gap_sig] = position;
				p.gap_unique_sigs++;

				/*PractRand::Tests::SampleSet ss;
				enum {I = 256};
				for (int i = 0; i < I; i++) {
					double f = (i + 0.5) / I;
					double random_gap = std::ceil((1ull << (e - (e == MAX_EXP) ? 1 : 0)) * -std::log(f));
					ss.add(std::log2((random_gap + position) * p.gap_expected_inverse));
				}
				p.gap_first_hits_total_variance += ss.get_variance();
				p.gap_first_hits += 1;
				p.gap_first_hits_example_product *= std::ceil((1ull << (e - (e == MAX_EXP) ? 1 : 0)) * -std::log(internal_rng->randlf(0.001, 0.999))) * p.gap_expected_inverse;

				//p.gap_warmed_up += 1;//
				//Uint64 warmup_distance = WARMUP_SETS << (e + GAP_SIG_BITS);
				//g += (g - warmup_distance) * 1.5;
				//g += std::pow(2.0, e + GAP_SIG_BITS); // * 0.63212056;//
				*/
			}
			else {
				p.gap_hits += 1;
				double g = (position - old_pos);// & 0x7FffFFffFFffFFffull;
				double normalized = g * p.gap_expected_inverse;
				p.gap_product *= normalized;
				if (!(p.gap_hits & 63)) {//testing suggests this is good up to 4095 - for safety margin I use 1023, and set autofail on any gap products out of range at that point
					// no wait... if GAP_SIG_BITS is adjusted the usable range changes (4095 was for 9 bits), though 127 seems to be usable for all useful values of GAP_SIG_BITS
					int L2;
					if (std::isinf(p.gap_product)) issue_error("FPMulti::process - gap product is infinite");
					if (std::isinf(p.gap_product)) autofail = true;
					if (0 == p.gap_product) issue_error("FPMulti::process - gap product is zero");
					if (p.gap_product == 0) autofail = true;
					if (std::isnan(p.gap_product)) issue_error("FPMulti::process - gap product is NaN");// should be impossible, I think
					if (autofail) return;
					p.gap_product = std::frexp(p.gap_product, &L2);
					p.gap_product_extracted_L2 += L2;
				}
			}
		}
	}
	/*
	if (1) {//frequency
		p.freq_count.increment(sig >> (BASE_SIG_BITS - GAP_SIG_BITS));
	}
	if (1) {//coupon
		unsigned long index = sig & ((1 << COUP_SIG_BITS) - 1);
		COUP_WORD b = COUP_WORD(1) << (index & (COUP_WORD_SIZE - 1));
		index >>= COUP_WORD_SIZE_L2;
		p.coup_mask[index] |= b;
		if (p.coup_mask[index] == ~COUP_WORD(0)) {//possible end of set
			bool done = true;
			for (int i = 0; i < COUP_MASK_SIZE; i++) if (~p.coup_mask[i]) done = false;
			if (done) {
				for (int i = 0; i < COUP_MASK_SIZE; i++) p.coup_mask[i] = 0;
				Sint64 delta = p.total_count - p.last_coup;
				p.last_coup = p.total_count;
				delta -= 1 << COUP_SIG_BITS;
				if (delta < 0) issue_error("FPM:C: impossible");
				index = delta >> (COUP_SIG_BITS - 5);
				if (index > 1023) autofail = true;
			}
		}
	}*/
}
void PractRand::Tests::FPMulti::init(RNGs::vRNG *known_good) {
	TestBaseclass::init(known_good);
	autofail = false;

	for (int e = 0; e <= MAX_EXP; e++) platter[e].reset(known_good, e);

	internal_rng = new PractRand::RNGs::Polymorphic::arbee(known_good);
}
void PractRand::Tests::FPMulti::deinit() {
	delete internal_rng;
	internal_rng = NULL;
}
std::string PractRand::Tests::FPMulti::get_name() const {
	//std::ostringstream str;
	//str << "FPM(" << skip_platters << "/" << (1 << stride_bits_L2) << ")";
	//return str.str();
	return "FPM";
}
void PractRand::Tests::FPMulti::get_results(std::vector<TestResult> &results) {
	Uint64 total_samples = blocks_tested * (TestBlock::SIZE / sizeof(Uint64));

	if (1) {// gap test preliminary work checking for autofail
		for (int e = 0; e <= MAX_EXP && !autofail; e++) {
			Platter &p = platter[e];
			//if (std::isinf(p.gap_product)) issue_error("FPMulti::get_results - gap product is infinite");
			if (std::isinf(p.gap_product)) { autofail = true; continue; }
			//if (0 == p.gap_product) issue_error("FPMulti::get_results - gap product is zero");
			if (p.gap_product == 0) { autofail = true; continue; }
			if (std::isnan(p.gap_product)) issue_error("FPMulti::get_results - gap product is NaN");// should be impossible, I think

			int L2;
			p.gap_product = std::frexp(p.gap_product, &L2);
			p.gap_product_extracted_L2 += L2;
		}
	}

	if (autofail) {
		results.push_back(TestResult(get_name() + ":!", autofail, autofail, TestResult::TYPE_PASSFAIL, 0.001));
		return;
	}

	if (1) {// full gap test
		enum { NUM_PRECALCED = 29 };
		static const double precalced_per_sample_variance[NUM_PRECALCED] = {
			// I printed lots of precision, but the actual calculations weren't that accurate.  I think it's all at least as good as single-precision though.  
			0.689767784941473,
			1.337738769110027,
			1.901334686734523,
			2.357736926112877,
			2.704552837299473,
			2.954032389076138,
			3.125391844474298,
			3.238662124118501,
			3.311200833954624,
			3.356456856372971,
			3.384086977313802,
			3.400654090320465,
			3.410437953635610,
			3.416141765833042,
			3.419430341441828,
			3.421308286318003,
			3.422371735456495,
			3.422969518546240,
			3.423303348278004,
			3.423488686362224,
			3.423591043978246,
			3.423647305527636,
			3.423678096892326,
			3.423694882477087,
			3.423704000023767,
			3.423708936075639,
			3.423711600191441,
			3.423713034019782,
			3.423713803681121
		};
		Uint64 total_gap_hits = 0;
		double total_gap_product_L2 = 0;
		double total_gap_product_L2_adjusted_variance = 0;
		Uint64 total_gap_first_hits = 0;
		double total_gap_first_hit_product_L2 = 0;
		double total_gap_first_hit_product_L2_adjusted_variance = 0;

		for (int e = 0; e <= MAX_EXP; e++) {
			int L = GAP_SIG_BITS + e - (e == MAX_EXP ? 1 : 0);
			Platter &p = platter[e];
			double gap_product_L2 = std::log2(p.gap_product);
			gap_product_L2 += p.gap_product_extracted_L2;
			if (p.gap_hits + p.gap_unique_sigs < 30) continue;

			//double L = 0 + e - (e == MAX_EXP ? 1 : 0);
			//double cLK = 0.7 - 0.8 / L + (4 + 32.0 / L) * std::pow(double(p.gap_hits), -3.0 / L) / 15;
			//double cLK = 0.7;
			//double cLK = 1.0;
			//double per_sample_variance = gap_log2_variance(std::pow(0.5, L));
			double per_sample_variance;
			if (L <= NUM_PRECALCED) per_sample_variance = precalced_per_sample_variance[L - 1];
			else per_sample_variance = precalced_per_sample_variance[NUM_PRECALCED - 1];
			//std::printf("%2d %2d %.15f\n", e, int(L), per_sample_variance);
			//double adjusted_variance = per_sample_variance * p.gap_hits * cLK * cLK;
			double adjusted_variance = per_sample_variance * (p.gap_hits * 1.0);// -p.gap_warmed_up * 0);

			total_gap_hits += p.gap_hits;
			total_gap_product_L2 += gap_product_L2;
			total_gap_product_L2_adjusted_variance += adjusted_variance;
			double avg = gap_product_L2 / p.gap_hits;
			double norm = gap_product_L2 / std::sqrt(adjusted_variance);
			if (adjusted_variance < 300) continue;

			std::ostringstream buf;
			buf << get_name() << ":G" << GAP_SIG_BITS << ":e" << e;
			results.push_back(TestResult(buf.str(), norm, math_normaldist_to_pvalue(norm), TestResult::TYPE_BAD_P, 0.02 * std::pow(0.65, e)));

			/*
				for Maurer's paper: (not that I'm actually following it that closely, but need these to be able to reference it quickly)
				Q = words skipped during warmup
				K = number of usable gap samples
				L = bits per word

				V = 2^L
				S^N = sequence of input bits
				S[n] = nth bit of input
				b[n] = nth word of input
				A[n](S^N) = gap value for nth word of input (distance since last occurance) - defined to be n if there is no gap prior occurance
				fTu(S^N) = log2 of geometric mean of all (usable) gaps
				c(L,K) = adjust factor that the standard deviation is multiplied by, probably fitted from empirical data
			*/
		}
		if (total_gap_product_L2_adjusted_variance > 5000) {
			double avg = total_gap_product_L2 / total_gap_hits;
			double norm = total_gap_product_L2 / std::sqrt(total_gap_product_L2_adjusted_variance);
			std::ostringstream buf;
			buf << get_name() << ":G" << GAP_SIG_BITS << ":comb";
			results.push_back(TestResult(buf.str(), norm, math_normaldist_to_pvalue(norm), TestResult::TYPE_BAD_P, 0.4));
		}
	}

	/*const int skip_platters = 0;
	for (int p = skip_platters; p <= MAX_EXP; p++) total_samples += platter[p].total_count;
	//double expected_total_samples2 = (blocks_tested * 8.0 * TestBlock::SIZE - 64) * std::pow(0.5, stride_bits_L2 + skip_platters) + 1;//not too sure about this calculation
	if (total_samples >= 16000 && true) {//platter local gap test results
		long double bins = 1 << GAP_SIG_BITS;
		long double base_prob = 1.0 / bins;
		long double base_prob_inv = 1.0 - base_prob;
		long double avg = 0;
		long double dev = 0;
		long double baseline = bins * bins - bins;
		for (int x = 64 << GAP_SIG_BITS; x >= 1; x--) {
			long double prob = std::pow(base_prob_inv, x - 1) / bins;
			double value = x - (1 << GAP_SIG_BITS);
			value *= value;
			value -= baseline;
			avg += value * prob;
			dev += value * value * prob;
		}
		dev = std::sqrt(dev - avg * avg);
		avg += baseline;
		for (int p = skip_platters; p <= MAX_EXP; p++) {
			//double expected_samples = 
			double samples = platter[p].total_count - platter[p].gap_negative_count1;
			if (samples < 8000) continue;
			double observed_avg = platter[p].gap_sum1 / samples;
			double diff = observed_avg - avg;
			double norm = diff * std::sqrt(samples) / dev;
			std::ostringstream name;
			name << get_name() << ":G1(" << p << ")";
			results.push_back(TestResult(name.str(), norm, norm, TestResult::TYPE_RAW_NORMAL, 0.02 * std::pow(0.65, p - skip_platters)));
		}
	}*/
	/*if (total_samples && true) {//frequency test results
		double freq_all_sum = 0;
		double freq_all_total = 0;
		double freq_all_bins = 0;
		double freq_all_tail_count = 0;
		double freq_all_tail_prob = 0;
		for (int p = skip_platters; p <= MAX_EXP; p++) {
			double expected_samples = total_samples * std::pow(0.5, p - skip_platters + 1 - (p == MAX_EXP));
			double samples = platter[p].total_count;
			int ebits = int(std::log(expected_samples) / std::log(2.0) * 0.5 - 0.4);
			if (ebits < 2 || !samples) {
				freq_all_tail_count += samples;
				freq_all_tail_prob += std::pow(0.5, p + 1 - (p == MAX_EXP) - skip_platters);
				continue;
			}
			std::vector<Uint64> counts_vec;
			const Uint64 *counts = platter[p].freq_count.get_array();
			if (ebits < FREQ_SIG_BITS) {
				counts_vec.resize(1 << FREQ_SIG_BITS);
				for (int i = 0; i < 1 << FREQ_SIG_BITS; i++) counts_vec[i] = counts[i];
				truncate_table_bits(&counts_vec[0], NULL, FREQ_SIG_BITS, ebits);
				counts = &counts_vec[0];
			}
			else ebits = FREQ_SIG_BITS;
			double freq_platter_sum = 0;
			for (int i = (1 << ebits) - 1; i >= 0; i--) if (counts[i]) freq_platter_sum += counts[i] * std::log(counts[i]);
			freq_all_bins += 1 << ebits;
			freq_all_total += samples;
			freq_all_sum += freq_platter_sum - samples * std::log(std::pow(0.5, p + 1 - (p == MAX_EXP) - skip_platters + ebits));
			if (samples < 200) continue;
			//double raw = PractRand::Tests::g_test_flat(1 << ebits, &counts[0]);
			//double norm = PractRand::Tests::math_chisquared_to_normal(raw, (1 << ebits) - 1);
			double raw2 = (freq_platter_sum - samples * std::log(samples * std::pow(0.5, ebits))) * 2.0;
			double norm2 = PractRand::Tests::math_chisquared_to_normal(raw2, (1 << ebits) - 1);
			std::ostringstream name;
			name << get_name() << ":F(" << p << "," << ebits << ")";
			//results.push_back(TestResult(name.str(), norm, norm, TestResult::TYPE_RAW_NORMAL, 0.02));
			results.push_back(TestResult(name.str(), norm2, norm2, TestResult::TYPE_RAW_NORMAL, 0.05 * std::pow(0.6, p - skip_platters)));
		}
		if (freq_all_total > 200) {//frequency, all platters combined
			if (freq_all_tail_prob) {
				freq_all_bins += 1;
				freq_all_total += freq_all_tail_count;
				if (freq_all_tail_count) freq_all_sum += freq_all_tail_count * std::log(freq_all_tail_count / freq_all_tail_prob);
			}
			freq_all_sum -= freq_all_total * std::log(freq_all_total);
			freq_all_sum *= 2;
			double all_norm = math_chisquared_to_normal(freq_all_sum, freq_all_bins - 1);
			double all_p = math_normaldist_to_pvalue(all_norm);
			TestCalibrationData *calib = NULL;//calibration_manager.get_calibration_data("FPF-14+6/16:overall", samples / 512.0 + 0.5);
			if (calib && total_samples >= 3000)
				results.push_back(TestResult(get_name() + ":F:all", all_norm, -calib->sample_to_suspicion(all_norm), TestResult::TYPE_GOOD_S, .25));
			else results.push_back(TestResult(get_name() + ":F:all", all_norm, all_norm, TestResult::TYPE_RAW_NORMAL, .25));
		}
		if (total_samples > 320) {//frequency, cross-platter
			std::vector<Uint64> cross_counts; cross_counts.resize(MAX_EXP + 1);
			for (int p = skip_platters; p <= MAX_EXP; p++) cross_counts[p] = platter[p].total_count;
			std::vector<double> cross_probs; cross_probs.resize(MAX_EXP + 1);
			for (int p = skip_platters; p <= MAX_EXP; p++) cross_probs[p] = std::pow(0.5, p + 1 - (p == MAX_EXP) - skip_platters);
			if (false) {
				for (int p = skip_platters; p <= MAX_EXP; p++) {
					std::printf("%2d:", p);
					if (cross_counts[p] >> 32) {
						std::printf("%8X", Uint32(cross_counts[p] >> 32));
						std::printf("%08X", Uint32(cross_counts[p]));
					}
					else std::printf("        %8X", Uint32(cross_counts[p]));
					std::printf("    %.15f\n", cross_probs[p]);
				}
			}
			//int bins = MAX_EXP + 1 - skip_platters;
			int bins = simplify_prob_table(MAX_EXP + 1 - skip_platters, total_samples / 10.0, &cross_probs[skip_platters], &cross_counts[skip_platters], true, true);
			if (false) {
				std::printf("---\n");
				for (int p = skip_platters; p < bins; p++) {
					std::printf("%2d:", p);
					if (cross_counts[p] >> 32) {
						std::printf("%8X", Uint32(cross_counts[p] >> 32));
						std::printf("%08X", Uint32(cross_counts[p]));
					}
					else std::printf("        %8X", Uint32(cross_counts[p]));
					std::printf("    %.15f\n", cross_probs[p]);
				}
			}
			double raw = PractRand::Tests::g_test(bins, &cross_probs[skip_platters], &cross_counts[skip_platters]);
			double norm = PractRand::Tests::math_chisquared_to_normal(raw, bins - 1);
			results.push_back(TestResult(get_name() + ":F:cross", norm, norm, TestResult::TYPE_RAW_NORMAL, 0.25));
		}
	}*/
}
void PractRand::Tests::FPMulti::test_blocks(TestBlock *data, int numblocks) {
	unsigned long end = numblocks << (TestBlock::SIZE_L2 - 3);
	Uint64 offset = (blocks_tested << (TestBlock::SIZE_L2 - 3)) + 1;
	//Uint64 *base_addr = &data[0].as64[-offset]; //we can optimize things slightly once we're more confident in this
	if (autofail) return;
	for (unsigned long i = 0; i < end; i++) {
		Uint64 raw = data[0].as64[i];
		unsigned long e = count_low_zeroes64(raw);
		unsigned long sig;
		if (e < MAX_EXP) {
			sig = (raw >> (e + 1)) & ((1 << BASE_SIG_BITS) - 1);
		}
		else {
			e = MAX_EXP;
			sig = (raw >> e) & ((1 << BASE_SIG_BITS) - 1);
		}
		process(i + offset, e, sig);
	}
#if 0
	//unsigned long stride_bits = 1 << stride_bits_L2;
	//Uint32 skip_mask = ((1 << skip_platters) - 1) << (32 - skip_platters);

	Uint32 skip_mask = (1 << skip_platters) - 1;
	base_position = blocks_tested * (8 * TestBlock::SIZE >> stride_bits_L2);

	if (stride_bits_L2 >= 5) {//long stride
		long stride32 = 1 << (stride_bits_L2 - 5);
		long max = numblocks * (TestBlock::SIZE / 4) - 1 - (stride32 - 1);
		for (long i = blocks_tested ? -1 : 0; i < max; i += stride32) {
			Uint32 cur = data->as32[i];
			if (cur & skip_mask) continue;
			unsigned long e = count_low_zeroes32(cur);
			unsigned long sig;
			if (e < 32 - BASE_SIG_BITS) process(platter[e], cur >> (e + 1),  i);
			else {
				//Uint64 cur2 = cur | (Uint64(reverse_bits32(data->as32[i + 1])) << 32);
				Uint64 cur2 = cur | (Uint64(data->as32[i + 1]) << 32);
				if (e == 32) e += count_low_zeroes32(Uint32(cur2 >> 32));
				if (e < MAX_EXP) process(platter[e], Uint32(cur2 >> (e + 1)), i);
				else process(platter[MAX_EXP], cur2 >> MAX_EXP, i);
			}
		}
	}
	else {//short stride
		long max = numblocks * (TestBlock::SIZE / 4) - 1;
		long start;
		Uint32 cur,next;
		if (blocks_tested) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
			cur = data->as32[-2];//64 bit samples on 32 bit words
			next = data->as32[-1];
#pragma GCC diagnostic pop
			start = -1;
		}
		else {
			cur = data->as32[0];//this misses the very first possible sample, but there's no way to do that with the shared code path
			next = data->as32[1];
			start = 1;
			//...so we have a special case just for the first possible sample (unoptimized)
			unsigned long e = count_low_zeroes32(cur);
			if (e == 32) e += count_low_zeroes32(next);
			if (e < MAX_EXP) process(platter[e], Uint32(data->as64[0] >> (e + 1)), 1);
			else process(platter[MAX_EXP], Uint32(data->as64[0] >> MAX_EXP), 1);
		}
		if (stride_bits_L2 == 4) {
			enum { STRIDE_BITS = 16 };
			for (long i = start; i < max; i ++) {
				Uint32 word = next;
				next = data->as32[i + 1];//we need 3 32-bit values to slide a 64 bit window across them

				unsigned long e, sig;

				cur >>= STRIDE_BITS;
				cur |= word << (32 - STRIDE_BITS);
				word >>= STRIDE_BITS;// it's easier to keep cur and word up to date as the window slides, but next rarely gets used
				if (!(cur & skip_mask)) {
					Uint32 n = (i << 1) + 0;
					e = count_low_zeroes32(cur);
					if (e < 32 - BASE_SIG_BITS) process(platter[e], cur >> (e + 1), n);
					else {
						Uint64 cur2 = cur | (Uint64(word) << 32) | (Uint64(next) << (64 - STRIDE_BITS));
						if (e == 32) e += count_low_zeroes32(Uint32(cur2 >> 32));
						if (e < MAX_EXP) process(platter[e], Uint32(cur2 >> (e + 1)), n);
						else process(platter[MAX_EXP], cur2 >> MAX_EXP, n);
					}
				}

				cur >>= STRIDE_BITS;
				cur |= word << (32 - STRIDE_BITS);
				//word >>= STRIDE_BITS; // should now be zero, no longer useful
				if (!(cur & skip_mask)) {
					Uint32 n = (i << 1) + 1;
					cur |= next << (32 - STRIDE_BITS);
					e = count_low_zeroes32(cur);
					if (e < 32 - BASE_SIG_BITS) process(platter[e], cur >> (e + 1), n);
					else {
						//Uint64 cur2 = cur | (Uint64(word) << 32) | (Uint64(next) << (64 - STRIDE_BITS * 2)); // word is now zero, ignore it
						Uint64 cur2 = cur | (Uint64(next) << (64 - STRIDE_BITS * 2));
						if (e == 32) e += count_low_zeroes32(Uint32(cur2 >> 32));
						if (e < MAX_EXP) process(platter[e], Uint32(cur2 >> (e + 1)), n);
						else process(platter[MAX_EXP], cur2 >> MAX_EXP, n);
					}
				}
			}
		}
		/*
		//global gaps not set up correctly for other strides, skip_mask not checked, other?
		else if (stride_bits_L2 == 3) {
			enum { STRIDE_BITS = 8 };
			for (long i = start; i < max32; i += 1) {
				Uint32 word = data->as32[i];
				for (unsigned long sub_word = 0; sub_word < 32 / STRIDE_BITS; sub_word++) {
					cur >>= STRIDE_BITS;
					cur |= word << (32 - STRIDE_BITS);
					word >>= STRIDE_BITS;//
					unsigned long e = count_low_zeroes32(cur);
					unsigned long sig;
					if (e < 32 - BASE_SIG_BITS) platter[e].process(cur >> (e + 1));
					else {
						Uint64 cur2 = cur | (Uint64(word) << 32) | (Uint64(data->as32[i + 1]) << (64 - (sub_word + 1) * STRIDE_BITS));//
						if (e == 32) e += count_low_zeroes32(Uint32(cur2 >> 32));
						if (e < MAX_EXP) platter[e].process(Uint32(cur2 >> (e + 1)));
						else platter[MAX_EXP].process(cur2 >> MAX_EXP);
					}
				}
			}
		}
		else if (stride_bits_L2 == 2) {
			enum { STRIDE_BITS = 4 };
			for (long i = start; i < max32; i += 1) {
				Uint32 word = data->as32[i];
				for (unsigned long sub_word = 0; sub_word < 32 / STRIDE_BITS; sub_word++) {
					cur >>= STRIDE_BITS;
					cur |= word << (32 - STRIDE_BITS);
					word >>= STRIDE_BITS;//
					unsigned long e = count_low_zeroes32(cur);
					unsigned long sig;
					if (e < 32 - BASE_SIG_BITS) platter[e].process(cur >> (e + 1));
					else {
						Uint64 cur2 = cur | (Uint64(word) << 32) | (Uint64(data->as32[i + 1]) << (64 - (sub_word + 1) * STRIDE_BITS));//
						if (e == 32) e += count_low_zeroes32(Uint32(cur2 >> 32));
						if (e < MAX_EXP) platter[e].process(Uint32(cur2 >> (e + 1)));
						else platter[MAX_EXP].process(cur2 >> MAX_EXP);
					}
				}
			}
		}
		*/
	}
#endif
	blocks_tested += numblocks;
}

PractRand::Tests::Birthday32::Birthday32() {
	if (BUFFER_SIZE * sizeof(buffer[0]) < TestBlock::SIZE) issue_error("Birthday32 - bad BUFFER_SIZE");
}
void PractRand::Tests::Birthday32::init(PractRand::RNGs::vRNG *known_good) {
	num_buffered = 0;
	for (int i = 0; i < MAX_DUPLICATES; i++) counts[i] = 0;
}
std::string Tests::Birthday32::get_name() const {
	return "BDay32";
}
void PractRand::Tests::Birthday32::get_results(std::vector<TestResult> &results) {
	Uint64 counts2[MAX_DUPLICATES];
	Uint64 total = 0;
	for (int i = 0; i < MAX_DUPLICATES; i++) {
		counts2[i] = counts[i];
		total += counts2[i];
	}
	if (total < 500) return;
	double probs[MAX_DUPLICATES];
	const int lambda_L2 = BUFFER_SIZE_L2 * 3 - 32 - 2;
	if (lambda_L2 < 1) issue_error("Birthday32 - bad configuration");
	const double lambda = 1 << lambda_L2;
	for (int i = 0; i < MAX_DUPLICATES + lambda * 4; i++) {
		double p = std::exp(-lambda) * std::pow(lambda, i) / math_factorial(i);
		if (i < MAX_DUPLICATES) probs[i] = p;
		else probs[MAX_DUPLICATES - 1] += p;
	}
	int cat = simplify_prob_table(MAX_DUPLICATES, total / 40.0, probs, counts2, true, false);
	double chisquared = g_test(cat, probs, counts2);
	double norm = math_chisquared_to_normal(chisquared, cat - 1);
	results.push_back(TestResult(get_name(), norm, norm, TestResult::TYPE_RAW_NORMAL, 0.125));
}
void PractRand::Tests::Birthday32::flush_buffer() {
	if (num_buffered != BUFFER_SIZE) issue_error("Birthday32::flush_buffer - buffer not full?");
	num_buffered = 0;
	enum {
		SORTHELP_SIZE_L2 = BUFFER_SIZE_L2 - 3,
		SORTHELP_SIZE = 1 << SORTHELP_SIZE_L2
	};
	Uint32 sort_helper[SORTHELP_SIZE + 1];
	std::memset(sort_helper, 0, (SORTHELP_SIZE + 1) * sizeof(Uint32));
	for (int i = 0; i < BUFFER_SIZE; i++) sort_helper[1 + (buffer[i] >> (32 - SORTHELP_SIZE_L2))]++;
	int running_total = 0;
	for (int i = 2; i <= SORTHELP_SIZE; i++) {
		sort_helper[i] += sort_helper[i-1];
	}
	Uint32 sorted_buffer[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++) {
		Uint32 value = buffer[i];
		int hi = value >> (32 - SORTHELP_SIZE_L2);
		sorted_buffer[sort_helper[hi]++] = value;
	}
	if (true) {
		int old = 0;
		for (int i = 0; i < SORTHELP_SIZE; i++) {
			int end = sort_helper[i];
			std::sort(&sorted_buffer[old], &sorted_buffer[end]);//is this reasonably optimized for very small sorts?
			old = end;
		}
	}
	Uint32 deltas[BUFFER_SIZE - 1];//possibly some sort of mask could work better instead?
	for (int i = 0; i < BUFFER_SIZE - 1; i++) deltas[i] = sorted_buffer[i + 1] - sorted_buffer[i];
	//tired of trying to write optimized sorting, just doing it the easy way now:
	std::sort(&deltas[0], &deltas[BUFFER_SIZE - 1]);
	int duplicates = 0;
	for (int i = 1; i < BUFFER_SIZE - 1; i++) if (deltas[i] == deltas[i - 1]) duplicates++;
	if (duplicates >= MAX_DUPLICATES) duplicates = MAX_DUPLICATES - 1;
	counts[duplicates]++;
}
void PractRand::Tests::Birthday32::test_blocks(TestBlock *data, int numblocks) {
	while (numblocks) {
		Uint32 *ptr = &buffer[num_buffered];
		if (false) buffer[num_buffered++] = data[0].as32[0];
		else {
			std::memcpy(&buffer[num_buffered], &data[0].as32[0], TestBlock::SIZE);
			num_buffered += TestBlock::SIZE / sizeof(Uint32);
		}		
		if (num_buffered == BUFFER_SIZE) flush_buffer();
		numblocks--;
		data++;
	}
}


void PractRand::Tests::Birthday64::_weird_sort64(Uint64 *base, long length) {
	if (!length) return;
	Uint64 ends[64 + 1];
	for (long i = 0; i <= 64; i++) ends[i] = 0;
	for (long i = 0; i < length; i++) ends[count_ones64(base[i])]++;
	for (long i = 1; i <= 64; i++) ends[i] += ends[i - 1];
	Uint64 starts[64 + 1];
	starts[0] = 0;
	for (long i = 1; i <= 64; i++) starts[i] = ends[i-1];
	long hw = 0;
	while (!ends[hw]) hw++;
	long i = 0;
	while (true) {
		Uint64 v = base[i];
		long bc = count_ones64(v);
		long i2 = starts[bc]++;
		if (i2 > i) {
			base[i] = base[i2];
			base[i2] = v;
			continue;
		}
		else {
			i++;
			while (i >= ends[hw]) {
				if (hw >= 64) break;
				hw++;
				i = starts[hw];
			}
		}
	}
	i = 0;
	for (hw = 0; hw <= 64; hw++) {
		if (i < ends[hw]) std::sort(&base[i], &base[ends[hw]]);
	}
}
void PractRand::Tests::Birthday64::_histogram_in_place_sort64(Uint64 *buffer, long length, long bits_already, Uint32 region_counts[1 << SORT_HELPER_BITS]) {
	if (length <= 1) return;
	if (length < (4 << SORT_HELPER_BITS) || 64 - bits_already - SORT_HELPER_BITS < SORT_HELPER_BITS) {
		std::sort(&buffer[0], &buffer[length]);
		return;
	}
	Uint64 total = 0;
	Uint32 region_bases[(1 << SORT_HELPER_BITS) + 1];
	Uint32 region_bases2[(1 << SORT_HELPER_BITS) + 1];
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	region_bases[0] = 0;
	region_bases2[0] = 0;
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) {
		region_bases[i + 1] = region_bases[i] + region_counts[i];
		region_bases2[i + 1] = region_bases[i + 1];
		total += region_counts[i];
	}
	if (total != length) issue_error("Birthday64::_histogram_in_place_sort64 - bad region counts");

	long region = 0;
	while (region < (1 << SORT_HELPER_BITS)) {
		long index = region_bases[region];
		if (!region_counts[region]) {
			region++;
			continue;
		}
		Uint64 value = buffer[index];
		long ri = (value >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		if (ri == region) {
			region_counts[region]--;
			region_bases[region]++;
			index++;
			continue;
		}

		while (true) {
			long index2 = region_bases[ri];
			Uint64 value2 = buffer[index2];
			long ri2 = (value2 >> shift) & ((1 << SORT_HELPER_BITS) - 1);
			if (ri2 == ri) {
				region_counts[ri]--;
				region_bases[ri]++;
				continue;
			}
			buffer[index2] = value;
			region_counts[ri]--;
			region_bases[ri]++;

			value = value2;
			ri = ri2;
			if (ri == region) break;
		}
		buffer[index] = value;
		region_counts[region]--;
		region_bases[region]++;
	}
	shift -= SORT_HELPER_BITS;
	for (int i = 0; i < 1 << SORT_HELPER_BITS; i++) {
		if (false) {
			std::sort(&buffer[region_bases2[i]], &buffer[region_bases2[i + 1]]);
		}
		else {
			int len = region_bases2[i + 1] - region_bases2[i];
			if (len >= (4 << SORT_HELPER_BITS)) {
				std::memset(region_counts, 0, sizeof(region_counts[0]) << SORT_HELPER_BITS);
				for (int x = region_bases2[i]; x < region_bases2[i + 1]; x++) {
					long ri = (buffer[x] >> shift) & ((1 << SORT_HELPER_BITS) - 1);
					region_counts[ri]++;
				}
				_histogram_in_place_sort64(&buffer[region_bases2[i]], len, bits_already + SORT_HELPER_BITS, region_counts);
			}
			else std::sort(&buffer[region_bases2[i]], &buffer[region_bases2[i + 1]]);
		}
	}
}
void PractRand::Tests::Birthday64::_histogram_in_place_sort64(Uint64 *base, long length) {
	if (length <= 1) return;
	/*if (bits_already >= 64 - SORT_HELPER_BITS) {// - no longer applies since bits_already is assumed to be zero
		std::sort(&base[0], &base[length]);
		return;
	}//*/
	Uint32 region_count[1 << SORT_HELPER_BITS];
	std::memset(region_count, 0, sizeof(region_count[0])* (1 << SORT_HELPER_BITS));
	long shift = 64 - SORT_HELPER_BITS;// -bits_already;
	//Uint64 any_1s = 0, all_1s = 0xFFffFFffFFffFFffull;//for debuging only, disable in release - no longer applies since bits_already is assumed to be zero
	for (Uint64 *cur = base; cur < &base[length]; cur++) {
		long ri = (*cur >> shift);// &((1 << SORT_HELPER_BITS) - 1);
		region_count[ri]++;
		//any_1s |= *cur; all_1s &= *cur;//for debuging only, disable in release - no longer applies since bits_already is assumed to be zero
	}
	//if ((any_1s ^ all_1s) >> bits_already) issue_error("_histogram_sort64 - bits_already not already sorted");//for debuging only, disable in release - no longer applies since bits_already is assumed to be zero
	_histogram_in_place_sort64(base, length, 0, region_count);
}
void PractRand::Tests::Birthday64::_histogram_sort64(Uint64 *buffer, long length, long bits_already, Uint32 region_counts[1 << SORT_HELPER_BITS]) {
	if (length <= 1) return;
	std::vector<Uint64> copied_buffer; copied_buffer.resize(length);
	std::memcpy(&copied_buffer[0], &buffer[0], sizeof(buffer[0]) * length);
	Uint32 region_bases[(1 << SORT_HELPER_BITS) + 1];
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	region_bases[0] = 0;
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) {
		region_bases[i + 1] = region_bases[i] + region_counts[i];
	}
	for (long i = 0; i < length; i++) {
		Uint64 value = copied_buffer[i];
		long ri = (value >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		buffer[region_bases[ri]++] = value;
	}
	copied_buffer.clear();
	shift -= SORT_HELPER_BITS;
	bits_already += SORT_HELPER_BITS;
	long begin = 0;
	for (long ri = 0; ri < (1 << SORT_HELPER_BITS); ri++) {
		long end = region_bases[ri];
		long begin = ri ? region_bases[ri - 1] : 0;
		long len = end - begin;
		if (len > (32 << SORT_HELPER_BITS)) {
			std::memset(region_counts, 0, sizeof(region_counts[0])* (1 << SORT_HELPER_BITS));
			for (long i = begin; i < end; i++) {
				Uint64 value = buffer[i];
				long ri2 = (value >> shift) & ((1 << SORT_HELPER_BITS) - 1);
				region_counts[ri2]++;
			}
			_histogram_sort64(&buffer[begin], end - begin, bits_already, region_counts);
		}
		else std::sort(&buffer[begin], &buffer[end]);
		begin = end;
	}
}
void PractRand::Tests::Birthday64::_histogram_sort64(Uint64 *base, long length) {
	Uint32 region_count[1 << SORT_HELPER_BITS];
	std::memset(region_count, 0, sizeof(region_count[0])* (1 << SORT_HELPER_BITS));
	long shift = 64 - SORT_HELPER_BITS;// -bits_already;
	for (Uint64 *cur = base; cur < &base[length]; cur++) {
		long ri = (*cur >> shift);
		region_count[ri]++;
	}
	_histogram_sort64(base, length, 0, region_count);
}
PractRand::Tests::Birthday64::Birthday64() {
	if (BUFFER_SIZE * sizeof(buffer[0]) < TestBlock::SIZE) issue_error("Birthday64 - bad BUFFER_SIZE");
}
void PractRand::Tests::Birthday64::init(PractRand::RNGs::vRNG *known_good) {
	num_buffered = 0;
	for (int i = 0; i < MAX_DUPLICATES; i++) counts[i] = 0;
}
std::string Tests::Birthday64::get_name() const {
	return "BDay64";
}
void PractRand::Tests::Birthday64::get_results(std::vector<TestResult> &results) {
	Uint64 counts2[MAX_DUPLICATES];
	Uint64 total = 0;
	for (int i = 0; i < MAX_DUPLICATES; i++) {
		counts2[i] = counts[i];
		total += counts2[i];
	}
	if (total < 200) return;
	double probs[MAX_DUPLICATES];
	const double lambda = std::pow(2.0, BUFFER_SIZE_L2 * 3 - 64 - 2);
	for (int i = 0; i < MAX_DUPLICATES + lambda * 4; i++) {
		double p = std::exp(-lambda) * std::pow(lambda, i) / math_factorial(i);
		if (i < MAX_DUPLICATES) probs[i] = p;
		else probs[MAX_DUPLICATES - 1] += p;
	}
	int cat = simplify_prob_table(MAX_DUPLICATES, total / 40.0, probs, counts2, true, false);
	double chisquared = g_test(cat, probs, counts2);
	double norm = math_chisquared_to_normal(chisquared, cat - 1);
	results.push_back(TestResult(get_name(), norm, norm, TestResult::TYPE_RAW_NORMAL, 0.125));
}
void PractRand::Tests::Birthday64::flush_buffer() {
	if (num_buffered != BUFFER_SIZE) issue_error("Birthday64::flush_buffer - buffer not full?");
	num_buffered = 0;
	//std::sort(&buffer[0], &buffer[BUFFER_SIZE]);
	_histogram_in_place_sort64(buffer, BUFFER_SIZE);// fastest in my tests
	//_histogram_sort64(buffer, BUFFER_SIZE);//
	//now it's sorted
	for (int i = 0; i < BUFFER_SIZE - 1; i++) {
		//if (buffer[i] > buffer[i + 1]) issue_error("sort broken");
		buffer[i] = buffer[i + 1] - buffer[i];
	}
	//now it's unsorted deltas of sorted values
	//std::sort(&buffer[0], &buffer[BUFFER_SIZE-1]);
	_histogram_in_place_sort64(buffer, BUFFER_SIZE - 1);//fastest even here, a case it's not optimized for
	//_histogram_sort64(buffer, BUFFER_SIZE - 1);
	//now it's sorted deltas of sorted values
	/*double mean = 0, low = 999999999999999999999999999999.0, high = 0;;
	for (int i = 0; i < BUFFER_SIZE - 1; i++) {
		mean += buffer[i];
		if (buffer[i] < low) low = buffer[i];
		if (buffer[i] > high) high = buffer[i];
	}
	mean /= BUFFER_SIZE - 1;//expected mean is 2**64 / (BUFFER_SIZE + 1), I think
	//std::printf("\nmean=%.6f   low=%.0f   high=%.0f", mean, low, high);
	//std::printf("\nlog2(mean)=%.6f   log2(low)=%.6f   log2(high)=%.6f\n", std::log2(mean), std::log2(low), std::log2(high));
	*/
	int duplicates = 0;
	//std::vector<std::pair<Uint64, Uint64> > repeated_values;
	for (int i = 0; i < BUFFER_SIZE - 2; i++) {
		if (buffer[i] == buffer[i + 1]) {
			//run found
			int first = i;
			while (i < BUFFER_SIZE - 2 && buffer[i] == buffer[i + 1]) i++;
			int run_len = i - first + 1;
			//repeated_values.push_back(std::pair<Uint64, Uint64>(buffer[i], run_len));
			duplicates += run_len - 1;
		}
	}
	/*
	for (int i = 0; i < repeated_values.size(); i++) {
		//if (repeated_values[i].second > 2 || repeated_values[i].first > 1ull << (41 + 2) || repeated_values[i].first < 1ull << (41 - 10)) 
		if (repeated_values[i].second > 2) 
			std::printf("  num = %2d   log2(value) = %.6f   value = %13.0f\n\n", int(repeated_values[i].second), double(std::log2(repeated_values[i].first)), double(repeated_values[i].first));
		//auto it = dup_value_counts.find()
		Uint64 total = dup_value_counts[repeated_values[i].first] += repeated_values[i].second + (1ull << 32);
		if (total > repeated_values[i].second) std::printf("  num = %2d   log2(value) = %.6f   value = %13.0f (combined from %d passes)\n\n", int(total), double(std::log2(repeated_values[i].first)), double(repeated_values[i].first), int(total >> 32));
	}
	*/
	if (duplicates >= MAX_DUPLICATES) duplicates = MAX_DUPLICATES - 1;
	counts[duplicates]++;
}
void PractRand::Tests::Birthday64::test_blocks(TestBlock *data, int numblocks) {
	while (numblocks) {
		Uint64 *ptr = &buffer[num_buffered];
		std::memcpy(&buffer[num_buffered], &data[0].as64[0], TestBlock::SIZE);
		num_buffered += TestBlock::SIZE / sizeof(Uint64);
		if (num_buffered == BUFFER_SIZE) flush_buffer();
		numblocks--;
		data++;
	}
}

void PractRand::Tests::BirthdayHelpers::histogram_sort_and_copy(i128 *buffer, i128 *dest, Uint64 length, long bits_already, Uint64 regions[1 << SORT_HELPER_BITS]) {
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	for (long i = 1; i < (1 << SORT_HELPER_BITS); i++) {
		regions[i] += regions[i - 1];
	}
	for (Uint64 i = 0; i < length; i++) {
		long bin = (buffer[i].high >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		dest[regions[bin]++] = buffer[i];
	}

	for (long i = (1 << SORT_HELPER_BITS) - 1; i > 0; i--) {
		regions[i] = regions[i - 1];
	}
	regions[0] = 0;
	bits_already += SORT_HELPER_BITS;
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) {
		Uint64 sublength = regions[i + 1] - regions[i];
		histogram_in_place_sort128(dest + regions[i], sublength, bits_already);
	}
	for (Uint64 i = 0; i < length; i++) buffer[i] = dest[i];
}
void PractRand::Tests::BirthdayHelpers::histogram_sort_and_copy(i128 *buffer, i128 *dest, Uint64 length, long bits_already) {
	Uint64 regions[(1 << SORT_HELPER_BITS) + 1];
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) regions[i] = 0;
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	for (Uint64 i = 0; i < length; i++) {
		long bin = (buffer[i].high >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		regions[1 + bin]++;
	}
	histogram_sort_and_copy(buffer, dest, length, bits_already, regions);
}
void PractRand::Tests::BirthdayHelpers::radix_sort_and_copy(i128 *buffer, i128 *dest, Uint64 length, long bits_already) {
	if (length <= 1) return;
	enum {COMBINED_PASSES = 3};
	Uint64 regions[COMBINED_PASSES << SORT_HELPER_BITS];
	if (true) {//count frequencies for each pass
		for (long i = 0; i < (COMBINED_PASSES << SORT_HELPER_BITS); i++) regions[i] = 0;
		long shift = 64 - SORT_HELPER_BITS - bits_already;
		Uint64 already_check_mask = ((1ull << bits_already) - 1) << (64 - bits_already);//debuging check, remove sometime
		Uint64 already_check_value = buffer[0].high & already_check_mask;//debuging check, remove sometime
		for (Uint64 i = 0; i < length; i++) {
			Uint64 value = buffer[i].high;
			if ((value & already_check_mask) != already_check_value) issue_error();//debuging check, remove sometime
			value <<= bits_already;

			for (long region_base = 0; region_base < (COMBINED_PASSES << SORT_HELPER_BITS); region_base += (1 << SORT_HELPER_BITS)) {
				long bin = value >> (64 - SORT_HELPER_BITS);
				regions[region_base + bin]++;

				value <<= SORT_HELPER_BITS;
			}
		}
	}
	if (true) {//convert from frequency counts to destination regions
		for (long region_base = 0; region_base < (COMBINED_PASSES << SORT_HELPER_BITS); region_base += (1 << SORT_HELPER_BITS)) {
			Uint64 sum = 0;
			for (long i = 1; i < (1 << SORT_HELPER_BITS); i++) {
				Uint64 new_sum = sum + regions[region_base + i];
				regions[region_base + i] = sum;
				sum = new_sum;
			}
		}
	}
	long shift = 64 - COMBINED_PASSES * SORT_HELPER_BITS - bits_already;
	for (long dimension = COMBINED_PASSES - 1; dimension >= 0; dimension--) {
		// now the actual sorting passes, least significant first
		Uint64 *region_base = &regions[dimension << SORT_HELPER_BITS];
		for (Uint64 i = 0; i < length; i++) {
			long bin = (buffer[i].high >> shift) & ((1 << SORT_HELPER_BITS) - 1);
			dest[region_base[bin]++] = buffer[i];
		}

		shift += SORT_HELPER_BITS;
		i128 *tmp = buffer; buffer = dest; dest = tmp;
	}
	// now the copy and recursive sorting in case that wasn't enough:
	Uint64 sorted_mask = Uint64(0) - Uint64((1ull << (64 - SORT_HELPER_BITS * COMBINED_PASSES - bits_already)) - 1);
	Uint64 run_value = buffer[0].high & sorted_mask;
	Uint64 run_length = 0;
	for (long i = 1; i < (1 << SORT_HELPER_BITS); i++) {
		Uint64 value = buffer[i].high & sorted_mask;
		if (value == run_value) run_length++;//run continuing
		else {
			run_value = value;
			if (!run_length) {// no run, needs a single-copy
				dest[i - 1] = buffer[i - 1];
			}
			else {//run ends, needs internal sorting and a multiple-copy
				Uint64 start_index = i - run_length - 1;
				std::sort(&buffer[start_index], &buffer[i]);
				std::copy(&buffer[start_index], &buffer[i], &dest[start_index]);
				run_length = 0;
			}
		}
	}
}
void PractRand::Tests::BirthdayHelpers::histogram_in_place_sort128(i128 *base, Uint64 length, long bits_already) {
	Uint64 freq_counts[1 << SORT_HELPER_BITS];
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) freq_counts[i] = 0;
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	for (Uint64 i = 0; i < length; i++) {
		freq_counts[(base[i].high >> shift) & ((1ull << SORT_HELPER_BITS) - 1)]++;
	}
	histogram_in_place_sort128(base, length, bits_already, freq_counts);
}
void PractRand::Tests::BirthdayHelpers::histogram_in_place_sort128(i128 *buffer, Uint64 length, long bits_already, Uint64 region_counts[1 << SORT_HELPER_BITS]) {
	if (length <= 1) return;
	if (length < (4 << SORT_HELPER_BITS) || bits_already >= 64 - SORT_HELPER_BITS) {
		std::sort(&buffer[0], &buffer[length]);
		return;
	}
	Uint64 total = 0;
	Uint32 region_bases[(1 << SORT_HELPER_BITS) + 1];
	Uint32 region_bases2[(1 << SORT_HELPER_BITS) + 1];
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	region_bases[0] = 0;
	region_bases2[0] = 0;
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) {
		region_bases[i + 1] = region_bases[i] + region_counts[i];
		region_bases2[i + 1] = region_bases[i + 1];
		total += region_counts[i];
	}
	if (total != length) issue_error("BirthdayHelpers::histogram_in_place_sort128 - bad region counts");

	long region = 0;
	while (region < (1 << SORT_HELPER_BITS)) {
		long index = region_bases[region];
		if (!region_counts[region]) {
			region++;
			continue;
		}
		i128 value = buffer[index];
		long ri = (value.high >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		if (ri == region) {
			region_counts[region]--;
			region_bases[region]++;
			index++;
			continue;
		}

		while (true) {
			if (!region_counts[ri]) {
				issue_error("region already empty");
			}
			long index2 = region_bases[ri];
			i128 value2 = buffer[index2];
			long ri2 = (value2.high >> shift) & ((1 << SORT_HELPER_BITS) - 1);
			if (ri2 == ri) {
				region_counts[ri]--;
				region_bases[ri]++;
				continue;
			}
			buffer[index2] = value;
			region_counts[ri]--;
			region_bases[ri]++;

			value = value2;
			ri = ri2;
			if (ri == region) break;
		}
		buffer[index] = value;
		region_counts[region]--;
		region_bases[region]++;
	}

	shift -= SORT_HELPER_BITS;
	bits_already += SORT_HELPER_BITS;
	for (int i = 0; i < 1 << SORT_HELPER_BITS; i++) {
		if (false || shift >= 64) {
			std::sort(&buffer[region_bases2[i]], &buffer[region_bases2[i + 1]]);
		}
		else {
			int len = region_bases2[i + 1] - region_bases2[i];
			if (len >= (4 << SORT_HELPER_BITS)) {
				std::memset(region_counts, 0, sizeof(region_counts[0]) << SORT_HELPER_BITS);
				for (int x = region_bases2[i]; x < region_bases2[i + 1]; x++) {
					long ri = (buffer[x].high >> shift) & ((1 << SORT_HELPER_BITS) - 1);
					region_counts[ri]++;
				}
				histogram_in_place_sort128(&buffer[region_bases2[i]], len, bits_already, region_counts);
			}
			else std::sort(&buffer[region_bases2[i]], &buffer[region_bases2[i + 1]]);
		}
	}
}
void PractRand::Tests::BirthdayHelpers::_sorted_deltas_of_sorted_values(i128 *base, long length_L2, Uint64 freq_counts[1 << SORT_HELPER_BITS]) {
	if (length_L2 < 1) issue_error();
	long length = 1 << length_L2;
	if (false) {// both sortings use regular algorithms
		std::sort(base, base + length);
		// it's now sorted
		for (int i = 0; i < length - 1; i++) base[i] = base[i + 1] - base[i];
		// it's now deltas of sorted values
		std::sort(base, base + (length - 1));
		// it's now sorted deltas of sorted values
	}
	else if (false) {// 2nd sorting uses regular algorithm, 1st uses a histogram sort
		histogram_in_place_sort128(base, length, 0, freq_counts);
		for (int i = 0; i < length - 1; i++) {
			if (base[i + 1] < base[i]) issue_error("BirthdayHelpers::_sorted_deltas_of_sorted_values - sort1 failed");//debugging only, remove
		}
		for (int i = 0; i < length - 1; i++) base[i] = base[i + 1] - base[i];
		// it's now deltas of sorted values
		std::sort(base, base + (length - 1));
		// it's now sorted deltas of sorted values
	}
	else if (true) {// both sortings use variants of histogram sort
		histogram_in_place_sort128(base, length, 0, freq_counts);
		//for (int i = 0; i < length - 1; i++) { if (base[i + 1] < base[i]) issue_error("BirthdayHelpers::_sorted_deltas_of_sorted_values - sort1 failed"); }//debugging only, remove
		// it's now sorted
		enum { SAFETY_MARGIN = 2 };
		std::memset(freq_counts, 0, sizeof(freq_counts[0]) << SORT_HELPER_BITS);
		std::vector<i128> spills;
		long shift = 64 - SORT_HELPER_BITS - length_L2 + SAFETY_MARGIN;
		for (Uint64 i = 0; i < length - 1; i++) {
			i128 delta = base[i + 1] - base[i];
			base[i] = delta;
			long ri = delta.high >> shift;
			if (ri >= (1 << SORT_HELPER_BITS)) spills.push_back(delta);
			else {
				freq_counts[ri]++;
				base[i - spills.size()] = delta;
			}
		}
		// it's now deltas of sorted values
		histogram_in_place_sort128(base, length - 1 - spills.size(), length_L2 - SAFETY_MARGIN, freq_counts);
		std::sort(spills.begin(), spills.end());
		for (int i = 0; i < spills.size(); i++) base[length - 1 - spills.size() + i] = spills[i];
		//for (int i = 0; i < length - 2; i++) {if (base[i + 1] < base[i]) issue_error("BirthdayHelpers::_sorted_deltas_of_sorted_values - sort2 failed");}//debugging only, remove
		// it's now sorted deltas of sorted values
	}
	else if (false) {// untested code, using bitmask for lower bits of 2nd sorting

	}
	else if (false) {// ...trying a radix-sort, despite not needing the extra copy of the data
		i128 *buffer2 = new i128[length];
		//histogram_in_place_sort128(base, length, 0, freq_counts);
		radix_sort_and_copy(base, buffer2, length, 0);
		for (int i = 0; i < length - 1; i++) {
			if (base[i + 1] < base[i]) issue_error("BirthdayHelpers::_sorted_deltas_of_sorted_values - sort1 failed");//debugging only, remove
		}
		// it's now sorted
		enum { SAFETY_MARGIN = 2 };
		std::memset(freq_counts, 0, sizeof(freq_counts[0]) << SORT_HELPER_BITS);
		std::vector<i128> spills;
		long shift = 64 - SORT_HELPER_BITS - length_L2 + SAFETY_MARGIN;
		for (Uint64 i = 0; i < length - 1; i++) {
			i128 delta = base[i + 1] - base[i];
			base[i] = delta;
			long ri = delta.high >> shift;
			if (ri >= (1 << SORT_HELPER_BITS)) spills.push_back(delta);
			else {
				freq_counts[ri]++;
				base[i - spills.size()] = delta;
			}
		}
		// it's now deltas of sorted values
		histogram_in_place_sort128(base, length - 1 - spills.size(), length_L2 - SAFETY_MARGIN, freq_counts);
		std::sort(spills.begin(), spills.end());
		for (int i = 0; i < spills.size(); i++) base[length - 1 - spills.size() + i] = spills[i];
		for (int i = 0; i < length - 2; i++) {
			if (base[i + 1] < base[i]) issue_error("BirthdayHelpers::_sorted_deltas_of_sorted_values - sort2 failed");//debugging only, remove
		}
		// it's now sorted deltas of sorted values
		delete[] buffer2;
	}
	else issue_error();
}
void PractRand::Tests::BirthdayHelpers::_sorted_deltas_of_sorted_values(i128 *base, long length_L2) {
	if (length_L2 < 1) issue_error();
	Uint64 length = 1 << length_L2;
	Uint64 freq_counts[1 << SORT_HELPER_BITS];
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) freq_counts[i] = 0;
	for (Uint64 i = 0; i < length; i++) {
		freq_counts[base[i].high >> (64 - SORT_HELPER_BITS)]++;
	}
	_sorted_deltas_of_sorted_values(base, length_L2, freq_counts);
}
void PractRand::Tests::BirthdayHelpers::histogram_in_place_sort64(Uint64 *buffer, Uint64 length, long bits_already, Uint64 region_counts[1 << SORT_HELPER_BITS]) {
	if (length <= 1) return;
	if (length < (4 << SORT_HELPER_BITS) || 64 - bits_already - SORT_HELPER_BITS < SORT_HELPER_BITS) {
		std::sort(&buffer[0], &buffer[length]);
		return;
	}
	Uint64 total = 0;
	Uint32 region_bases[(1 << SORT_HELPER_BITS) + 1];
	Uint32 region_bases2[(1 << SORT_HELPER_BITS) + 1];
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	region_bases[0] = 0;
	region_bases2[0] = 0;
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) {
		region_bases[i + 1] = region_bases[i] + region_counts[i];
		region_bases2[i + 1] = region_bases[i + 1];
		total += region_counts[i];
	}
	if (total != length) issue_error("BirthdayHelpers::histogram_in_place_sort64 - bad region counts");

	long region = 0;
	while (true) {
		long index = region_bases[region];
		if (!region_counts[region]) {
			region++;
			if (region >= (1 << SORT_HELPER_BITS)) break;
			continue;
		}
		Uint64 value = buffer[index];
		long ri = (value >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		if (ri == region) {
			region_counts[region]--;
			region_bases[region]++;
			continue;
		}

		while (true) {
			long index2 = region_bases[ri];
			Uint64 value2 = buffer[index2];
			long ri2 = (value2 >> shift) & ((1 << SORT_HELPER_BITS) - 1);
			if (ri2 == ri) {
				region_counts[ri]--;
				region_bases[ri]++;
				continue;
			}
			buffer[index2] = value;
			region_counts[ri]--;
			region_bases[ri]++;

			value = value2;
			ri = ri2;
			if (ri == region) break;
		}
		buffer[index] = value;
		region_counts[region]--;
		region_bases[region]++;
	}
	shift -= SORT_HELPER_BITS;
	for (int i = 0; i < 1 << SORT_HELPER_BITS; i++) {
		if (false) {
			std::sort(&buffer[region_bases2[i]], &buffer[region_bases2[i + 1]]);
		}
		else {
			int len = region_bases2[i + 1] - region_bases2[i];
			if (len >= (4 << SORT_HELPER_BITS)) {
				std::memset(region_counts, 0, sizeof(region_counts[0]) << SORT_HELPER_BITS);
				for (int x = region_bases2[i]; x < region_bases2[i + 1]; x++) {
					long ri = (buffer[x] >> shift) & ((1 << SORT_HELPER_BITS) - 1);
					region_counts[ri]++;
				}
				histogram_in_place_sort64(&buffer[region_bases2[i]], len, bits_already + SORT_HELPER_BITS, region_counts);
			}
			else std::sort(&buffer[region_bases2[i]], &buffer[region_bases2[i + 1]]);
		}
	}
}
void PractRand::Tests::BirthdayHelpers::histogram_in_place_sort64(Uint64 *buffer, Uint64 length, long bits_already) {
	if (length <= 1) return;
	Uint64 region_count[1 << SORT_HELPER_BITS];
	std::memset(region_count, 0, sizeof(region_count[0]) * (1 << SORT_HELPER_BITS));
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	for (Uint64 *cur = buffer, *end = &buffer[length]; cur < end; cur++) {
		long ri = (*cur >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		region_count[ri]++;
	}
	histogram_in_place_sort64(buffer, length, bits_already, region_count);
}
void PractRand::Tests::BirthdayHelpers::histogram_sort_and_copy64(Uint64 *buffer, Uint64 *dest, Uint64 length, long bits_already, Uint64 freq_counts[1 << SORT_HELPER_BITS]) {
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	for (long i = 1; i < (1 << SORT_HELPER_BITS); i++) {
		freq_counts[i] += freq_counts[i - 1];
	}
	for (Uint64 i = 0; i < length; i++) {
		long bin = (buffer[i] >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		dest[freq_counts[bin]++] = buffer[i];
	}

	for (long i = (1 << SORT_HELPER_BITS) - 1; i > 0; i--) {
		freq_counts[i] = freq_counts[i - 1];
	}
	freq_counts[0] = 0;
	bits_already += SORT_HELPER_BITS;
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) {
		Uint64 sublength = freq_counts[i + 1] - freq_counts[i];
		histogram_in_place_sort64(dest + freq_counts[i], sublength, bits_already);
	}
	for (Uint64 i = 0; i < length; i++) buffer[i] = dest[i];
}
void PractRand::Tests::BirthdayHelpers::histogram_sort_and_copy64(Uint64 *buffer, Uint64 *dest, Uint64 length, long bits_already) {
	Uint64 regions[(1 << SORT_HELPER_BITS) + 1];
	for (long i = 0; i < (1 << SORT_HELPER_BITS); i++) regions[i] = 0;
	long shift = 64 - SORT_HELPER_BITS - bits_already;
	for (Uint64 i = 0; i < length; i++) {
		long bin = (buffer[i] >> shift) & ((1 << SORT_HELPER_BITS) - 1);
		regions[1 + bin]++;
	}
	histogram_sort_and_copy64(buffer, dest, length, bits_already, regions);
}
void PractRand::Tests::BirthdayHelpers::radix_sort_and_copy64(Uint64 *buffer, Uint64 *dest, Uint64 length, long bits_already) {
	if (length <= 1) return;
	enum { COMBINED_PASSES = 3 };
	Uint64 regions[COMBINED_PASSES << SORT_HELPER_BITS];
	if (true) {//count frequencies for each pass
		for (long i = 0; i < (COMBINED_PASSES << SORT_HELPER_BITS); i++) regions[i] = 0;
		long shift = 64 - SORT_HELPER_BITS - bits_already;
		Uint64 already_check_mask = ((1ull << bits_already) - 1) << (64 - bits_already);//debuging check, remove sometime
		Uint64 already_check_value = buffer[0] & already_check_mask;//debuging check, remove sometime
		for (Uint64 i = 0; i < length; i++) {
			Uint64 value = buffer[i];
			if ((value & already_check_mask) != already_check_value) issue_error();//debuging check, remove sometime
			value <<= bits_already;

			for (long region_base = 0; region_base < (COMBINED_PASSES << SORT_HELPER_BITS); region_base += (1 << SORT_HELPER_BITS)) {
				long bin = value >> (64 - SORT_HELPER_BITS);
				regions[region_base + bin]++;

				value <<= SORT_HELPER_BITS;
			}
		}
	}
	if (true) {//convert from frequency counts to destination regions
		for (long region_base = 0; region_base < (COMBINED_PASSES << SORT_HELPER_BITS); region_base += (1 << SORT_HELPER_BITS)) {
			Uint64 sum = 0;
			for (long i = 1; i < (1 << SORT_HELPER_BITS); i++) {
				Uint64 new_sum = sum + regions[region_base + i];
				regions[region_base + i] = sum;
				sum = new_sum;
			}
		}
	}
	long shift = 64 - COMBINED_PASSES * SORT_HELPER_BITS - bits_already;
	for (long dimension = COMBINED_PASSES - 1; dimension >= 0; dimension--) {
		// now the actual sorting passes, least significant first
		Uint64 *region_base = &regions[dimension << SORT_HELPER_BITS];
		for (Uint64 i = 0; i < length; i++) {
			long bin = (buffer[i] >> shift) & ((1 << SORT_HELPER_BITS) - 1);
			dest[region_base[bin]++] = buffer[i];
		}

		shift += SORT_HELPER_BITS;
		Uint64 *tmp = buffer; buffer = dest; dest = tmp;
	}
	// now the copy and recursive sorting in case that wasn't enough:
	Uint64 sorted_mask = Uint64(0) - Uint64((1ull << (64 - SORT_HELPER_BITS * COMBINED_PASSES - bits_already)) - 1);
	Uint64 run_value = buffer[0] & sorted_mask;
	Uint64 run_length = 0;
	for (long i = 1; i < (1 << SORT_HELPER_BITS); i++) {
		Uint64 value = buffer[i] & sorted_mask;
		if (value == run_value) run_length++;//run continuing
		else {
			run_value = value;
			if (!run_length) {// no run, needs a single-copy
				dest[i - 1] = buffer[i - 1];
			}
			else {//run ends, needs internal sorting and a multiple-copy
				Uint64 start_index = i - run_length - 1;
				std::sort(&buffer[start_index], &buffer[i]);
				std::copy(&buffer[start_index], &buffer[i], &dest[start_index]);
				run_length = 0;
			}
		}
	}
}
bool PractRand::Tests::BirthdayHelpers::count_low16_duplicates(Uint64 *buffer, Uint64 length, Uint64 &_count) {
	enum {
		BITS = 16,
		WORD_BITS_L2 = 6,
		WORD_BITS = 1 << WORD_BITS_L2,
		ARRAY_BITS = 1ull << BITS,
	};
	Uint64 bitmask[ARRAY_BITS / WORD_BITS];// that's 8 kilobytes for 16 bits
	std::memset(bitmask, 0, ARRAY_BITS / 8);
	long count = 0;
	Uint64 ored;
	Uint64 *end = buffer + length;
	for (; buffer != end; buffer++) {Uint64 current = *buffer;
	//for (Uint64 i = 0; i < length; i++) {Uint64 current = buffer[i];
		ored |= current;
		long low16 = current  & (ARRAY_BITS - 1);
		long shift = low16 & 63;
		long index = low16 >> WORD_BITS_L2;
		Uint64 shifted = 1ull << shift;
		Uint64 &indexed = bitmask[index];
		if (indexed & shifted) count++;
		indexed |= shifted;
	}
	_count = count;
	return (ored >> BITS) == 0;
}



PractRand::Tests::BirthdayLambda1::~BirthdayLambda1() {
	delete[] buffer;
}
PractRand::Tests::BirthdayLambda1::BirthdayLambda1(int buffer_size_L2_) : buffer_size_L2(buffer_size_L2_) {
	bits_to_use = buffer_size_L2 * 3 - 2;
	//const double lambda = std::pow(2.0, buffer_size_L2 * 3 - BITS_TO_USE - 2);
	//if (lambda != 1.0) issue_error("BirthdayLambda1 - bad configuration");
	duplicates = expected_duplicates = 0;
	if (bits_to_use > 128 || bits_to_use < 1) issue_error("BirthdayLamda1 - bad bits_to_use");
	if ((1ull << buffer_size_L2) * sizeof(buffer[0]) < TestBlock::SIZE) issue_error("BirthdayLamda1 - bad buffer size");

	buffer = new i128[1 << buffer_size_L2];
	//buffer.resize(1 << buffer_size_L2);
	if (!buffer) issue_error("BirthdayLamda1 - failed to allocate buffer");
}
void PractRand::Tests::BirthdayLambda1::init(PractRand::RNGs::vRNG *known_good) {
	num_buffered = 0;
	for (int i = 0; i < (1 << SORT_HELPER_BITS); i++) sort_helper_counts[i] = 0;
	autofail = false;
}
std::string Tests::BirthdayLambda1::get_name() const {
	std::stringstream buf;
	buf << "BDayL1(" << buffer_size_L2 << ")";
	return buf.str();
}
static double largest_spacing_cdf(Uint64 N, double value) {
	long double invN = 1.0 / N;
	if (value < invN) return 0;
	if (value > 1) return 1;
	long double invX = 1.0 / value;
	Uint64 mink, maxk;
	bool invert_result;
	if (invX > N / 2) { mink = int(invX + 1); maxk = N; invert_result = true; }
	else { mink = 1; maxk = int(invX); invert_result = false; }
	long double p = 0, pp = 0;
	for (Uint64 k = mink; k <= maxk; k++) {
		long double term = (((k - 1) & 1) ? -1 : 0);
		term *= std::pow(1 - k * value, N - 1);
		if (k > N) issue_error();//multiply by zero otherwise, but I don't think that should happen?
		using PractRand::Tests::math_factorial_log;
		long double choose_n_of_k = math_factorial_log(N) - math_factorial_log(k) - math_factorial_log(N - k);
		choose_n_of_k = std::exp(choose_n_of_k);
		term *= choose_n_of_k;
		if (pp) {
			p += (pp + term);
			pp = 0;
		}
		else pp = term;
	}
	p += pp;
	if (invert_result) p = 1 - p;
	return p;
}
void PractRand::Tests::BirthdayLambda1::get_results(std::vector<TestResult> &results) {
	if (autofail) {
		results.push_back(TestResult(get_name() + ":!", -1, 1, TestResult::TYPE_PASSFAIL, 0.125));
		return;
	}
	if (!expected_duplicates) return;


	double norm = (duplicates - expected_duplicates) / std::sqrt(expected_duplicates);

	//results.push_back(TestResult(get_name() + ":gamma", norm, math_low))

	if (expected_duplicates < 100) {
		if (norm <= 0) {
			double p = 0;
			for (int i = 0; i < duplicates; i++) p += math_poisson_pmf(expected_duplicates, i);
			p += 0.5 * math_poisson_pmf(expected_duplicates, duplicates);
			std::vector<double> probs; probs.resize(duplicates + 1);
			for (int i = 0; i <= duplicates; i++) probs[i] = math_poisson_pmf(expected_duplicates, i);
			results.push_back(TestResult(get_name(), norm, 1-p, TestResult::TYPE_BAD_P, 0.125));
		}
		else {
			double p = 0;
			double high = duplicates > expected_duplicates ? duplicates : expected_duplicates;
			high += 8 * std::sqrt(expected_duplicates);
			for (int i = high + 3; i > duplicates; i--) p += math_poisson_pmf(expected_duplicates, i);
			p += 0.5 * math_poisson_pmf(expected_duplicates, duplicates);
			results.push_back(TestResult(get_name(), norm, p, TestResult::TYPE_BAD_P, 0.125));
		}
	}
	else {
		results.push_back(TestResult(get_name(), norm, math_normaldist_to_pvalue(-norm), TestResult::TYPE_BAD_P, 0.125));
		// would be better passed through calibration, but this will do for now
	}
}
Uint64 PractRand::Tests::BirthdayLambda1::flush_buffer() {
	const Uint64 buffer_size = 1ull << buffer_size_L2;
	if (num_buffered != buffer_size) issue_error("BirthdayLamda1::flush_buffer - buffer not full?");
	num_buffered = 0;
	if (autofail) return 0;

	//BirthdayHelpers::_sorted_deltas_of_sorted_values(&buffer[0], buffer_size_L2, sort_helper_counts);
	BirthdayHelpers::histogram_in_place_sort128(buffer, buffer_size, 0, sort_helper_counts);
	// it's now sorted
	enum { SAFETY_MARGIN = 2 };
	std::memset(sort_helper_counts, 0, sizeof(sort_helper_counts[0]) << SORT_HELPER_BITS);
	std::vector<i128> spills;
	i128 largest; largest.high = largest.low = 0;
	long shift = 64 - SORT_HELPER_BITS - buffer_size_L2 + SAFETY_MARGIN;
	for (Uint64 i = 0; i < buffer_size - 1; i++) {
		i128 delta = buffer[i + 1] - buffer[i];
		buffer[i] = delta;
		if (DO_LARGEST_SPACING && largest < delta) {
			largest = delta;
		}
		long ri = delta.high >> shift;
		if (ri >= (1 << SORT_HELPER_BITS)) spills.push_back(delta);
		else {
			sort_helper_counts[ri]++;
			buffer[i - spills.size()] = delta;
		}
	}
	// it's now deltas of sorted values
	BirthdayHelpers::histogram_in_place_sort128(buffer, buffer_size - 1 - spills.size(), buffer_size_L2 - SAFETY_MARGIN, sort_helper_counts);
	std::sort(spills.begin(), spills.end());
	for (int i = 0; i < spills.size(); i++) buffer[buffer_size - 1 - spills.size() + i] = spills[i];
	//for (int i = 0; i < length - 2; i++) {if (base[i + 1] < base[i]) issue_error("BirthdayHelpers::_sorted_deltas_of_sorted_values - sort2 failed");}//debugging only, remove
	// it's now sorted deltas of sorted values
	std::memset(sort_helper_counts, 0, sizeof(sort_helper_counts[0]) << BirthdayHelpers::SORT_HELPER_BITS); // have to reset the histograms even if we never use them, to prevent overflow if nothing else

	if (DO_LARGEST_SPACING) {
		double d_largest = (largest.high + (largest.low / 18446744073709551616.0)) / 18446744073709551616.0;
		double largest_spacing_p = largest_spacing_cdf(buffer_size, d_largest);
		if (largest_spacing_p > longest_spacing) longest_spacing = largest_spacing_p;
	}

	//std::vector<std::pair<i128, Uint64> > repeated_values;
	Uint64 rv = 0;
	for (int i = 0; i < buffer_size - 2; i++) {
		if (buffer[i] == buffer[i + 1]) {
			//run found
			int first = i;
			while (i < buffer_size - 2 && buffer[i] == buffer[i + 1]) i++;
			int run_len = i - first + 1;
			//repeated_values.push_back(std::pair<Uint64, Uint64>(sorted_buffer[i], run_len));
			rv += run_len - 1;
		}
	}
	duplicates += rv;
	expected_duplicates += std::pow(2.0, buffer_size_L2 * 3 - (bits_to_use + 2));
	return rv;
}
void PractRand::Tests::BirthdayLambda1::test_blocks(TestBlock *data, int numblocks) {
	if (autofail) return;

	Uint64 mask_high = Uint64(Sint64(-1)), mask_low;
	if (bits_to_use < 64) {
		mask_low = 0;
		mask_high <<= (64 - bits_to_use);
	}
	else if (bits_to_use == 64) mask_low = 0;
	else if (bits_to_use < 128) {
		mask_low = mask_high << (128 - bits_to_use);
	}
	else if (bits_to_use == 128) mask_low = mask_high;
	else issue_error();

	while (numblocks) {
		i128 *dest = &buffer[num_buffered];
		Uint64 *cur = &data[0].as64[0];
		Uint64 *end = &data[1].as64[0];
		enum { LOW = 0, HIGH = 1 };//that's kind of endian-ist, but the tests generally don't bother dealing with such issues
		for (; cur != end; cur += 2, dest++) {
			dest->low = cur[LOW] & mask_low;
			dest->high = cur[HIGH] & mask_high;
			int ri = dest->high >> (64 - SORT_HELPER_BITS);
			sort_helper_counts[ri]++;
		}
		num_buffered += TestBlock::SIZE / sizeof(i128);
		if (num_buffered >> buffer_size_L2) {
			flush_buffer();
		}
		numblocks--;
		data++;
	}
}


PractRand::Tests::BirthdaySystematic128::BirthdaySystematic128(int bufsize_L2_) : BirthdayLambda1(bufsize_L2_) {}
void PractRand::Tests::BirthdaySystematic128::init(PractRand::RNGs::vRNG *known_good) {
	BirthdayLambda1::init(known_good);

	already_sorted = 0;
	score = 0;
	incomplete_expected_duplicates = 0;
	//incomplete_duplicates = 0;
}
std::string Tests::BirthdaySystematic128::get_name() const {
	std::ostringstream buf;
	buf << "BDayS128(" << buffer_size_L2 << ")";
	return buf.str();
}
void PractRand::Tests::BirthdaySystematic128::do_incomplete_buffer() {
	const Uint64 buffer_size = 1ull << buffer_size_L2;
	const Uint64 half_buffer_size = buffer_size >> 1;
	if (expected_duplicates) issue_error("BirthdaySystematic128 - do_incomplete_buffer should not be called after a full sample");
	if (num_buffered > half_buffer_size) issue_error("BirthdaySystematic128 - do_incomplete_buffer should not be called with a buffer this full");

	already_sorted = 0;
	Uint64 num_unsorted = num_buffered - already_sorted;
	if (num_unsorted) {
		if (already_sorted) {
			BirthdayHelpers::histogram_in_place_sort128(&buffer[already_sorted], num_unsorted);
			Uint64 p1 = 0, p2 = already_sorted, p3 = half_buffer_size, max = half_buffer_size + num_buffered;
			while (true) {
				if (p1 == already_sorted) {
					while (p2 < num_buffered) buffer[p3++] = buffer[p2++];
					break;
				}
				else if (p2 == num_buffered) {
					while (p1 < already_sorted) buffer[p3++] = buffer[p1++];
					break;
				}
				else {
					while (true) {
						if (buffer[p1] < buffer[p2]) {
							buffer[p3++] = buffer[p1++];
							if (p1 == already_sorted) break;
						}
						else {
							buffer[p3++] = buffer[p2++];
							if (p2 == num_buffered) break;
						}
					}
				}
			}
			std::copy(&buffer[half_buffer_size], &buffer[half_buffer_size + num_buffered], &buffer[0]);
		}
		else {
			Uint64 region_counts[1 << BirthdayHelpers::SORT_HELPER_BITS];
			std::copy(&sort_helper_counts[0], &sort_helper_counts[1 << SORT_HELPER_BITS], &region_counts[0]);
			BirthdayHelpers::histogram_in_place_sort128(buffer, num_buffered, 0, region_counts);
			std::copy(&buffer[0], &buffer[num_buffered], &buffer[half_buffer_size]);
		}
		already_sorted = num_buffered;
	}
	else {
		std::copy(&buffer[0], &buffer[num_buffered], &buffer[half_buffer_size]);
	}
	//now we should all be sorted, plus the second half of the buffer should contain a redundant copy of the data - that's why we can't do this if more than half the buffer is full already
//	Uint64 effective_num_buffered = num_buffered;
//	int effective_bufsize_L2 = 10;
	if (num_buffered < 64) return;
	//while ((2ull << effective_bufsize_L2) < num_buffered) effective_bufsize_L2++;
	double log2_of_buffer_size = std::log(double(num_buffered)) / std::log(2.0);
	long bits_per_sample = std::floor(3 * log2_of_buffer_size - 2);
	if (bits_per_sample > bits_to_use) issue_error();
	//const Uint64 effective_buffer_size = 1ull << effective_bufsize_L2;
	Uint64 high_mask = 0xFFffFFffFFffFFffull, low_mask;
	if (bits_per_sample == 128) low_mask = high_mask;
	else if (bits_per_sample > 64) low_mask = high_mask << (128 - bits_per_sample);
	else if (bits_per_sample == 64) low_mask = 0;
	else if (bits_per_sample < 64) {
		low_mask = 0;
		high_mask <<= (64 - bits_per_sample);
	}
	buffer[half_buffer_size].low &= low_mask;
	buffer[half_buffer_size].high &= high_mask;
	for (Uint64 i = 1; i < num_buffered; i++) {
		//buffer[half_buffer_size + i].low = buffer[i].low & low_mask;
		//buffer[half_buffer_size + i].high = buffer[i].high & high_mask;
		buffer[half_buffer_size + i].low &= low_mask;
		buffer[half_buffer_size + i].high &= high_mask;
		buffer[half_buffer_size + i - 1] = buffer[half_buffer_size + i] - buffer[half_buffer_size + i - 1];
	}
	BirthdayHelpers::histogram_in_place_sort128(&buffer[half_buffer_size], num_buffered - 1);
	Uint64 dup = 0;
	for (Uint64 i = 1; i < num_buffered - 1; i++) {
		if (buffer[half_buffer_size + i] == buffer[half_buffer_size + i - 1]) dup++;
	}
	incomplete_duplicates = dup;
	incomplete_expected_duplicates = std::pow(2.0, 3 * log2_of_buffer_size - 2 - bits_per_sample);
	score = evaluate_score(incomplete_expected_duplicates, dup);
	return;
}
void PractRand::Tests::BirthdaySystematic128::get_results(std::vector<TestResult> &results) {
	if (autofail) {
		results.push_back(TestResult(get_name() + ":!", 0, 1.0, TestResult::TYPE_PASSFAIL, 0.01));
		return;
	}
	double total_expected_duplicates = BirthdayLambda1::expected_duplicates;
	double total_actual_duplicates = BirthdayLambda1::duplicates;
	if (!total_expected_duplicates) {// we have not yet completed a full pass
		if (!num_buffered) return;;
		const Uint64 half_buffer_size = 1ull << (buffer_size_L2 - 1);
		if (already_sorted < half_buffer_size && already_sorted < num_buffered) {
			// we have not yet completed a half-pass, and we lack an up-to-date incomplete assessment, so do a fresh one
			do_incomplete_buffer();
		}
		total_expected_duplicates += incomplete_expected_duplicates;
		total_actual_duplicates += incomplete_duplicates;
	}
	if (!total_expected_duplicates) return;

	double norm = (total_actual_duplicates - total_expected_duplicates) / std::sqrt(double(total_expected_duplicates));
	std::ostringstream buf;
	buf << get_name();
	if (BirthdayLambda1::expected_duplicates) buf << ":all";
	else buf << ":inc";
	if (total_expected_duplicates < 100) {
		if (norm <= 0) {
			double p = 0;
			for (int i = 0; i < total_actual_duplicates; i++) p += math_poisson_pmf(total_expected_duplicates, i);
			p += 0.5 * math_poisson_pmf(total_expected_duplicates, total_actual_duplicates);
			results.push_back(TestResult(buf.str() + "1", norm, 1 - p, TestResult::TYPE_BAD_P, 0.125));
		}
		else {
			double p = 0;
			double high = total_actual_duplicates > total_expected_duplicates ? total_actual_duplicates : total_expected_duplicates;
			high += 8 * std::sqrt(total_expected_duplicates);

			for (int i = high + 3; i > total_actual_duplicates; i--) p += math_poisson_pmf(total_expected_duplicates, i);
			p += 0.5 * math_poisson_pmf(total_expected_duplicates, total_actual_duplicates);
			results.push_back(TestResult(buf.str() + "1", norm, p, TestResult::TYPE_BAD_P, 0.125));
		}
	}
	else {
		results.push_back(TestResult(buf.str() + "1", norm, math_normaldist_to_pvalue(-norm), TestResult::TYPE_BAD_P, 0.125));
		// would be better passed through calibration, but this will do for now
	}

	buf << "2";
	double norm2 = score / std::sqrt(double(total_expected_duplicates));
	results.push_back(TestResult(buf.str(), norm2, math_normaldist_to_pvalue(-norm2), TestResult::TYPE_BAD_P, 0.125));
}
double PractRand::Tests::BirthdaySystematic128::evaluate_score(double lambda, Uint64 num_duplicates) {
	long SIZE = lambda * 2 + std::sqrt(lambda) * 5 + 5;
	std::vector<double> probs; probs.resize(SIZE);
	for (int i = 0; i < SIZE; i++) probs[i] = math_poisson_pmf(lambda, i);
	double tmp = 0;
	for (int i = 0; i < SIZE; i++) tmp += probs[i];
	for (int i = 0; i < SIZE; i++) probs[i] /= tmp;
	std::vector<double> score_roots; score_roots.resize(SIZE);
	for (int i = 0; i < SIZE; i++) score_roots[i] = -std::log(probs[i]);
	double mean = 0;
	for (int i = 0; i < SIZE; i++) mean += score_roots[i] * probs[i];
	for (int i = 0; i < SIZE; i++) score_roots[i] -= mean;
	double dev = 0;
	for (int i = 0; i < SIZE; i++) dev += score_roots[i] * score_roots[i] * probs[i];
	dev = std::sqrt(dev);
	for (int i = 0; i < SIZE; i++) probs[i] /= dev;

	return (-std::log(math_poisson_pmf(lambda, num_duplicates)) - mean) / dev;
}
Uint64 PractRand::Tests::BirthdaySystematic128::flush_buffer() {
	Uint64 dups = BirthdayLambda1::flush_buffer();
	if (autofail) return dups;
	
	if (expected_duplicates == 1) score = 0;
	score += evaluate_score(1.0, dups);//scoring method 2
	return dups;
}
void PractRand::Tests::BirthdaySystematic128::test_blocks(TestBlock *data, int numblocks) {
	if (autofail) return;
	Uint64 mask_high = Uint64(Sint64(-1)), mask_low;
	if (bits_to_use < 64) {
		mask_low = 0;
		mask_high <<= (64 - bits_to_use);
	}
	else if (bits_to_use == 64) mask_low = 0;
	else if (bits_to_use < 128) {
		mask_low = mask_high << (128 - bits_to_use);
	}
	else if (bits_to_use == 128) mask_low = mask_high;
	else {
		issue_error();
		mask_low = 0;
	}
	while (numblocks) {
		i128 *dest = &buffer[num_buffered];
		Uint64 *cur = &data[0].as64[0];
		Uint64 *end = &data[1].as64[0];
		enum { LOW = 0, HIGH = 1 };//that's kind of endian-ist, but the tests generally don't bother dealing with such issues
		for (; cur != end; cur += 2, dest++) {
			dest->low = cur[LOW] & mask_low;
			dest->high = cur[HIGH] & mask_high;
			int ri = dest->high >> (64 - SORT_HELPER_BITS);
			sort_helper_counts[ri]++;
		}
		num_buffered += TestBlock::SIZE / sizeof(i128);
		if (num_buffered >> buffer_size_L2) {
			flush_buffer();
			mask_high = Uint64(Sint64(-1));
			if (bits_to_use < 64) {
				mask_low = 0;
				mask_high <<= (64 - bits_to_use);
			}
			else if (bits_to_use == 64) mask_low = 0;
			else if (bits_to_use < 128) {
				mask_low = mask_high << (128 - bits_to_use);
			}
		}
		else if (expected_duplicates == 0 && num_buffered == (1ull << (buffer_size_L2 - 1))) {
			do_incomplete_buffer();
		}

		numblocks--;
		data++;
	}
}


PractRand::Tests::BirthdayAlt::BirthdayAlt(int buffer_size_L2_, int filter_bits_) : buffer_size_L2(buffer_size_L2_), filter_bits(filter_bits_) {
	if (buffer_size_L2 < 6 || buffer_size_L2 > 29) issue_error("BirthdayAlt - bad buffer_size_L2");
	buffer = NULL;
}
PractRand::Tests::BirthdayAlt::~BirthdayAlt() {
	delete[] buffer;
	buffer = NULL;
}
void PractRand::Tests::BirthdayAlt::init(PractRand::RNGs::vRNG *known_good) {
	if (!buffer) buffer = new i128[1 << buffer_size_L2];
	num_buffered = 0;
	autofail = false;

	for (int i = 0; i < (1 << SORT_HELPER_BITS); i++) sort_helper_counts[i] = 0;

	score_sum_log = 0;
	score_sum_log2 = 0;
	score_sum_log_sqr = 0;
	count = 0;
}
std::string Tests::BirthdayAlt::get_name() const {
	std::ostringstream buf;
	buf << "BDayX(" << buffer_size_L2 << ")";
	return buf.str();
}
void PractRand::Tests::BirthdayAlt::_lookup_constants(int table_size_L2,long double *_offset, long double *_deviation, long double *_sample_size) {
	struct PerSizeEmpiricalData {
		double mean;
		double dev;
		Uint64 samples;
	};
	if (table_size_L2 < 6) issue_error("BirthdayAlt::_lookup_constants: table_size_L2 too low");
	if (table_size_L2 >= 32) issue_error("BirthdayAlt::_lookup_constants: table_size_L2 too high");
	static const PerSizeEmpiricalData table[32] = {
		// haven't found a good formula for these values yet, but the number of possible cases is limited and we can tolerate minor errors
		// and the units are in an ugly format, will try to fix that shortly
		{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },//0-5
		{ 0.004760612, 0.0102824, 1ull << 26 },//6
		{ 0.005065438, 0.0103139, 1ull << 26 },//7
		{ 0.005284246, 0.0104208, 1ull << 24 },//8
		{ 0.005454068, 0.0105728, 1ull << 24 },//9
		{ 0.005596846, 0.0107491, 1ull << 26 },//10
		{ 0.005726084, 0.0109407, 1ull << 26 },//11
		{ 0.005849366, 0.0111451, 1ull << 26 },//12
		{ 0.005971410, 0.0113600, 1ull << 26 },//13
		{ 0.006094866, 0.0115875, 1ull << 26 },//14
		{ 0.006221473, 0.0118215, 1ull << 24 },//15
		{ 0.006352280, 0.0120681, 1ull << 24 },//16
		{ 0.00648806,  0.012326, 1ull << 20 },//17
		{ 0.00662950,  0.012593, 1ull << 20 },//18
		{ 0.00677697,  0.012870, 1ull << 20 },//19
		{ 0.00693114,  0.013128, 1ull << 17 },//20
		{ 0.00709237,  0.013502, 1ull << 17 },//21
		{ 0.007261390, 0.0137676, 1ull << 19 },//22
		{ 0.00743840,  0.01419, 1ull << 14 },//23
		{ 0.00762432,  0.01453, 1ull << 14 },//24
		{ 0.007819814, 0.0148607, 1ull << 18 },//25
		{ 0.00802560,  0.01510, 1ull << 12 },//26
		{ 0.00824261,  0.01566, 1ull << 12 },//27
		{ 0.00847146,  0.01613, 1ull << 10 },//28
		{ 0.00871354,  0.01682, 1ull << 10 },//29
		{ 0, 0, 0 },//30
		{ 0, 0, 0 },//31
	};
	long double offset, dev, samples;
	if (_offset) {
		if (table[table_size_L2].mean > 0) *_offset = table[table_size_L2].mean;
		else ;
	}
	if (_deviation) *_deviation = table[table_size_L2].dev;
	if (_sample_size) *_sample_size = table[table_size_L2].samples;
}

void PractRand::Tests::BirthdayAlt::get_results(std::vector<TestResult> &results) {
	if (!count) return;
	long buffer_size = 1 << buffer_size_L2;
	long double dev, _sample_size, uncertainty;
	_lookup_constants(buffer_size_L2, NULL, &dev, &_sample_size);

	double score = score_sum_log / std::sqrt(double(count)) / dev;

	results.push_back(TestResult(get_name(), score, score, TestResult::TYPE_RAW_NORMAL, 0.125));
}
void PractRand::Tests::BirthdayAlt::flush_buffer() {
	long buffer_size = 1 << buffer_size_L2;
	if (num_buffered != buffer_size) issue_error("BirthdayAlt::flush_buffer - buffer not full?");
	num_buffered = 0;
	if (autofail) return;
	BirthdayHelpers::_sorted_deltas_of_sorted_values(buffer, buffer_size_L2, sort_helper_counts);
	for (int i = 0; i < (1 << SORT_HELPER_BITS); i++) sort_helper_counts[i] = 0;

	long double expected_log_offset, expected_log_samples, deviation, uncertainty;
	_lookup_constants(buffer_size_L2, &expected_log_offset, &deviation, &expected_log_samples);

	long double sum_log = 0, sum_log2 = 0;
	long double expected_delta1 = std::pow(2.0, 128 - buffer_size_L2);
	//long double expected_log1 = std::log(expected_delta1);
	long double expected_delta2 = expected_delta1 / (buffer_size - 1);
	long double expected_log2 = std::log(expected_delta2);
	long double adjusted_expected_log2 = expected_log2 * (1 + expected_log_offset);
	long double adjusted_expected_delta2 = std::exp(adjusted_expected_log2);

	for (int i = 0; i < buffer_size - 2; i++) {
		i128 diff = buffer[i + 1] - buffer[i];
		long double delta = diff.low + 1.0;
		delta += diff.high * 18446744073709551616.0;
		//long double L = std::log(delta);
		//sum_log -= L - expected_log2;//instead of subtracting expected_log2, change to dividing by expected_delta2
		//long double L = std::log(delta / expected_delta2) - expected_log_offset * expected_log2;
		long double L = std::log(delta / adjusted_expected_delta2);
		sum_log -= L;
		sum_log2 += L * L;
	}
	sum_log /= expected_log2;
	sum_log /= std::sqrt(double(buffer_size - 2));
	//sum_log += expected_log_offset * std::sqrt(double(buffer_size - 2));

	score_sum_log += sum_log;
	score_sum_log2 += sum_log2;
	score_sum_log_sqr += sum_log * sum_log;
	count++;
	if (count == 1) {
		for (int i = 6; i <= 29; i++) {
			long double old;
			_lookup_constants(i, &old, NULL, NULL);
			Uint64 bufsize = 1ull << i;
			long double preadj = std::pow(2.0, 128.0 - i) / (bufsize - 1);
			//long double revised = ;
		}
	}
	if (count > 1 && !(count & (count-1))) {
		std::printf("\n");
		std::printf("sum_log obs ofs: %+.10Lf\n", (long double)score_sum_log / count / std::sqrt(double(buffer_size - 2)) + expected_log_offset);
		std::printf("sum_log exp ofs: %+.10Lf\n", expected_log_offset);
		std::printf("sum_log   delta: %+.12Lf\n", (long double)score_sum_log / count / std::sqrt((long double)(buffer_size - 2)));
		std::printf("revised obs ofs: %+.10Lf\n", score_sum_log / count / std::sqrt(double(buffer_size - 2)) + adjusted_expected_log2 - expected_log2);
		std::printf("revised exp ofs: %+.10Lf\n", expected_log2 * expected_log_offset);
		std::printf("sum_log obs dev: %.9Lf\n", (long double)std::sqrt(score_sum_log_sqr / count - (score_sum_log / count * score_sum_log / count)));
		std::printf("sum_log exp dev: %.9Lf\n", deviation);
		std::printf("count: %.0Lf\n\n", (long double)(count));
		//std::exit(0);
	}
}
void PractRand::Tests::BirthdayAlt::test_blocks(TestBlock *data, int numblocks) {
	blocks_tested += numblocks;
	if (autofail) return;
	if (!filter_bits) {
		while (numblocks) {
			i128 *dest = &buffer[num_buffered];
			Uint64 *cur = &data[0].as64[0];
			Uint64 *end = &data[1].as64[0];
			enum { LOW = 0, HIGH = 1 };//that's kind of endian-ist, but the tests generally don't bother dealing with such issues
			for (; cur != end; cur += 2, dest++) {
				dest->low = cur[LOW];
				dest->high = cur[HIGH];
				int ri = dest->high >> (64 - SORT_HELPER_BITS);
				sort_helper_counts[ri]++;
			}

			num_buffered += TestBlock::SIZE / sizeof(i128);
			if (num_buffered >> buffer_size_L2) {
				flush_buffer();
			}
			numblocks--;
			data++;
		}
	}
	else {
		long shift = 64 - SORT_HELPER_BITS - filter_bits;
		i128 *dest = &buffer[num_buffered];
		i128 *dest_end = &buffer[num_buffered];
		Uint64 *cur = &data[0].as64[0];
		Uint64 *end = &data[numblocks].as64[0];
		enum { LOW = 0, HIGH = 1 };//that's kind of endian-ist, but the tests generally don't bother dealing with such issues
		for (; cur != end; cur += 2) {
			dest->low = cur[LOW];
			dest->high = cur[HIGH];
			Uint64 ri = dest->high >> shift;
			if (ri >= ((1 << SORT_HELPER_BITS) - 1)) continue;
			sort_helper_counts[ri]++;
			dest++;
			if (dest == dest_end) {
				num_buffered = 1ull << buffer_size_L2;
				flush_buffer();
				dest = &buffer[0];
			}
		}
		num_buffered = dest - &buffer[0];
	}
}









PractRand::Tests::Pat5::Pat5()
//:
//	lifespan(1<<16),
{
}
void PractRand::Tests::Pat5::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	for (int pi = 0; pi < (1 << PATTERN_INDEX_BITS); pi++) patterns[pi].total_count = -1;
	blocks_tested = 0;
}
std::string PractRand::Tests::Pat5::get_name() const {
	std::ostringstream tmp;
	tmp << "Pat5";
	return tmp.str();
}
void PractRand::Tests::Pat5::test_blocks(TestBlock *data, int numblocks) {
	enum { CENTER = (PATTERN_WIDTH - 1) / 2 };
	int max = numblocks * (TestBlock::SIZE / sizeof(Word)) - CENTER;
	int min = blocks_tested ? 0 : CENTER;
	for (long i = min; i < max; i++) {
		Word word = data->as32[i];
		if (!(word & ((1<<ZERO_FILTER_BITS)-1))) {
			int pi = word >> ((8 * sizeof(word)) - PATTERN_INDEX_BITS);
			if (patterns[pi].total_count == -1) {
				if (count_ones32(data->as32[i] ^ data->as32[i + 1]) > 4) {
					patterns[pi].total_count = 0;
					for (int j = 0; j < PATTERN_WIDTH; j++) patterns[pi].base_pattern[j] = data->as32[i + j - CENTER];
				}
			}
			else {
				patterns[pi].total_count++;
				int d1 = count_ones32(patterns[pi].base_pattern[CENTER] ^ word);
				unsigned int d2 = 0, d3 = 0, d4 = 0;
				enum { D1 = 1, D2 = D1 + NUM_SECONDARY_WORDS, D3 = D2 + NUM_TERTIARY_WORDS, D4 = D3 + NUM_QUATERNARY_WORDS };
				for (int j = D1; j < D2; j++) {
					d2 += count_ones32(patterns[pi].base_pattern[CENTER - j] ^ data->as32[i - j]);
					d2 += count_ones32(patterns[pi].base_pattern[CENTER + j] ^ data->as32[i + j]);
				}
				for (int j = D2; j < D3; j++) {
					d3 += count_ones32(patterns[pi].base_pattern[CENTER - j] ^ data->as32[i - j]);
					d3 += count_ones32(patterns[pi].base_pattern[CENTER + j] ^ data->as32[i + j]);
				}
				for (int j = D3; j < D4; j++) {
					d4 += count_ones32(patterns[pi].base_pattern[CENTER - j] ^ data->as32[i - j]);
					d4 += count_ones32(patterns[pi].base_pattern[CENTER + j] ^ data->as32[i + j]);
				}
				d1 >>= PRIMARY_WORD_DISTANCE_EXTRA_BITS;
				d2 >>= SECONDARY_WORD_DISTANCE_EXTRA_BITS;
				d3 >>= TERTIARY_WORD_DISTANCE_EXTRA_BITS;
				d4 >>= QUATERNARY_WORD_DISTANCE_EXTRA_BITS;
				enum { MAX_d1 = (1 << PRIMARY_WORD_DISTANCE_BITS) - 1 };
				enum { MAX_d2 = (1 << SECONDARY_WORD_DISTANCE_BITS) - 1 };
				enum { MAX_d3 = (1 << TERTIARY_WORD_DISTANCE_BITS) - 1 };
				enum { MAX_d4 = (1 << QUATERNARY_WORD_DISTANCE_BITS) - 1 };
				if (d1 > MAX_d1) d1 = MAX_d1;
				if (d2 > MAX_d2) d2 = MAX_d2;
				if (d3 > MAX_d3) d3 = MAX_d3;
				if (d4 > MAX_d4) d4 = MAX_d4;
				enum { D1SH = 0, D2SH = D1SH + PRIMARY_WORD_DISTANCE_BITS, D3SH = D2SH + SECONDARY_WORD_DISTANCE_BITS, D4SH = D3SH + TERTIARY_WORD_DISTANCE_BITS };
				unsigned int ci = (d1 << D1SH) +(d2 << D2SH) + (d3 << D3SH) + (d4 << D4SH);
				ci += pi << (TABLE_SIZE_L2 - PATTERN_INDEX_BITS);
				counts.increment(ci);
			}
		}
	}
	blocks_tested += numblocks;
}
static std::vector<double> get_Pat5_prob_sub_table(int base_bits, int shift, int final_bits) {
	std::vector<double> pdf, cdf, result;
	Tests::get_hamming_weight_chances(base_bits, pdf, cdf);
	result.resize(1 << final_bits, 0.0);
	int max = (1 << final_bits) - 1;
	for (int i = 0; i <= base_bits; i++) {
		double chance = (i <= base_bits / 2) ? pdf[i] : pdf[base_bits - i];
		int transform = i >> shift;
		if (transform > max) transform = max;
		result[transform] += chance;
	}
	return result;
}
void PractRand::Tests::Pat5::get_results(std::vector<TestResult> &results) {
	if (!blocks_tested) return;
	const Uint64 *_counts = counts.get_array();
	Uint64 total_opportunities = blocks_tested * TestBlock::SIZE / sizeof(Word) - PATTERN_WIDTH + 1;
	enum { TOTAL_SIZE = 1 << TABLE_SIZE_L2, TOTAL_PATTERNS = 1 << PATTERN_INDEX_BITS };
	enum { BASE_SIZE = TOTAL_SIZE / TOTAL_PATTERNS };
	std::vector<double> base_probs; base_probs.resize(BASE_SIZE);
	if (true) {
		std::vector<double> primary_probs, secondary_probs, tertiary_probs, quaternary_probs;
		std::vector<double> cdf;
		enum { PRIMARY_BITS = WORD_BITS - ZERO_FILTER_BITS - PATTERN_INDEX_BITS };
		primary_probs = get_Pat5_prob_sub_table(PRIMARY_BITS, PRIMARY_WORD_DISTANCE_EXTRA_BITS, PRIMARY_WORD_DISTANCE_BITS);
		secondary_probs = get_Pat5_prob_sub_table(WORD_BITS * 2 * NUM_SECONDARY_WORDS, SECONDARY_WORD_DISTANCE_EXTRA_BITS, SECONDARY_WORD_DISTANCE_BITS);
		tertiary_probs = get_Pat5_prob_sub_table(WORD_BITS * 2 * NUM_TERTIARY_WORDS, TERTIARY_WORD_DISTANCE_EXTRA_BITS, TERTIARY_WORD_DISTANCE_BITS);
		quaternary_probs = get_Pat5_prob_sub_table(WORD_BITS * 2 * NUM_QUATERNARY_WORDS, QUATERNARY_WORD_DISTANCE_EXTRA_BITS, QUATERNARY_WORD_DISTANCE_BITS);
		for (int i = 0; i < BASE_SIZE; i++) {
			double chance = 1.0;
			enum { D1SH = 0, D2SH = D1SH + PRIMARY_WORD_DISTANCE_BITS, D3SH = D2SH + SECONDARY_WORD_DISTANCE_BITS, D4SH = D3SH + TERTIARY_WORD_DISTANCE_BITS };
			chance *= primary_probs[(i >> D1SH) & ((1 << PRIMARY_WORD_DISTANCE_BITS) - 1)];
			chance *= secondary_probs[(i >> D2SH) & ((1 << SECONDARY_WORD_DISTANCE_BITS) - 1)];
			chance *= tertiary_probs[(i >> D3SH) & ((1 << TERTIARY_WORD_DISTANCE_BITS) - 1)];
			chance *= quaternary_probs[(i >> D4SH) & ((1 << QUATERNARY_WORD_DISTANCE_BITS) - 1)];
			base_probs[i] = chance;
		}
	}
	if (true) {
		enum { INCLUDE_NON_MATCHES = 0 };// 0 or 1
		std::vector<Uint64> counts2; counts2.resize(TOTAL_SIZE + INCLUDE_NON_MATCHES);
		std::vector<double> probs2; probs2.resize(TOTAL_SIZE + INCLUDE_NON_MATCHES);
		std::copy(&_counts[0], &_counts[TOTAL_SIZE], &counts2[0]);
		double any_match_prob = std::pow(0.5, INCLUDE_NON_MATCHES ? ZERO_FILTER_BITS : 0);
		double specific_match_prob = any_match_prob * std::pow(0.5, PATTERN_INDEX_BITS);
		if (INCLUDE_NON_MATCHES) probs2[TOTAL_SIZE] = 1 - any_match_prob;
		for (int i = 0; i < TOTAL_SIZE; i++) probs2[i] = base_probs[i & (BASE_SIZE - 1)] * specific_match_prob;
		Sint64 total_matches = 0;
		for (int i = 0; i < TOTAL_PATTERNS; i++) total_matches += patterns[i].total_count;
		if (total_matches < 100) return;
		if (INCLUDE_NON_MATCHES) counts2[TOTAL_SIZE] = total_opportunities - total_matches;
		if (total_opportunities > 300) {
			double rarity = Tests::rarity_test(TOTAL_SIZE + INCLUDE_NON_MATCHES, &probs2[0], &counts2[0], false, true);
			std::ostringstream ss; ss << get_name() << "(*,r)"; results.push_back(TestResult(ss.str(), rarity, 0, TestResult::TYPE_RAW_NORMAL, 0.125));
		}
		if (total_opportunities > 3000) {
			int n = Tests::simplify_prob_table(TOTAL_SIZE + INCLUDE_NON_MATCHES, (INCLUDE_NON_MATCHES ? total_opportunities : total_matches) / 40.0, &probs2[0], &counts2[0], false, false);
			double raw = Tests::g_test(n, &probs2[0], &counts2[0]);
			double norm = Tests::math_chisquared_to_normal(raw, n - 1);
			std::ostringstream ss; ss << get_name() << "(*,g)"; results.push_back(TestResult(ss.str(), norm, 0, TestResult::TYPE_RAW_NORMAL, 0.125));
		}
	}
	for (int pi = 0; pi < TOTAL_PATTERNS; pi++) {
		if (patterns[pi].total_count < 30) continue;
		std::vector<double> local_probs = base_probs;
		std::vector<Uint64> local_counts; local_counts.resize(BASE_SIZE);
		std::copy(&_counts[pi * BASE_SIZE], &_counts[pi * BASE_SIZE + BASE_SIZE], &local_counts[0]);
		if (true) {
			double rarity = Tests::rarity_test(BASE_SIZE, &local_probs[0], &local_counts[0], false, true);
			std::ostringstream ss; ss << get_name() << "(" << pi << ",r)"; results.push_back(TestResult(ss.str(), rarity, 0, TestResult::TYPE_RAW_NORMAL, 0.001 / TOTAL_PATTERNS));
		}
		if (patterns[pi].total_count > 300) {
			int n = Tests::simplify_prob_table(BASE_SIZE, patterns[pi].total_count / 40.0, &local_probs[0], &local_counts[0], false, false);
			double raw = Tests::g_test(n, &local_probs[0], &local_counts[0]);
			double norm = Tests::math_chisquared_to_normal(raw, n - 1);
			std::ostringstream ss; ss << get_name() << "(" << pi << ",g)"; results.push_back(TestResult(ss.str(), norm, 0, TestResult::TYPE_RAW_NORMAL, 0.001 / TOTAL_PATTERNS));
		}
	}

	//finishing
/*	int reduced_size = simplify_prob_table(size,
		blocks_tested * (TestBlock::SIZE >> unitsL) / 25.0,
		&probs[0], &tmp_counts[0], true, true);
	double r = g_test(reduced_size, &probs[0], &tmp_counts[0]);
	r = math_chisquared_to_normal(r, reduced_size - 1);
	double weight = std::pow(2.0, 1.0 - unitsL / 2.0);
	if (unitsL != 0) weight *= 0.75;
	if (size < 1024 * 128) weight *= 0.5;
	Uint64 min_len = calibration_manager.get_minimum_length(get_name());
	if (min_len && min_len <= blocks_tested) {
		TestCalibrationData *calib = calibration_manager.get_calibration_data(get_name(), blocks_tested);
		double suspicion = calib->sample_to_suspicion(r) * -1;//negation to make the normal failure type occur at 0 instead of 1
		results.push_back(TestResult(get_name(), r, suspicion, TestResult::TYPE_GOOD_S, weight));
	}
	else if (blocks_tested > unitsL * 1024 * 1024 * 16) {
		results.push_back(TestResult(get_name(), r, 0, TestResult::TYPE_RAW_NORMAL, weight / 5));
	}
	else {
		results.push_back(TestResult(get_name(), r, 0, TestResult::TYPE_RAW, .01));
	}*/
}





void PractRand::Tests::CoupGap::init( RNGs::vRNG *known_good ) {
	autofail = 0;
	blocks_tested = 0;
	symbols_ready = 0;
	for (int i = 0; i < 256; i++) {
		sym_has_appeared[i] = false;
		next_younger_sym[i] = (i+1) & 255;
		youngest_sym = 255;
		oldest_sym = 0;
//		last_sym_pos[i] = 0;
	}
	count_syms_by_oldest_sym.reset_counts();
//	count_gaps_by_oldest_sym.reset_counts();
}
std::string PractRand::Tests::CoupGap::get_name( ) const {
	return std::string("CoupGap");
}
void PractRand::Tests::CoupGap::get_results(std::vector<TestResult> &results) {
	if (autofail) {
		results.push_back(TestResult(this->get_name() + ":!", autofail, autofail, TestResult::TYPE_PASSFAIL, 0.0000001));
		return;
	}

	std::vector<double> probs;

	if (blocks_tested >= 256) {
		probs.resize(65536);
		for (int i = 0; i < 65536; i++) probs[i] = 1 / 65536.0;
		const Uint64 *counts_ = count_syms_by_oldest_sym.get_array();
		const double *probs_ = &probs[0];
		double raw = g_test(65536, probs_, counts_);
		raw = (raw - 3 * 65536) / (256 * 32);
		TestCalibrationData *calib = calibration_manager.get_calibration_data("CoupGap:SxO", blocks_tested);
		double suspicion = calib->sample_to_suspicion(raw) * -1;//negation to make the normal failure type occur at 0 instead of 1
		results.push_back(TestResult(get_name() + ":SxO", raw, suspicion, TestResult::TYPE_GOOD_S, 0.25));
	}
}
void PractRand::Tests::CoupGap::test_blocks(TestBlock *data, int numblocks) {
	if (autofail) return;
	int i;
	Uint32 ofs = Uint32(blocks_tested) * TestBlock::SIZE;
	int max = TestBlock::SIZE * numblocks;
	for (i = 0; i < max; i++, ofs++) {
		unsigned long sym = data[0].as8[i];
//		Uint32 last_pos = last_sym_pos[sym];

		if (symbols_ready == 256) {
//			Uint32 oldest_age = ofs - last_sym_pos[oldest_sym] - 256;
//			if (oldest_age > MAX_OLDEST_AGE-1) oldest_age = MAX_OLDEST_AGE-1;
//			Uint32 current_age = ofs - last_pos - 1;
//			if (current_age > MAX_CURRENT_AGE-1) current_age = MAX_CURRENT_AGE-1;
			count_syms_by_oldest_sym.increment(sym * 256 + oldest_sym);
//			count_gaps_by_oldest_sym.increment(oldest_sym + (current_age << 8));
		}
		else if (!sym_has_appeared[sym]) {
			sym_has_appeared[sym] = true;
			symbols_ready ++;
		}

		//update linked list
//		last_sym_pos[sym] = ofs;
		if (oldest_sym == sym) {
			oldest_sym = next_younger_sym[sym];
		}
		next_younger_sym[youngest_sym] = Uint8(sym);
		youngest_sym = sym;
	}
	Uint64 oblocks = blocks_tested;
	blocks_tested += numblocks;
//	if ((oblocks>>17) != (blocks_tested>>17)) {//once every 128 megabytes or so... prevent overflow
//		for (i = 0; i < 256; i++) {
//			int n = last_sym_pos[i];
//			if (n == -1) continue;
//			if ((ofs - n) > 12345678) {
//				autofail += 1;
//			}
//		}
//	}
}






PractRand::Tests::NearSeq::NearSeq() {
	//if (false);
	//else if (BITS_PER_BLOCK == 8) verify_NearSeq_byte_code8(NearSeq_byte_code8);
	//else if (BITS_PER_BLOCK == 4) verify_NearSeq_byte_code4x2(NearSeq_byte_code4x2);
	//else if (BITS_PER_BLOCK == 5) verify_NearSeq_byte_code5x2(NearSeq_byte_code5x2);
	//else issue_error("NearSeq: what block size?");
	lookup_table = NULL;
}
void PractRand::Tests::NearSeq::init(PractRand::RNGs::vRNG *known_good) {
	if (SEQUENCE_BITS < CORE_SEQUENCE_BITS) issue_error("NearSeq - bad settings");
	static const Word lookup[2] = { (1 << BITS_PER_BLOCK) - 1, 0x0 };//backwards, as we're inverting values here for maximum hamming distance
	for (int bi = 0; bi < NUM_BUCKETS; bi++) {
		for (int x = 0; x < SEQUENCE_WORDS; x++) buckets[bi].sequence[x] = 0;
		for (int x = 0; x < BLOCKS_PER_CORE; x++) {
			int dest_bit_pos = x * BITS_PER_BLOCK;
			int dest_word_start = dest_bit_pos >> WORD_BITS_L2;
			int dest_word_end = (dest_bit_pos + BITS_PER_BLOCK - 1) >> WORD_BITS_L2;
			int dest_word_offset = dest_bit_pos & (WORD_BITS - 1);
			//here we're trying to initialize all buckets to impossible values, so that we can recognize them as not-yet-populated later
			buckets[bi].sequence[SEQUENCE_WORD_OFFSET + dest_word_start] |= lookup[(bi >> x) & 1] << dest_word_offset;
			if (dest_word_end != dest_word_start) {//maximum of 1 extra word, since BITS_PER_BLOCK is guaranteed to be less than or equal to WORD_BITS
				buckets[bi].sequence[SEQUENCE_WORD_OFFSET + dest_word_end] |= lookup[(bi >> x) & 1] >> (WORD_BITS - dest_word_offset);
			}
		}
	}

	for (int x = 0; x < MAX_CORE_DISTANCES; x++) core_distances[x] = 0;
	for (int x = 0; x < MAX_CORE_DISTANCES; x++) sum_extra_distances[x] = 0;

	int table_size = 1 << BITS_PER_BLOCK;
	lookup_table = new Uint8[table_size];
	lookup_table2 = new Uint8[table_size];
	for (int i = 0; i < table_size; i++) {
		int bits = count_ones16(i);
		if (bits > BITS_PER_BLOCK / 2) {
			bits = BITS_PER_BLOCK - bits;
			lookup_table[i] = 1;
		}
		else lookup_table[i] = 0;
		if (bits > MAX_ERRORS_PER_BLOCK) lookup_table[i] |= 128;

		if (bits > 7) bits = 7;
		lookup_table2[i] = 1 << bits;
	}

	TestBaseclass::init(known_good);
}
void PractRand::Tests::NearSeq::deinit() {
	delete[] lookup_table;
	delete[] lookup_table2;
	TestBaseclass::deinit();
}
std::string PractRand::Tests::NearSeq::get_name() const {
	return "NearS";
	//std::ostringstream name;
	//name << "BRank(" << rate_hl2 << ")";
	//return name.str();
}
void PractRand::Tests::NearSeq::get_results(std::vector<TestResult> &results) {
	//if (blocks_tested < 1 << 24) return;

	Uint64 total_count = 0;
	for (int i = 0; i < MAX_CORE_DISTANCES; i++) total_count += core_distances[i];
	if (!total_count) return;
	Uint64 target_threshold = (total_count * 9) / 10;
	Uint64 total_so_far = 0;
	int run_high = MAX_CORE_DISTANCES - 1;
	while (true) {
		Uint64 run_sum = sum_extra_distances[run_high];
		Uint64 run_count = core_distances[run_high];
		if (!run_count) {
			run_high--;
			continue;
		}
		int run_low = run_high;
		for (int i = run_high - 1; i >= 0 && run_count < target_threshold; i--) {
			if (core_distances[i]) {
				run_low = i;
				run_count += core_distances[i];
				run_sum += sum_extra_distances[i];
			}
		}
		Uint64 total_bits = run_count * (SEQUENCE_BITS - CORE_SEQUENCE_BITS);
		double variance = 0.5 * 0.5 * total_bits;
		double value = ((Sint64(run_sum << 1) - Sint64(total_bits)) * 0.5) / std::sqrt(variance);
		std::ostringstream os;
		os << get_name() << ":[" << run_low << "-" << run_high << "](" << run_count << ")";
		results.push_back(TestResult(os.str(), value, math_normaldist_to_suspicion(-value), total_bits > 4000 ? TestResult::TYPE_GOOD_S : TestResult::TYPE_BAD_S, 0.1));

		total_so_far += run_count;
		if (total_so_far == total_count) return;
		target_threshold = ((total_count - total_so_far) * 9) / 10;
		run_high = run_low - 1;
	}

	/*double leftover_sum = 0;
	Uint64 leftover_count = 0;
	int lowest = 0;
	for (int i = 0; i < MAX_CORE_DISTANCES; i++) {
		enum {EXTRA_BITS = SEQUENCE_BITS - CORE_SEQUENCE_BITS };
		double sum = leftover_sum + sum_extra_distances[i];//binomial distribution
		Uint64 count = leftover_count + core_distances[i];
		double total_bits = count * EXTRA_BITS;
		double variance = 0.5 * 0.5 * total_bits;
		if (variance < 1) variance = 1;
		double mean = total_bits * 0.5;
		double value = (sum - mean) / std::sqrt(variance);
		if (total_bits > 10000) {
			std::ostringstream os;
			os << get_name() << ":[" << lowest << "-" << i << "](" << count << ")";
			results.push_back(TestResult(os.str(), value, math_normaldist_to_suspicion(-value), TestResult::TYPE_GOOD_S, 0.1));
			leftover_sum = 0;
			leftover_count = 0;
			lowest = i + 1;
		}
		else if (total_bits > 0) {
			leftover_sum = sum;
			leftover_count = count;
		}
		else {
			lowest = i + 1;
		}
	}*/
	return;
}
int PractRand::Tests::NearSeq::is_core_good(const Word *core) const {
	int index = 0;
	int worst = 0;
	Word w = core[0];
	int usable_bits_left = WORD_BITS;
	while (index < BLOCKS_PER_CORE) {
		if (usable_bits_left < BITS_PER_BLOCK) {
			int ibit = index * BITS_PER_BLOCK + usable_bits_left;
			int wi = ibit >> WORD_BITS_L2;
			Word w2 = core[wi];
			int boff = ibit & (WORD_BITS - 1);
			if (boff) {
				w2 >>= boff;
				w2 |= core[wi + 1] << (WORD_BITS - boff);
			}
			w |= w2 << usable_bits_left;
			usable_bits_left = WORD_BITS;
		}
		worst |= lookup_table2[w & ((1 << BITS_PER_BLOCK) - 1)];
		w >>= BITS_PER_BLOCK;
		usable_bits_left -= BITS_PER_BLOCK;
		index++;
	}
	return worst < (1 << (GOOD_ERRORS_PER_BLOCK + 1));
}
int PractRand::Tests::NearSeq::core_to_index(const Word *core) const {
	int index = 0;
	Uint8 flags = 0;
	Word w = core[0];
	enum {LOOP1_MAX = WORD_BITS > CORE_SEQUENCE_BITS ? BLOCKS_PER_CORE : WORD_BITS / BITS_PER_BLOCK };
	int obit = 0;
	for (; obit < LOOP1_MAX; obit++) {
		flags |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)];
		index |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)] << obit;
		w >>= BITS_PER_BLOCK;
	}
	if (flags & 128) return -1;
	if (CORE_SEQUENCE_BITS <= WORD_BITS) return index;
	if (WORD_BITS % BITS_PER_BLOCK) {//they don't evenly divide
		int usable_bits_left = BITS_PER_BLOCK % WORD_BITS;
		while (obit < BLOCKS_PER_CORE) {
			if (usable_bits_left < BITS_PER_BLOCK) {
				if (flags & 128) return -1;
				int ibit = obit * BITS_PER_BLOCK + usable_bits_left;
				int wi = ibit >> WORD_BITS_L2;
				Word w2 = core[wi];
				int boff = ibit & (WORD_BITS - 1);
				if (boff) {
					w2 >>= boff;
					w2 |= core[wi + 1] << (WORD_BITS - boff);
				}
				w |= w2 << usable_bits_left;
				usable_bits_left = WORD_BITS;
			}
			flags |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)];
			index |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)] << obit++;
			w >>= BITS_PER_BLOCK;
			usable_bits_left -= BITS_PER_BLOCK;
		}
		if (flags & 128) return -1;
		return index;
	}
	else { // they do evenly divide
		enum {BLOCKS_PER_WORD = WORD_BITS / BITS_PER_BLOCK};
		enum {MAX_FULL_WORDS = BLOCKS_PER_CORE / BLOCKS_PER_WORD};
		for (int word = 1; word < MAX_FULL_WORDS; word++) {
			w = core[word];
			for (int i = 0; i < BLOCKS_PER_WORD; i++) {
				flags |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)];
				index |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)] << obit++;
				w >>= BITS_PER_BLOCK;
			}
			if (flags & 128) return -1;
		}
		if (!(BLOCKS_PER_CORE % BLOCKS_PER_WORD)) return index;
		w = core[(CORE_SEQUENCE_BITS - 1) / WORD_BITS];
		for (int i = 0; i < (BLOCKS_PER_CORE % BLOCKS_PER_WORD); i++) {
			flags |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)];
			index |= lookup_table[w & ((1 << BITS_PER_BLOCK) - 1)] << obit++;
			w >>= BITS_PER_BLOCK;
		}
		if (flags & 128) return -1;
		return index;
	}
}
int  PractRand::Tests::NearSeq::get_core_distance(const Word *core, int bucket_index) const {
	typedef int(*CBFUNC)(Word);
	CBFUNC _count_ones = (WORD_BITS == 64) ? (CBFUNC)count_ones64 : ((WORD_BITS == 32) ? (CBFUNC)count_ones32 : ((WORD_BITS == 16) ? (CBFUNC)count_ones16 : NULL));// ((WORD_BITS == 8) ? counts_bits8 : NULL)));
	int bits_left = CORE_SEQUENCE_BITS;
	int core_distance = 0;
	int pos = 0;
	while (bits_left >= WORD_BITS) {
		core_distance += _count_ones(core[pos] ^ buckets[bucket_index].sequence[pos + SEQUENCE_WORD_OFFSET]);
		bits_left -= WORD_BITS;
		pos++;
	}
	if (CORE_SEQUENCE_BITS % WORD_BITS) {
		Word delta = core[pos] ^ buckets[bucket_index].sequence[pos + SEQUENCE_WORD_OFFSET];
		delta &= Word((1ull << (CORE_SEQUENCE_BITS % WORD_BITS)) - 1);
		core_distance += _count_ones(delta);
	}
	return core_distance;
}
int  PractRand::Tests::NearSeq::get_extra_distance(const Word *core, int bucket_index) const {
	typedef int(*CBFUNC)(Word);
	CBFUNC _count_ones = (WORD_BITS == 64) ? (CBFUNC)count_ones64 : ((WORD_BITS == 32) ? (CBFUNC)count_ones32 : ((WORD_BITS == 16) ? (CBFUNC)count_ones16 : NULL));// ((WORD_BITS == 8) ? counts_bits8 : NULL)));
	int bits_left = CORE_SEQUENCE_BITS;
	int extra_distance = 0;
	int pos = 0;
	//early words
	for (int i = 0; i < SEQUENCE_WORD_OFFSET; i++) extra_distance += _count_ones(core[i - SEQUENCE_WORD_OFFSET] ^ buckets[bucket_index].sequence[i]);
	if (CORE_SEQUENCE_BITS % WORD_BITS) {
		enum { INDEX = CORE_SEQUENCE_BITS / WORD_BITS };
		extra_distance += _count_ones((core[INDEX] ^ buckets[bucket_index].sequence[INDEX + SEQUENCE_WORD_OFFSET]) >> (CORE_SEQUENCE_BITS % WORD_BITS));
	}
	for (int i = (CORE_SEQUENCE_BITS + WORD_BITS - 1) / WORD_BITS; i < SEQUENCE_BITS / WORD_BITS - SEQUENCE_WORD_OFFSET; i++) 
		extra_distance += _count_ones(core[i] ^ buckets[bucket_index].sequence[i + SEQUENCE_WORD_OFFSET]);
	return extra_distance;
}
void PractRand::Tests::NearSeq::test_blocks(TestBlock *data, int numblocks) {
	int start = blocks_tested ? -SEQUENCE_WORD_OFFSET : SEQUENCE_WORD_OFFSET;
	int end = (numblocks * (TestBlock::SIZE * 8 / WORD_BITS)) - (SEQUENCE_BITS / WORD_BITS - SEQUENCE_WORD_OFFSET);

	static const Word lookup_block_value[2] = { (1 << BITS_PER_BLOCK) - 1, 0x0 };//backwards

	for (int pos = start; pos < end; pos++) {
		/* chances at various settings
		setting:			64/8		64/8		64/4		32/4		32/8		32/8		16/8		16/8		16/4		desired
		codable/ideal		3/2			2/1			1/0			1/0			2/1			1/0			2/1			1/0			1/0			--
		bucketing bits		8			8			16			8			4			4			2			2			4			6 - 16
		is_word_codable		12.9		20.5 K		1845		43			143			41 K		12			202			6.5			1, or 20 - 50 K
		is_near_ideal		20.5 K		1.67 B		281 T		16.7 M		40.9 K		268 M		202			16 K		4.0 K		10 K - 1 B
		looks nice?						.						*			.			.						

		not perfectly aligned, no tails:

		setting:            60/6		30/6		60/5		63/3		30/3		30/5		63/7		63/9		63/9
		codable/ideal		2/1			2/1			2/1			1/0			1/0			2/1			2/1			3/2			3/1
		bucketing bits		10			5			12			21			10			6			9			7			7
		is_word_codable     42.4		6.5			1			1			1			1			1.2 K		115			115
		is_near_ideal       4.0 M		2.0 K		129 K		4.4 T		2 M			719			268 M		165 K		7.2 B
		looks nice?			*			.			*			.									**			*

		maybe some larger area, or at least not as round numbers? no tails.

		setting:            120/10	120/10	120/12	120/12	96/12	99/9	156/12	156/12	112/8	40/4
		codable/ideal		3/2		4/3		4/3		4/2		4/2		3/2		5/4		4/3		3/2		1/0
		bucketing bits		12		12		10		10		8		11		13		13		14		10
		is_word_codable     367 K	16.8	13 K	13 K	~2 K	1.7 K	27.8	224 K	87.5	110
		is_near_ideal       341 B	367 K	227 M	137 T	204 B	159 M	224 K	73 B	35 M	1 B
		looks nice?					.		*						*						*		*
		*/
		Word *core;
		if (false); 
		else if (WORD_BITS == 8) core = (Word*)&data[0].as8[pos];
		else if (WORD_BITS == 16) core = (Word*)&data[0].as16[pos];
		else if (WORD_BITS == 32) core = (Word*)&data[0].as32[pos];
		else if (WORD_BITS == 64) core = (Word*)&data[0].as64[pos];
		else issue_error("NearS - what word size???");
		int bucket_index = core_to_index(core);
		if (bucket_index < 0) continue;
		Bucket &bucket = buckets[bucket_index];
		if ((bucket.sequence[SEQUENCE_WORD_OFFSET] & ((1 << BITS_PER_BLOCK) - 1)) == lookup_block_value[bucket_index & 1]) {
			if (!is_core_good(core)) continue;
			for (int x = 0; x < SEQUENCE_WORDS; x++) bucket.sequence[x] = core[x - SEQUENCE_WORD_OFFSET];
		}
		else {//a populated bucket
			int core_distance = get_core_distance(core, bucket_index);
			int extra_distance = get_extra_distance(core, bucket_index);
			if (core_distance >= MAX_CORE_DISTANCES) core_distance = MAX_CORE_DISTANCES - 1;
			core_distances[core_distance]++;
			sum_extra_distances[core_distance] += extra_distance;
		}
	}
	blocks_tested += numblocks;
}

static std::vector<double> convolve_signals(std::vector<double> s1, std::vector<double> s2) {
	std::vector<double> rv;
	if (!s1.size() || !s2.size()) return rv;
	rv.resize(s1.size() + s2.size() - 1);
	for (int x = 0; x < s1.size(); x++) {
		for (int y = 0; y < s2.size(); y++) {
			rv[x + y] += s1[x] * s2[y];
		}
	}
	return rv;
}
PractRand::Tests::NearSeq2::NearSeq2() {
	//if (false);
	//else if (BITS_PER_BLOCK == 8) verify_NearSeq_byte_code8(NearSeq_byte_code8);
	//else if (BITS_PER_BLOCK == 4) verify_NearSeq_byte_code4x2(NearSeq_byte_code4x2);
	//else if (BITS_PER_BLOCK == 5) verify_NearSeq_byte_code5x2(NearSeq_byte_code5x2);
	//else issue_error("NearSeq: what block size?");
	lookup_table1 = NULL;
	lookup_table2 = NULL;
	lookup_table3 = NULL;

	if (MAX_HDIST_PER_BLOCK * 2 >= BITS_PER_BLOCK) issue_error();
	for (int i = 0; i <= MAX_HDIST_PER_BLOCK; i++) base_chances[i] = math_nChooseR(BITS_PER_BLOCK, i);
	total_base_chances = 0;
	for (int i = 0; i <= MAX_HDIST_PER_BLOCK; i++) total_base_chances += base_chances[i];
	prob_of_valid_block = total_base_chances * std::pow(0.5, BITS_PER_BLOCK - 1);
	prob_of_valid_core = std::pow(prob_of_valid_block, BLOCKS_PER_CORE);
	for (int i = 0; i <= MAX_HDIST_PER_BLOCK; i++) base_probs[i] = base_chances[i] * std::pow(0.5, BITS_PER_BLOCK - 1);
	for (int i = 0; i <= MAX_HDIST_PER_BLOCK; i++) normalized_base_probs[i] = base_chances[i] / (long double)total_base_chances;

	if (true) {
		std::vector<double> base;
		base.resize(MAX_HDIST_PER_BLOCK + 1);
		for (int i = 0; i <= MAX_HDIST_PER_BLOCK; i++) base[i] = base_probs[i];
		std::vector<double> current = base;
		for (int i = 1; i < BLOCKS_PER_CORE; i++) {
			current = convolve_signals(current, base);
		}
		if (current.size() != MAX_TOTAL_HDIST + 1) issue_error();
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) full_probs[i] = current[i];
	}
	if (true) {
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_scores[i] = -std::log(full_probs[i] / prob_of_valid_core);
		double avg = 0;
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) avg += hdist_scores[i] * full_probs[i] / prob_of_valid_core;
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_scores[i] -= avg;
		double avg_sqr = 0;
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) avg_sqr += hdist_scores[i] * hdist_scores[i] * full_probs[i] / prob_of_valid_core;
		double stddev = std::sqrt(avg_sqr);
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_scores[i] /= stddev;
	}

	if (EXTRA1_FULL_WORDS & 1) issue_error("NearSeq2 - odd number of extra words");
	if (EXTRA1_PARTIAL_WORD_BITS < 0 || EXTRA1_PARTIAL_WORD_BITS >= WORD_BITS) issue_error("NearSeq2 - EXTRA1_PARTIAL_WORD_BITS value outside of range");
	if (BITS_PER_BLOCK > 64) issue_error("NearSeq2 - blocks too large");
	if (MAX_HDIST_PER_BLOCK * 2 > BITS_PER_BLOCK) issue_error("NearSeq2 - noise tolerance exceeds maximum possible");
}
void PractRand::Tests::NearSeq2::init(PractRand::RNGs::vRNG *known_good) {
	//_total_cores = 0;
	//_total_valid_cores = 0;
	//_total_invalid_cores = 0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
	if (!lookup_table1) {
		lookup_table3 = new Uint8[MAX_TOTAL_HDIST + 1];
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) {
			int bin = (MAX_TOTAL_HDIST - i) >> HDIST_BIN_SHIFT;
			if (bin >= HDIST_BINS) bin = HDIST_BINS - 1;
			lookup_table3[i] = HDIST_BINS - 1 - bin;
		}
		for (int i = 0; i < HDIST_BINS; i++) hdist_bin_probs[i] = 0;
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_bin_probs[get_hdist_bin(i)] += full_probs[i];

		if (BITS_PER_BLOCK <= MAX_LOOKUP_L2) {//index directly with block value
			lookup_table1 = new Sint8[1 << BITS_PER_BLOCK];
			lookup_table2 = new Uint8[1 << BITS_PER_BLOCK];
			for (int i = 0; i < (1 << BITS_PER_BLOCK); i++) {
				int h = count_ones32(i);
				int v1, v2;
				if (h >= BITS_PER_BLOCK - MAX_HDIST_PER_BLOCK) {
					v1 = 1;
					v2 = BITS_PER_BLOCK - h;
				}
				else if (h > MAX_HDIST_PER_BLOCK) {
					v1 = 1 << 7;
					v2 = 255;
				}
				else {
					v1 = 0;
					v2 = h;
				}
				lookup_table1[i] = v1;
				lookup_table2[i] = v2;
			}
		}
		else {//block value too large, index with block hamming weight instead
			lookup_table1 = new Sint8[BITS_PER_BLOCK + 1];
			lookup_table2 = new Uint8[BITS_PER_BLOCK + 1];
			for (int h = 0; h <= BITS_PER_BLOCK; h++) {
				int v1, v2;
				if (h >= BITS_PER_BLOCK - MAX_HDIST_PER_BLOCK) {
					v1 = 1;
					v2 = BITS_PER_BLOCK - h;
				}
				else if (h > MAX_HDIST_PER_BLOCK) {
					v1 = 1 << 7;
					v2 = 255;
				}
				else {
					v1 = 0;
					v2 = h;
				}
				lookup_table1[h] = v1;
				lookup_table2[h] = v2;
			}
		}
	}
#pragma GCC diagnostic pop
	for (int i = 0; i < NUM_BUCKETS; i++) buckets[i].reset();
	TestBaseclass::init(known_good);
}
void PractRand::Tests::NearSeq2::Bucket::reset() {
	for (int i = 0; i <= MAX_TOTAL_HDIST; i++) core_hdist[i] = 0;
	for (int b = 0; b < HDIST_BINS; b++) for (int i = 0; i < EXTRA1_BITS; i++) extra_counts[b][i] = 0;
}
void PractRand::Tests::NearSeq2::deinit() {
	delete[] lookup_table1;
	delete[] lookup_table2;
	lookup_table1 = NULL;
	lookup_table2 = NULL;
	TestBaseclass::deinit();
}
std::string PractRand::Tests::NearSeq2::get_name() const {
	return "NearS2";
}
bool PractRand::Tests::NearSeq2::is_word_bad(Word word) const {
	// this should only be used in cases where there is an integer number of blocks per word, and zero extra partial word bits
	if (BITS_PER_BLOCK > WORD_BITS) issue_error();
	if (WORD_BITS % BITS_PER_BLOCK) issue_error();
	if (EXTRA1_PARTIAL_WORD_BITS) issue_error();
	Sint8 is_bad = 0;
	for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
		Word block = word;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
		if (BITS_PER_BLOCK < WORD_BITS) block &= Word((1ull << BITS_PER_BLOCK) - 1);
		if (BITS_PER_BLOCK * 2 <= WORD_BITS) word >>= BITS_PER_BLOCK;
#pragma GCC diagnostic pop
		is_bad |= lookup1(block);
	}
	return is_bad < 0;
}
bool PractRand::Tests::NearSeq2::is_core_bad(const Word *core) const {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
	Sint8 is_bad = 0;
	if (CORE_WORDS == 1) {
		Word w = core[0];
		for (int i = 0; i < BLOCKS_PER_CORE; i++) {
			is_bad |= lookup1(w);
			w >>= BITS_PER_BLOCK;
			if (CHECK_VALIDITY_EARLY && is_bad < 0) return true;
		}
		if (is_bad < 0) return true;
		else return false;
	}
	else if (!(WORD_BITS % BITS_PER_BLOCK)) {//blocks align to word boundaries
		int index;
		Word w = core[0];
		for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
			if (lookup1(w) < 0) return true;
			w >>= BITS_PER_BLOCK;
		}
		for (index = 1; index < CORE_WORDS - 1; index++) {
			w = core[index];
			for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
				is_bad |= lookup1(w);
				w >>= BITS_PER_BLOCK;
			}
			if (is_bad < 0) return true;
		}
		w = core[index];
		for (int i = 0; i < BLOCKS_PER_CORE - (CORE_WORDS - 1) * (WORD_BITS / BITS_PER_BLOCK); i++) {
			is_bad |= lookup1(w);
			w >>= BITS_PER_BLOCK;
		}
		if (is_bad < 0) return true;
		else return false;
	}
	else {//blocks do NOT align to word boundaries
		//
		// ...this is ugly, so if possible choose parameterizations that don't hit this
		// or, failing that, paramterizations that will reject most cores in the first word
		Word w = core[0];
		for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {//first word
			is_bad |= lookup1(w);
			w >>= BITS_PER_BLOCK;
		}
		if (is_bad & 128) return true;
		int index = 1;
		enum { WORD_LEFTOVERS = WORD_BITS % BITS_PER_BLOCK };
		int usable_bits = WORD_LEFTOVERS;
		while (index < CORE_WORDS - 1) {//middle words
			Word w2 = core[index++];
			w |= w2 << usable_bits;
			for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
				is_bad |= lookup1(w);
				w >>= BITS_PER_BLOCK;
			}
			if (is_bad < 0) return true;
			usable_bits += WORD_LEFTOVERS;
			w = w2 >> (WORD_BITS - usable_bits);
			if (usable_bits >= BITS_PER_BLOCK) {
				is_bad |= lookup1(w);
				w >>= BITS_PER_BLOCK;
				usable_bits -= BITS_PER_BLOCK;
			}
		}
		//last word
		enum { INDEX = CORE_WORDS - 1 };
		if (index != INDEX) issue_error("NearSeq2::is_core_good - internal error");
		enum { USABLE_BITS = ((CORE_WORDS - 1) * WORD_BITS) % BITS_PER_BLOCK };
		if (usable_bits != USABLE_BITS) issue_error("NearSeq2::is_core_good - internal error 2");
		enum { FINAL_WORD_BITS = WORD_BITS - EXTRA1_PARTIAL_WORD_BITS };
		Word w2 = core[INDEX];
		w |= w2 << USABLE_BITS;
		if (USABLE_BITS + FINAL_WORD_BITS <= WORD_BITS - WORD_LEFTOVERS) {
			for (int i = 0; i < (USABLE_BITS + FINAL_WORD_BITS) / BITS_PER_BLOCK; i++) {
				is_bad |= lookup1(w);
				w >>= BITS_PER_BLOCK;
			}
		}
		else {
			for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
				is_bad |= lookup1(w);
				w >>= BITS_PER_BLOCK;
			}
			w = w2 >> (WORD_BITS - USABLE_BITS - WORD_LEFTOVERS);
			is_bad |= lookup1(w);
		}
		if (is_bad < 0) return true;
		else return false;
	}
#pragma GCC diagnostic pop
}
void PractRand::Tests::NearSeq2::core_analysis(const Word *core, int &index, int &ham) const {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
	long core_bucket = 0;
	long bucket_bit = 0;
	long h = 0;
	if (CORE_WORDS == 1) {
		Word w = core[0];
		for (int i = 0; i < BLOCKS_PER_CORE; i++) {
			analyze_block(w, core_bucket, bucket_bit++, h);
			w >>= BITS_PER_BLOCK;
		}
	}
	else if (!(WORD_BITS % BITS_PER_BLOCK)) {//blocks align to word boundaries
		int index;
		Word w = core[0];
		for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
			analyze_block(w, core_bucket, bucket_bit++, h);
			if (BITS_PER_BLOCK < WORD_BITS) w >>= BITS_PER_BLOCK;
		}
		for (index = 1; index < CORE_WORDS - 1; index++) {
			w = core[index];
			for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
				analyze_block(w, core_bucket, bucket_bit++, h);
				if (BITS_PER_BLOCK < WORD_BITS) w >>= BITS_PER_BLOCK;
			}
		}
		w = core[index];
		for (int i = 0; i < BLOCKS_PER_CORE - (CORE_WORDS - 1) * (WORD_BITS / BITS_PER_BLOCK); i++) {
			analyze_block(w, core_bucket, bucket_bit++, h);
			if (BITS_PER_BLOCK < WORD_BITS) w >>= BITS_PER_BLOCK;
		}
	}
	else {//blocks do NOT align to word boundaries
		//
		// ...this is ugly, so if possible choose parameterizations that don't hit this
		// or, failing that, paramterizations that will reject most cores in the first word
		Word w = core[0];
		for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {//first word
			analyze_block(w, core_bucket, bucket_bit++, h);
			w >>= BITS_PER_BLOCK;
		}
		int index = 1;
		enum { WORD_LEFTOVERS = WORD_BITS % BITS_PER_BLOCK };
		int usable_bits = WORD_LEFTOVERS;
		while (index < CORE_WORDS - 1) {//middle words
			Word w2 = core[index++];
			w |= w2 << usable_bits;
			for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
				analyze_block(w, core_bucket, bucket_bit++, h);
				w >>= BITS_PER_BLOCK;
			}
			usable_bits += WORD_LEFTOVERS;
			w = w2 >> (WORD_BITS - usable_bits);
			if (usable_bits >= BITS_PER_BLOCK) {
				analyze_block(w, core_bucket, bucket_bit++, h);
				w >>= BITS_PER_BLOCK;
				usable_bits -= BITS_PER_BLOCK;
			}
		}
		//last word
		enum { INDEX = CORE_WORDS - 1 };
		if (index != INDEX) issue_error("NearSeq2::analyze_core - internal error");
		enum { USABLE_BITS = ((CORE_WORDS - 1) * WORD_BITS) % BITS_PER_BLOCK };
		if (usable_bits != USABLE_BITS) issue_error("NearSeq2::analyze_core - internal error 2");
		enum { FINAL_WORD_BITS = WORD_BITS - EXTRA1_PARTIAL_WORD_BITS };
		Word w2 = core[INDEX];
		w |= w2 << USABLE_BITS;
		if (USABLE_BITS + FINAL_WORD_BITS <= WORD_BITS - WORD_LEFTOVERS) {
			for (int i = 0; i < (USABLE_BITS + FINAL_WORD_BITS) / BITS_PER_BLOCK; i++) {
				analyze_block(w, core_bucket, bucket_bit++, h);
				w >>= BITS_PER_BLOCK;
			}
		}
		else {
			for (int i = 0; i < WORD_BITS / BITS_PER_BLOCK; i++) {
				analyze_block(w, core_bucket, bucket_bit++, h);
				w >>= BITS_PER_BLOCK;
			}
			w = w2 >> (WORD_BITS - USABLE_BITS - WORD_LEFTOVERS);
			analyze_block(w, core_bucket, bucket_bit++, h);
		}
	}
	index = core_bucket;
	ham = h;
#pragma GCC diagnostic pop
}
int PractRand::Tests::NearSeq2::get_hdist_bin(int hdist) const {
	// some parameterizations have a wide variety of possible total core hamming weights
	// to keep the size of Bucket::extra_counts under control, I index them by a function of the hamming weight with more limited range
	//return hdist * HDIST_BINS / (MAX_TOTAL_HDIST + 1);
	//this could be done with a lookup table for more speed, but I think it's not called much so there's no point
	return lookup_table3[hdist];
}
void PractRand::Tests::NearSeq2::count_ones_distribution(Word bits, Uint64 *counts, int num) {
	Uint64 *end = counts + num;
	while (counts < end) {
		if (bits & 1) (*counts)++;
		bits >>= 1;
		counts++;
	}
}
void PractRand::Tests::NearSeq2::get_results(std::vector<TestResult> &results) {
	if (!blocks_tested) return;

	/*
		things to check:
		1. overall distribution of hamming weights within each bucket
			report as seperate p-values or unify in to a single p-value? lets try unifying
		2. overall distribution of extra-bits within each bucket, on a per-hamming-weight-bin basis
			report as seperate p-values or unify in to a single p-value? both by bit position and by bucket?  lets try unifying buckets and bit positions, but not hamming weights
		3. distribution between buckets, and invalid cores
	*/

	Uint64 total_cores = blocks_tested * (TestBlock::SIZE / sizeof(Word)) - (CORE_WORDS + EXTRA1_FULL_WORDS - 1);
	//if (total_cores != _total_cores) issue_error();
	Uint64 total_valid_cores = 0; for (int b = 0; b < NUM_BUCKETS; b++) for (int i = 0; i <= MAX_TOTAL_HDIST; i++) total_valid_cores += buckets[b].core_hdist[i];
	//if (total_valid_cores != _total_valid_cores) issue_error();
	Uint64 total_invalid_cores = total_cores - total_valid_cores;
	//if (total_invalid_cores != _total_invalid_cores) issue_error();

	double expected_valid_cores = total_cores * prob_of_valid_core;
	double expected_invalid_cores = total_cores * (1 - prob_of_valid_core);

	double valid_fraction = total_valid_cores / double(total_cores);

	if (false && expected_invalid_cores > 30 && expected_valid_cores > 30) {
		G_TEST cores_valid;
		cores_valid.add_category(total_invalid_cores, 1 - prob_of_valid_core);
		cores_valid.add_category(total_valid_cores, prob_of_valid_core);
		cores_valid.finalize();

		Uint64 counts[2] = { total_invalid_cores, total_valid_cores };
		double probs[2] = { 1 - prob_of_valid_core, prob_of_valid_core };
		double raw1 = cores_valid.get_result();
		double raw2 = g_test(2, probs, counts);
		double n1 = math_chisquared_to_normal(cores_valid.get_result(), cores_valid.get_DoF());
		double n2 = math_chisquared_to_normal(raw2, 1);
		double p1 = math_chisquared_to_pvalue(cores_valid.get_result(), cores_valid.get_DoF());
		double p2 = math_chisquared_to_pvalue(raw2, 1);

		std::ostringstream os;
		os << get_name() << ":valid1";
		results.push_back(TestResult(os.str(), n1, 1 - p1, TestResult::TYPE_GOOD_P, 0.01));
	}
	if (expected_invalid_cores > 30 && expected_valid_cores > 30) {//dumb, but it was catching some PRNGs ; after changes, they should now be caught by oa_hdist instead, remove this once we're sure
		double n = valid_fraction - prob_of_valid_core;
		n *= std::sqrt(total_cores / prob_of_valid_core / (1 - prob_of_valid_core));
		double p = math_normaldist_to_pvalue(-n);

		std::ostringstream os;
		os << get_name() << ":valid2";
		results.push_back(TestResult(os.str(), n, p, TestResult::TYPE_BAD_P, 0.001));
	}

	if (expected_valid_cores > 60 && total_cores > 300) {// this merges all buckets together and looks at the hdist distribution for the merger
		Uint64 hdist_counts[MAX_TOTAL_HDIST + 2];
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_counts[i] = 0;
		for (int b = 0; b < NUM_BUCKETS; b++) {
			for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_counts[i] += buckets[b].core_hdist[i];
		}
		hdist_counts[MAX_TOTAL_HDIST + 1] = total_invalid_cores;
		double hdist_probs[MAX_TOTAL_HDIST + 2];
		for (int i = 0; i <= MAX_TOTAL_HDIST; i++) hdist_probs[i] = full_probs[i];
		hdist_probs[MAX_TOTAL_HDIST + 1] = 1 - prob_of_valid_core;

		int cat = simplify_prob_table(MAX_TOTAL_HDIST + 2, total_cores / 30.0, hdist_probs, hdist_counts, true, true);

		double raw = g_test(cat, hdist_probs, hdist_counts);
		double n = math_chisquared_to_normal(raw, cat - 1);
		double p = 1 - math_chisquared_to_pvalue(raw, cat - 1);

		std::ostringstream os;
		os << get_name() << ":oa_hdist";
		results.push_back(TestResult(os.str(), n, p, TestResult::TYPE_BAD_P, 0.1));
	}

	if (total_valid_cores > NUM_BUCKETS * 5) {// here we look at distribution between buckets at each hdist (including all lower hdists), skipping some hdists to prevent excessive numbers of pvalues
		Uint64 total_count = 0;
		Uint64 counts[NUM_BUCKETS];
		Uint64 threshold = NUM_BUCKETS * 5;
		for (int i = 0; i < NUM_BUCKETS; i++) counts[i] = 0;
		for (int h = 0; h <= MAX_TOTAL_HDIST; h++) {
			for (int b = 0; b < NUM_BUCKETS; b++) {
				total_count += buckets[b].core_hdist[h];
				counts[b] += buckets[b].core_hdist[h];
			}
			if (total_count < threshold) continue;
			threshold *= 5;

			double raw = g_test_flat(NUM_BUCKETS, counts);
			double n = math_chisquared_to_normal(raw, NUM_BUCKETS - 1);
			double p = math_chisquared_to_pvalue(raw, NUM_BUCKETS - 1);
			std::ostringstream os;
			os << get_name() << ":bd@h" << h;
			results.push_back(TestResult(os.str(), n, 1 - p, TestResult::TYPE_BAD_P, 0.01));
		}
	}
	if (total_valid_cores > NUM_BUCKETS * 30) {// here we look at the hdist distribution within each bucket
		double worst_n = 0;
		int worst_b;
		double worst_p;
		double total_score = 0;
		for (int b = 0; b < NUM_BUCKETS; b++) {
			Uint64 total = 0;
			Uint64 counts[MAX_TOTAL_HDIST + 1];
			double bucket_score = 0;
			for (int h = 0; h <= MAX_TOTAL_HDIST; h++) {
				counts[h] = buckets[b].core_hdist[h];
				total += buckets[b].core_hdist[h];
				double score = buckets[b].core_hdist[h] * hdist_scores[h];
				bucket_score += score;
			}
			total_score += bucket_score * bucket_score;

			if (total < 32) continue;
			double probs[MAX_TOTAL_HDIST + 1];
			for (int h = 0; h <= MAX_TOTAL_HDIST; h++) probs[h] = full_probs[h] / prob_of_valid_core;

			int cat = simplify_prob_table(MAX_TOTAL_HDIST + 1, total / 16.0, probs, counts, true, true);
			double raw = g_test(cat, probs, counts);
			double n = math_chisquared_to_normal(raw, cat - 1);
			double p = math_chisquared_to_pvalue(raw, cat - 1);

			if (std::fabs(n) > std::fabs(worst_n)) {
				worst_n = n;
				worst_b = b;
				worst_p = p;
			}
		}
		if (total_valid_cores > NUM_BUCKETS * 30 && true) {
			double seminormalized = total_score / total_valid_cores - 0.5;
			double normalized = (seminormalized - (NUM_BUCKETS - 1)) / std::sqrt(double(NUM_BUCKETS) - 1);
			double p = math_chisquared_to_pvalue(seminormalized, NUM_BUCKETS - 1);
			std::ostringstream os;
			os << get_name() << ":WB&";
			results.push_back(TestResult(os.str(), normalized, p, TestResult::TYPE_BAD_P, 0.01));
		}
		if (total_valid_cores > NUM_BUCKETS * 50 && worst_n && true) {
			if (worst_p < 0.5) worst_p *= NUM_BUCKETS;
			else worst_p = 1 - ((1 - worst_p) * NUM_BUCKETS);
			if (worst_p < 0 || worst_p > 1) worst_p = 0.5;
			std::ostringstream os;
			os << get_name() << ":WB" << worst_b;
			results.push_back(TestResult(os.str(), worst_n, worst_p, TestResult::TYPE_RAW_NORMAL, 0.01));
		}
	}
	if (total_valid_cores > 10000 * NUM_BUCKETS) {// now we look at the "extra bits"
		struct Entry {
			double n;
			int bucket;
			int hbin;
			int ebit;
		};
		Entry low, high;
		low.n = 999999999;
		high.n = -999999999;
		Uint64 num_entries = 0;
		double chi_squared = 0;
		double chi_squared_DoF = 0;
		for (int bi = 0; bi < NUM_BUCKETS; bi++) {
			for (int hb = 0; hb < HDIST_BINS; hb++) {
				Uint64 sum_counts = 0;
				//long h = 0;
				//for (; get_hdist_bin(h) < hb; h++) ;
				//for (; h <= MAX_TOTAL_HDIST && get_hdist_bin(h) == hb; h++) sum_counts += buckets[bi].core_hdist[h];
				for (long h = 0; h <= MAX_TOTAL_HDIST; h++) if (get_hdist_bin(h) == hb) sum_counts += buckets[bi].core_hdist[h];
				if (sum_counts < 9000) continue;
				num_entries += EXTRA1_BITS;
				for (int eb = 0; eb < EXTRA1_BITS; eb++) {
					Uint64 counts[2] = { buckets[bi].extra_counts[hb][eb], sum_counts - buckets[bi].extra_counts[hb][eb] };
					chi_squared += g_test_flat(2, counts);
					chi_squared_DoF += 1;
					double raw = buckets[bi].extra_counts[hb][eb] - 0.5 * sum_counts;
					double n = raw / std::sqrt(0.25 * sum_counts);
					double p = math_normaldist_to_pvalue(n);
					if (n < low.n) {
						low.n = n;
						low.bucket = bi;
						low.hbin = hb;
						low.ebit = eb;
					}
					if (n > high.n) {
						high.n = n;
						high.bucket = bi;
						high.hbin = hb;
						high.ebit = eb;
					}
				}
			}
		}
		if (chi_squared_DoF && true) {
			std::ostringstream os2;
			os2 << get_name() << ":EB";
			double n = math_chisquared_to_normal(chi_squared, chi_squared_DoF);
			double p = math_chisquared_to_pvalue(chi_squared, chi_squared_DoF);
			results.push_back(TestResult(os2.str(), n, p, TestResult::TYPE_BAD_P, 0.1));
		}
		if (num_entries && true) {
			std::ostringstream os1;
			os1 << get_name() << ":EB:L:" << low.bucket << "_" << low.hbin << "_" << low.ebit;
			double lp = math_normaldist_to_pvalue(low.n) * num_entries;
			if (lp > 0.5) lp = 0.5;
			results.push_back(TestResult(os1.str(), low.n, lp, TestResult::TYPE_BAD_P, 0.01));
			std::ostringstream os2;
			os2 << get_name() << ":EB:H:" << high.bucket << "_" << high.hbin << "_" << high.ebit;
			double hp = math_normaldist_to_pvalue(-high.n) * num_entries;
			if (hp > 0.5) hp = 0.5;
			results.push_back(TestResult(os2.str(), high.n, hp, TestResult::TYPE_BAD_P, 0.01));
		}
	}

	return;
}

void PractRand::Tests::NearSeq2::test_blocks(TestBlock *data, int numblocks) {
	if (WORD_BITS % BITS_PER_BLOCK == 0 && EXTRA1_PARTIAL_WORD_BITS == 0) {
		long start = blocks_tested ? -(SEQUENCE_WORD_OFFSET) : SEQUENCE_WORD_OFFSET + CORE_WORDS - 1;
		long end = numblocks * (TestBlock::SIZE / sizeof(Word)) - (SEQUENCE_WORD_OFFSET);

		unsigned long invalid_words_window = 0;
		if (blocks_tested) {
			for (long pos = start - (CORE_WORDS - 1); pos < start; pos++) {
				Word w;
				if (0) ;
				else if (WORD_BITS == 8) w = data[0].as8[pos];
				else if (WORD_BITS == 16) w = data[0].as16[pos];
				else if (WORD_BITS == 32) w = data[0].as32[pos];
				else if (WORD_BITS == 64) w = data[0].as64[pos];
				else issue_error("NearS2 - what word size???");
				invalid_words_window <<= 1;
				invalid_words_window |= is_word_bad(w) ? 1 : 0;
			}
		}
		else invalid_words_window = CORE_WORDS > 1 ? 1 : 0;

		for (long pos = start; pos < end; pos++) {
			//_total_cores++;
			Word w;
			if (0) ;
			else if (WORD_BITS == 8) w = data[0].as8[pos];
			else if (WORD_BITS == 16) w = data[0].as16[pos];
			else if (WORD_BITS == 32) w = data[0].as32[pos];
			else if (WORD_BITS == 64) w = data[0].as64[pos];
			else issue_error("NearS2 - what word size???");
			invalid_words_window <<= 1;
			invalid_words_window |= is_word_bad(w) ? 1 : 0;
			invalid_words_window &= (1ull << NUM_BUCKETS_L2) - 1;
			if (invalid_words_window) {
				//_total_invalid_cores++;
				continue;
			}
			//_total_valid_cores++;
			// optimized for valid cores being relatively rare, like 1 in 40 or less - instead of keeping a running window of the match and validity, we just keep a running window of validity and reconstruct the rest when we need it

			Word *core;
			if (false) ;
			else if (WORD_BITS == 8) core = (Word*)&data[0].as8[pos];
			else if (WORD_BITS == 16) core = (Word*)&data[0].as16[pos];
			else if (WORD_BITS == 32) core = (Word*)&data[0].as32[pos];
			else if (WORD_BITS == 64) core = (Word*)&data[0].as64[pos];
			else issue_error("NearS2 - what word size???");
			core -= CORE_WORDS - 1;

			if (is_core_bad(core)) issue_error();
			int bucket_index, hdist;
			core_analysis(core, bucket_index, hdist);
			if (bucket_index < 0 || bucket_index > NUM_BUCKETS) issue_error("NearS2::text_blocks bucket_index out of range, bad analysis");

			Bucket &bucket = buckets[bucket_index];
			bucket.core_hdist[hdist]++;
			int hdist_bin = get_hdist_bin(hdist);
			Uint64 *extra_pos = &bucket.extra_counts[hdist_bin][0];
			for (int i = 1; i <= EXTRA1_FULL_WORDS / 2; i++) {
				count_ones_distribution(core[-i], extra_pos);
				extra_pos += WORD_BITS;
			}
			for (int i = 0; i < EXTRA1_FULL_WORDS / 2; i++) {
				count_ones_distribution(core[CORE_WORDS + i], extra_pos);
				extra_pos += WORD_BITS;
			}
			if (EXTRA1_PARTIAL_WORD_BITS) issue_error("NearS2 - can't have extra word bits when word is a multiple of block size");
		}
	}
	else {
		long start = blocks_tested ? -(CORE_WORDS + SEQUENCE_WORD_OFFSET - 1) : SEQUENCE_WORD_OFFSET + CORE_WORDS - 1;
		long end = numblocks * (TestBlock::SIZE / sizeof(Word)) - (CORE_WORDS + SEQUENCE_WORD_OFFSET - 1);
		for (long pos = start; pos < end; pos++) {
			Word *core;
			if (false);
			else if (WORD_BITS == 8) core = (Word*)&data[0].as8[pos];
			else if (WORD_BITS == 16) core = (Word*)&data[0].as16[pos];
			else if (WORD_BITS == 32) core = (Word*)&data[0].as32[pos];
			else if (WORD_BITS == 64) core = (Word*)&data[0].as64[pos];
			else issue_error("NearS2 - what word size???");

			//_total_cores++;
			if (is_core_bad(core)) {
				//_total_invalid_cores++;
				continue;
			}
			//_total_valid_cores++;

			int bucket_index, hdist;
			core_analysis(core, bucket_index, hdist);
			if (bucket_index < 0 || bucket_index > NUM_BUCKETS) issue_error("NearS2::test_blocks bucket_index out of range, bad analysis");

			Bucket &bucket = buckets[bucket_index];
			bucket.core_hdist[hdist]++;
			int hdist_bin = get_hdist_bin(hdist);
			Uint64 *extra_pos = &bucket.extra_counts[hdist_bin][0];
			for (int i = 1; i <= EXTRA1_FULL_WORDS / 2; i++) {
				count_ones_distribution(core[-i], extra_pos);
				extra_pos += WORD_BITS;
			}
			for (int i = 0; i < EXTRA1_FULL_WORDS / 2; i++) {
				count_ones_distribution(core[CORE_WORDS + i], extra_pos);
				extra_pos += WORD_BITS;
			}
			if (EXTRA1_PARTIAL_WORD_BITS) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
				count_ones_distribution(core[CORE_WORDS - 1] >> (WORD_BITS - EXTRA1_PARTIAL_WORD_BITS), extra_pos, EXTRA1_PARTIAL_WORD_BITS);
#pragma GCC diagnostic pop
			}
		}
	}
	blocks_tested += numblocks;
}

PractRand::Tests::NearSeq3::NearSeq3() {
	//typedef char _compiletime_assert_[(1ull << (NUM_SIZES-1)) <= (TestBlock::SIZE / sizeof(Word)) ? 1 : -1];
	if (sizeof(Word) << (NUM_SIZES-1) > TestBlock::SIZE / 2) issue_error("NearSeq3 - don't have code to safely handle these sizes yet");
	/*
		initially, max out at a half block
		if we like the results, then expand it to support maybe up to a megabyte
	*/
}
PractRand::Tests::NearSeq3::~NearSeq3() {
	if (NUM_SIZES < 4) issue_error();//optimization assumes there's at least that many
	if (THRESHOLD_TIGHTENING_FRACTION >= 9000) issue_error();//is parts per ten-thousand
}
void PractRand::Tests::NearSeq3::init(PractRand::RNGs::vRNG *known_good) {
	TestBaseclass::init(known_good);
	for (long si = 0; si < NUM_SIZES; si++) {
		PerSizeData &psd = per_size_data[si];
		long block_size = 1 << si;
		psd.size_index = si;
		long long num_bits = block_size * 8 * sizeof(Word);
		//static char initial_threshold_lookup[2] = { 5, 7 };
		//psd.threshold = (num_bits >> 1) - (initial_threshold_lookup[si & 1] << (si >> 1));
		psd.threshold = (num_bits >> 1) - std::sqrt(BITS_PER_WORD * 1.0 + si * 24.0) * std::pow(2.0, si * 0.5);
		psd.count = 0;
		psd.max_count = BUNCH_SIZE;
		psd.partial_hw = 0;
		psd.num_bits = num_bits;
		psd.even = true;

		psd.old_thresholds.resize(0);
		psd.hw_core.resize(psd.max_count);
		psd.hw_pre_.resize(psd.max_count);
		psd.hw_post.resize(psd.max_count);
		psd.raw_core.resize(psd.max_count << si);
		psd.raw_pre_.resize(psd.max_count << si);
		psd.raw_post.resize(psd.max_count << si);
	}
}
void PractRand::Tests::NearSeq3::handle_excess_sizes(long si, long hw, Word *data) {
	si += 1;
	while (si < NUM_SIZES) {
		PerSizeData &psd = per_size_data[si];
		if (psd.even) {
			psd.even = false;
			psd.partial_hw = hw;
			return;
		}
		else {
			psd.even = true;
			hw += psd.partial_hw;
			if (hw < psd.threshold) psd.handle_sample(hw, data - (1 << si));
		}
		//data -= 1 << si;
		si += 1;
	}
}
double PractRand::Tests::NearSeq3::_basic_hw_scoring_function(std::vector<double> &pdf, std::vector<double> &cdf, long hw) {
	if (hw >= pdf.size()) hw = ((pdf.size() - 1) << 1) - hw;
	//return -std::log(pdf[hw]);
	return -std::log(cdf[hw]);
	//return -std::log(cdf[hw] - pdf[hw] * 0.5);
	//return -(std::log(cdf[hw]) + std::log(1 - std::log(cdf[hw])));
}
long PractRand::Tests::NearSeq3::PerSizeData::do_scoring(ScoringData &data, bool filter) {
	long block_size = 1 << size_index;
	long num_bits = BITS_PER_WORD << size_index;
	long midpoint = num_bits >> 1;
	data.threshold = threshold;
	data.count = count;
	data.total_pairs = 0;
	data.flushed = 0;

	if (num_bits <= 16384) {//really this should be able to work for larger num_bits, maybe up to the point where memory usage becomes an issue?
		// gathering basic PDFs & CDFs
		std::vector<double> pdf, cdf;
		get_hamming_weight_chances(num_bits, pdf, cdf);
		std::vector<double> dw_pdf, dw_cdf;
		get_hamming_weight_chances(num_bits*2, dw_pdf, dw_cdf);

		// picking next threshold
		double old_rate = cdf[threshold];
		double new_rate = 0;
		long new_threshold = threshold - 1;
		while (new_threshold >= 0) {
			new_rate = cdf[new_threshold];
			if (new_rate <= (THRESHOLD_TIGHTENING_FRACTION / 1000.0) * old_rate) break;
			new_threshold -= 1;
		}
		long filter_threshold = filter ? new_threshold : -1;
		long filter_rate = filter ? cdf[filter_threshold] : 0;

		// basic callibration of basic scoring
		double mean = 0;
		double mean_sqr = 0;
		double sum_p = 0;
		for (long i = 0; i <= midpoint; i++) {
			double p = pdf[i];
			double score = p ? _basic_hw_scoring_function(pdf, cdf, i) : 0;
			if (i == midpoint) {
				mean *= 2;
				sum_p *= 2;
			}
			mean += score * p;
			sum_p += p;
		}
		for (long i = 0; i <= midpoint; i++) {
			double p = pdf[i];
			double score = p ? (_basic_hw_scoring_function(pdf, cdf, i) - mean) : 0;
			if (i == midpoint) mean_sqr *= 2;
			mean_sqr += score * score * p;
			//sum_p += p;
		}
		double scale = 1.0 / std::sqrt(mean_sqr);

		// basic callibration of double-width scoring
		double dw_mean = 0;
		double dw_mean_sqr = 0;
		double dw_sum_p = 0;
		long dw_midpoint = midpoint * 2;
		for (long i = 0; i <= dw_midpoint; i++) {
			double p = dw_pdf[i];
			double score = p ? _basic_hw_scoring_function(dw_pdf, dw_cdf, i) : 0;
			if (i == dw_midpoint) {
				dw_mean *= 2;
				dw_sum_p *= 2;
			}
			dw_mean += score * p;
			dw_sum_p += p;
		}
		for (long i = 0; i <= dw_midpoint; i++) {
			double p = dw_pdf[i];
			double score = p ? (_basic_hw_scoring_function(dw_pdf, dw_cdf, i) - dw_mean) : 0;
			if (i == dw_midpoint) dw_mean_sqr *= 2;
			dw_mean_sqr += score * score * p;
			//dw_sum_p += p;
		}
		double dw_scale = 1.0 / std::sqrt(dw_mean_sqr);

		// basic HW scoring, and counting how many samples get filtered
		data.score_hw_pre = 0;
		long filtered_count = 0;
		for (long i = 0; i < count; i++) {
			if (hw_core[i] <= filter_threshold) continue;
			long hw = hw_pre_[i];
			data.score_hw_pre += _basic_hw_scoring_function(pdf, cdf, hw) - mean;
			filtered_count++;
		}
		data.count = filtered_count;
		data.score_hw_pre *= scale;
		data.score_hw_post = 0;
		for (long i = 0; i < count; i++) {
			if (hw_core[i] <= filter_threshold) continue;
			long hw = hw_post[i];
			data.score_hw_post += _basic_hw_scoring_function(pdf, cdf, hw) - mean;
		}
		data.score_hw_post *= scale;

		// double-width HW scoring
		data.score_hw_both = 0;
		for (long i = 0; i < count; i++) {
			if (hw_core[i] <= filter_threshold) continue;
			long hw = hw_pre_[i] + hw_post[i];
			data.score_hw_both += _basic_hw_scoring_function(dw_pdf, dw_cdf, hw) - dw_mean;
		}
		data.score_hw_both *= dw_scale;

		// the more complicated case of core HW scoring
		double biased_mean = 0;
		double biased_mean_sqr = 0;
		double biased_sum_p = 0;
		//double threshold_cdf = cdf[threshold];
		double threshold_cdf_inv = 1.0 / (cdf[threshold] - filter_rate);
		for (long i = filter_threshold + 1; i <= threshold; i++) {
			double p = pdf[i] * threshold_cdf_inv;
			double score = p ? _basic_hw_scoring_function(pdf, cdf, i) : 0;
			biased_mean += score * p;
			biased_sum_p += p;
		}
		for (long i = filter_threshold + 1; i <= threshold; i++) {
			double p = pdf[i] * threshold_cdf_inv;
			double score = p ? (_basic_hw_scoring_function(pdf, cdf, i) - biased_mean) : 0;
			biased_mean_sqr += score * score * p;
			//sum_p += p;
		}
		double biased_scale = 1.0 / std::sqrt(biased_mean_sqr * count);
		data.score_hw_core = 0;
		for (long i = 0; i < count; i++) {
			long hw = hw_core[i];
			if (hw <= new_threshold) continue;
			data.score_hw_core += _basic_hw_scoring_function(pdf, cdf, hw) - biased_mean;
		}
		data.score_hw_core *= biased_scale;

		// the original hypothesis is that regions of PRNG output will have low hamming distance from other parts of the same stream, lower than they should be
		// in particular, based upon observations of chaotic PRNGs, there will be a central region of lowest hamming distance, surrounded on both sides with small regions where hamming distance is still lower than it should be but greater than the central region
		// directly testing this for the whole stream seems infeasible, so instead I look for regions of increasinly low hamming weight
		// theoretically this converges towards having samples where the core regions are all zeroes, but that functionally takes forever, more realistically it's just a gentle decrease in the hamming weights of core regions
		// any two such regions should naturally have relatively low hamming distance from each other, because they are both (relatively) close to all zeros
		// earlier tests checked that the number of low hamming weight regions was reasonable, that the distribution of their hamming weights was reasonable, and that the distributions of surrounding regions did too
		// here I test the distribution of hamming distances *between* different samples surrounding regions
		// only the closest pairs matter, but I check every pair because I can figure out the expected distribution more easily for that
		// I'd like to text the distribution of hamming distances between entire samples, but the math for that looks uglier
		long total_pairs = 0;
		double score_pd_pre = 0;
		double score_pd_post = 0;
		double score_pd_both = 0;
		// todo: add a half-both option, that looks at half-sized regions surrounding the core
		//double score_pd_core = 0;  // todo: maybe someday
		//double score_pd_all = 0;  // todo: maybe someday
		for (long i1 = 0; i1 < count; i1++) {
			if (hw_core[i1] <= filter_threshold) continue;
			const long base_pos = i1 << size_index;
			for (long i2 = 0; i2 < i1; i2++) {
				if (hw_core[i2] <= filter_threshold) continue;
				const long base_pos2 = i2 << size_index;
				long dist_pre = 0;
				long dist_post = 0;
				long dist_core = 0;
				for (long i3 = 0; i3 < block_size; i3++) dist_pre += count_ones_(raw_pre_[base_pos + i3] ^ raw_pre_[base_pos2 + i3]);
				for (long i3 = 0; i3 < block_size; i3++) dist_post += count_ones_(raw_post[base_pos + i3] ^ raw_post[base_pos2 + i3]);
				for (long i3 = 0; i3 < block_size; i3++) dist_core += count_ones_(raw_core[base_pos + i3] ^ raw_core[base_pos2 + i3]);
				long dist_both = dist_pre + dist_post;
				long dist_all = dist_both + dist_core;

				score_pd_pre += _basic_hw_scoring_function(pdf, cdf, dist_pre) - mean;
				score_pd_post += _basic_hw_scoring_function(pdf, cdf, dist_post) - mean;
				score_pd_both += _basic_hw_scoring_function(dw_pdf, dw_cdf, dist_both) - dw_mean;
				//nothing doable for core
				//nothing doable for all
				total_pairs += 1;
			}
		}
		data.score_pd_pre = score_pd_pre * scale;
		data.score_pd_post = score_pd_post * scale;
		data.score_pd_both = score_pd_both * dw_scale;
		data.total_pairs = total_pairs;

		// todo: I would also like to comparing hamming distances between the XOR results of closest pairs for each each sample.  
		//		Faintly analogous to Birthday spacing tests, though my earlier attempts at doing things analogous to Birthday Spacings yielded poor results
		//attempting now
		//pausing: more low-hanging fruit elsewhere in this function
		/*long total_pairpairs = 0;
		std::vector<Sint32> closest_distance;
		std::vector<Sint32> emperical;
		std::vector<Word> closest_pre; closest_pre.resize(block_size * count);
		std::vector<Word> closest_core;
		std::vector<Word> closest_post; closest_pre.resize(block_size * count);
		for (long i1 = 0; i1 < count; i1++) {
			long closest_index = -1;
			long closest_distance = 999999999;
			if (hw_core[i1] <= filter_threshold) continue;
			const long base_pos = i1 << size_index;
			for (long i2 = 0; i2 < i1; i2++) {
				if (hw_core[i2] <= filter_threshold) continue;
				const long base_pos2 = i2 << size_index;
				long dist_core = 0;
				long dist_pre = 0;
				long dist_post = 0;
				for (long i3 = 0; i3 < block_size; i3++) dist_core += count_ones_(raw_core[base_pos + i3] ^ raw_core[base_pos2 + i3]);
				for (long i3 = 0; i3 < block_size; i3++) dist_pre += count_ones_(raw_pre_[base_pos + i3] ^ raw_pre_[base_pos2 + i3]);
				for (long i3 = 0; i3 < block_size; i3++) dist_post += count_ones_(raw_post[base_pos + i3] ^ raw_post[base_pos2 + i3]);
				long dist_both = dist_pre + dist_post;
				long dist_all = dist_both + dist_core;

				score_pd_pre += _basic_hw_scoring_function(pdf, cdf, dist_pre) - mean;
				score_pd_post += _basic_hw_scoring_function(pdf, cdf, dist_post) - mean;
				total_pairs += 1;
			}
		}//*/

		return new_threshold;
	}
	else {
		// when calibrating, start in the middle, proceed towards the edge until cdf[i] drops below 1/million
		issue_error("PractRand::Tests::NearSeq3::PSD::archive - this path isn't written yet");
		return -1;
	}
}
void PractRand::Tests::NearSeq3::PerSizeData::tighten_threshold() {
	long block_size = 1 << size_index;
	old_thresholds.resize(old_thresholds.size() + 1);
	ScoringData &data = old_thresholds.back();
	long new_threshold = do_scoring(data, !DUMP_ALL_ON_THRESHOLD_CHANGE);
	if (new_threshold < 0) issue_error("NearSeq3::PSD::tighten_threshold() - threshold dropped below 0");

	if (!DUMP_ALL_ON_THRESHOLD_CHANGE) {
		long kept = 0;
		long src_index = 0;
		long dest_index = 0;
		for (long i = 0; i < count; i++) {
			if (hw_core[i] <= new_threshold) {
				if (kept != i) {
					for (long i2 = 0; i2 < block_size; i2++) raw_core[dest_index + i2] = raw_core[src_index + i2];
					for (long i2 = 0; i2 < block_size; i2++) raw_pre_[dest_index + i2] = raw_pre_[src_index + i2];
					for (long i2 = 0; i2 < block_size; i2++) raw_post[dest_index + i2] = raw_post[src_index + i2];
					//std::copy(&raw_core[src_index], &raw_core[src_index + block_size], &raw_core[dest_index]);
					//std::copy(&raw_pre_[src_index], &raw_pre_[src_index + block_size], &raw_pre_[dest_index]);
					//std::copy(&raw_post[src_index], &raw_post[src_index + block_size], &raw_post[dest_index]);
					hw_core[kept] = hw_core[i];
					hw_pre_[kept] = hw_pre_[i];
					hw_post[kept] = hw_post[i];
				}
				src_index += block_size;
				dest_index += block_size;
				kept++;
			}
			else {
				src_index += block_size;
			}
		}
		data.flushed = count - kept;
		if (data.count != data.flushed) issue_error("NearSeq3::tighten_threshold() - data.count != data.kept");
		count = kept;
	}
	else {
		data.flushed = count;
		if (data.count != count) issue_error("NearSeq3::tighten_threshold() - data.count != old_count");
		count = 0;
	}
	threshold = new_threshold;
}
void PractRand::Tests::NearSeq3::PerSizeData::handle_sample(long weight, Word *position) {
	if (weight > threshold) return; // this *should* only be called if weight already meets the threshold, but for optimization reasons that might not be fully up to date with recent adjustments to the threshold
	while (count >= max_count) tighten_threshold();
	long block_size = 1 << size_index;

	if (true) {
		//long memsize = sizeof(Word) * block_size;
		long base_pos = count << size_index;
		std::copy(position, position + block_size, &raw_core[base_pos]);
		std::copy(position - block_size, position, &raw_pre_[base_pos]);
		std::copy(position + block_size, position + block_size + block_size, &raw_post[base_pos]);
	}

	if (true) {
		long base_pos = count << size_index;
		long hw = 0;
		for (long i2 = 0; i2 < block_size; i2++) hw += count_ones_(raw_core[base_pos + i2]);
		hw_core[count] = weight;
		if (hw != weight) issue_error();
		hw = 0;
		for (long i2 = 0; i2 < block_size; i2++) hw += count_ones_(raw_pre_[base_pos + i2]);
		hw_pre_[count] = hw;
		hw = 0;
		for (long i2 = 0; i2 < block_size; i2++) hw += count_ones_(raw_post[base_pos + i2]);
		hw_post[count] = hw;
	}

	count += 1;
}
void PractRand::Tests::NearSeq3::test_blocks(TestBlock *data, int numblocks) {
	//Uint64 offset = blocks_tested;
	//Uint64 margin_needed = 

	enum { HALF_BLOCK = TestBlock::SIZE / sizeof(Word) / 2 };
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
	Word *p = &data[0].as64[blocks_tested ? -HALF_BLOCK : HALF_BLOCK];
	Word *end = &data[numblocks].as64[-HALF_BLOCK];
#pragma GCC diagnostic pop
	blocks_tested += numblocks;

	if (p >= end) return;

	enum {
		FAST_SIZES = 5
	};
	Uint32 threshold[FAST_SIZES];
	for (long i = 0; i < FAST_SIZES; i++) threshold[i] = per_size_data[i].threshold;

	while (p != end) {
		long countA, count0, count1, count2, count3;

#define DO_BOTTOM_LEVEL(VARIABLE,POSITION) VARIABLE = count_ones_(p[POSITION]); if (VARIABLE <= threshold[0]) per_size_data[0].handle_sample(VARIABLE, &p[POSITION]);
#define DO_MID_LEVEL(NVAR,OVAR1,OVAR2,LEVEL,POSITION) NVAR = OVAR1 + OVAR2; if (NVAR <= threshold[LEVEL]) per_size_data[LEVEL].handle_sample(NVAR, &p[POSITION]);
//#define DO_BOTTOM_LEVEL(VARIABLE,POSITION) VARIABLE = count_ones_(p[POSITION]); if (VARIABLE <= per_size_data[0].threshold) per_size_data[0].handle_sample(VARIABLE, &p[POSITION]);
//#define DO_MID_LEVEL(NVAR,OVAR1,OVAR2,LEVEL,POSITION) NVAR = OVAR1 + OVAR2; if (NVAR <= per_size_data[LEVEL].threshold) per_size_data[LEVEL].handle_sample(NVAR, &p[POSITION]);
		DO_BOTTOM_LEVEL(count0, 0);
		DO_BOTTOM_LEVEL(countA, 1);
		;;;;; DO_MID_LEVEL(count1, count0, countA, 1, 0)
		DO_BOTTOM_LEVEL(count0, 2);
		DO_BOTTOM_LEVEL(countA, 3);
		;;;;; DO_MID_LEVEL(countA, count0, countA, 1, 2)
		;;;;;;;;;; DO_MID_LEVEL(count2, count1, countA, 2, 0)
		DO_BOTTOM_LEVEL(count0, 4);
		DO_BOTTOM_LEVEL(countA, 5);
		;;;;; DO_MID_LEVEL(count1, count0, countA, 1, 4)
		DO_BOTTOM_LEVEL(count0, 6);
		DO_BOTTOM_LEVEL(countA, 7);
		;;;;; DO_MID_LEVEL(countA, count0, countA, 1, 6)
		;;;;;;;;;; DO_MID_LEVEL(countA, count1, countA, 2, 4)
		;;;;;;;;;;;;;;; DO_MID_LEVEL(count3, count2, countA, 3, 0)
		DO_BOTTOM_LEVEL(count0, 8);
		DO_BOTTOM_LEVEL(countA, 9);
		;;;;; DO_MID_LEVEL(count1, count0, countA, 1, 8)
		DO_BOTTOM_LEVEL(count0, 10);
		DO_BOTTOM_LEVEL(countA, 11);
		;;;;; DO_MID_LEVEL(countA, count0, countA, 1, 10)
		;;;;;;;;;; DO_MID_LEVEL(count2, count1, countA, 2, 8)
		DO_BOTTOM_LEVEL(count0, 12);
		DO_BOTTOM_LEVEL(countA, 13);
		;;;;; DO_MID_LEVEL(count1, count0, countA, 1, 12)
		DO_BOTTOM_LEVEL(count0, 14);
		DO_BOTTOM_LEVEL(countA, 15);
		;;;;; DO_MID_LEVEL(countA, count0, countA, 1, 14)
		;;;;;;;;;; DO_MID_LEVEL(countA, count1, countA, 2, 12)
		;;;;;;;;;;;;;;; DO_MID_LEVEL(countA, count2, countA, 3, 8)
		;;;;;;;;;;;;;;;;;;;; DO_MID_LEVEL(countA, count3, countA, 4, 0)
		handle_excess_sizes(4, countA, &p[16]);
		p += 16;
	}
}
void PractRand::Tests::NearSeq3::get_results(std::vector<TestResult> &results) {
	if (blocks_tested < 2) return;
	for (long si = 0; si < NUM_SIZES; si++) {
		PerSizeData &psd = per_size_data[si];
		std::string sis = std::to_string(si);
		//per_size_data[si].get_results(results);
/*		if (!DUMP_ALL_ON_THRESHOLD_CHANGE) {
			for (long oi = 0; oi < psd.old_thresholds.size(); oi++) {
				ScoringData &data = psd.old_thresholds[oi];
				results.push_back(TestResult("NS3[" + sis + ":hw-:" + std::to_string(oi) + "]", data.score_hw_pre, -math_normaldist_to_suspicion(data.score_hw_pre), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				results.push_back(TestResult("NS3[" + sis + ":hw0:" + std::to_string(oi) + "]", data.score_hw_core, -math_normaldist_to_suspicion(data.score_hw_core), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				results.push_back(TestResult("NS3[" + sis + ":hw+:" + std::to_string(oi) + "]", data.score_hw_post, -math_normaldist_to_suspicion(data.score_hw_post), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				results.push_back(TestResult("NS3[" + sis + ":hw*:" + std::to_string(oi) + "]", (data.score_hw_post + data.score_hw_core + data.score_hw_pre) / std::sqrt(3.0), -math_normaldist_to_suspicion((data.score_hw_post + data.score_hw_pre) / std::sqrt(2.0)), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
			}
			if (psd.count) {
				ScoringData data;
				psd.do_scoring(data, false);
				long oi = psd.old_thresholds.size();
				results.push_back(TestResult("NS3[" + sis + ":hw-:" + std::to_string(oi) + "]:" + std::to_string(data.count), data.score_hw_pre, -math_normaldist_to_suspicion(data.score_hw_pre), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				results.push_back(TestResult("NS3[" + sis + ":hw0:" + std::to_string(oi) + "]:" + std::to_string(data.count), data.score_hw_core, -math_normaldist_to_suspicion(data.score_hw_core), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				results.push_back(TestResult("NS3[" + sis + ":hw+:" + std::to_string(oi) + "]:" + std::to_string(data.count), data.score_hw_post, -math_normaldist_to_suspicion(data.score_hw_post), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				results.push_back(TestResult("NS3[" + sis + ":hw*:" + std::to_string(oi) + "]:" + std::to_string(data.count), (data.score_hw_post + data.score_hw_core + data.score_hw_pre) / std::sqrt(3.0), -math_normaldist_to_suspicion((data.score_hw_post + data.score_hw_core + data.score_hw_pre) / std::sqrt(3.0)), TestResult::TYPE_BAD_S, 1.0 / (1.0 + si)));
				// todo: :hw0 for the core, but it requires threshold-dependent normalization so is a lot more complicated
			}
		}
		else {*/
		if (true) {//all modes  should now use this, right?
			double score_hw_pre = 0;
			double score_hw_core = 0;
			double score_hw_both = 0;
			double score_hw_post = 0;
			double score_pd_pre = 0;
			double score_pd_post = 0;
			//double score_pd_core = 0;
			double score_pd_both = 0;
			double hw_variance = 0;
			double pd_variance = 0;
			const double step_mult = std::pow(THRESHOLD_TIGHTENING_FRACTION / 1000.0, 1.5);//pretty sure the exponent should be above 0.5, but not sure how far above
			const double step_mult2 = step_mult * step_mult;
			for (long oi = 0; oi <= psd.old_thresholds.size(); oi++) {
				ScoringData data;
				if (oi < psd.old_thresholds.size()) data = psd.old_thresholds[oi];
				else psd.do_scoring(data, false);

				hw_variance *= step_mult2;
				score_hw_pre *= step_mult; score_hw_core *= step_mult; score_hw_post *= step_mult; score_hw_both *= step_mult;
				pd_variance *= step_mult2;
				score_pd_pre *= step_mult; score_pd_post *= step_mult; score_pd_both *= step_mult;

				hw_variance += data.count;
				score_hw_pre += data.score_hw_pre;
				score_hw_core += data.score_hw_core;
				score_hw_post += data.score_hw_post;
				score_hw_both += data.score_hw_both;

				pd_variance += data.total_pairs;
				score_pd_pre += data.score_pd_pre;
				score_pd_post += data.score_pd_post;
				//score_pd_core += data.score_pd_core;
				score_pd_both += data.score_pd_both;
			}
			std::string history_string = "";// "{" + std::to_string(psd.old_thresholds.size()) + /*"," + std::to_string(long(hw_variance)) +*/ "}";
			if (true && !DUMP_ALL_ON_THRESHOLD_CHANGE) {
				std::vector<double> pdf, cdf;
				get_hamming_weight_chances(BITS_PER_WORD << si, pdf, cdf);
				double chances = (blocks_tested - 1) * (1024.0 / sizeof(Word)) * std::pow(0.5, si);
				double e = cdf[psd.threshold] * chances;
				double o = psd.count;
				double n = (o - e) / std::sqrt(e) / 2;
				//results.push_back(TestResult("NS3[" + sis + ":f]", n, math_normaldist_to_suspicion(n), TestResult::TYPE_BAD_S, 0.5 / (1.0 + si)));
				//results.push_back(TestResult("NS3[" + sis + ":f]" + "{" + std::to_string(psd.count) + "/" + std::to_string(long(e)) + "}", n, math_normaldist_to_suspicion(n), TestResult::TYPE_BAD_S, 0.5 / (1.0 + si)));
			}
			if (!hw_variance) continue;
			double hw_scale = 1.0 / std::sqrt(hw_variance);
			score_hw_pre *= hw_scale;
			score_hw_core *= hw_scale;
			score_hw_post *= hw_scale;
			score_hw_both *= hw_scale;
			//results.push_back(TestResult("NS3[" + sis + ":hw:pre]" + history_string, score_hw_pre, -math_normaldist_to_suspicion(score_hw_pre), TestResult::TYPE_BAD_S, 0.05 / std::sqrt(1.0 + si)));
			results.push_back(TestResult("NS3[" + sis + ":hw:core]" + history_string, score_hw_core, -math_normaldist_to_suspicion(score_hw_core), TestResult::TYPE_BAD_S, 0.2 / std::sqrt(1.0 + si)));
			//results.push_back(TestResult("NS3[" + sis + ":hw:post]" + history_string, score_hw_post, -math_normaldist_to_suspicion(score_hw_post), TestResult::TYPE_BAD_S, 0.05 / std::sqrt(1.0 + si)));
			results.push_back(TestResult("NS3[" + sis + ":hw:both]" + history_string, score_hw_both, -math_normaldist_to_suspicion(score_hw_both), TestResult::TYPE_BAD_S, 0.5 / std::sqrt(1.0 + si)));
			results.push_back(TestResult("NS3[" + sis + ":hw:all-]" + history_string, (score_hw_both + 0.5 * score_hw_core) / std::sqrt(1.25), -math_normaldist_to_suspicion(score_hw_both + 0.5 * score_hw_core) / std::sqrt(1.25), TestResult::TYPE_BAD_S, 0.2 / std::sqrt(1.0 + si)));

			if (!pd_variance) continue;
			history_string = "";// "{" + std::to_string(psd.old_thresholds.size()) + /*"," + std::to_string(long(pd_variance)) +*/ "}";
			double pd_scale = 1.0 / std::sqrt(pd_variance);
			score_pd_pre *= pd_scale;
			score_pd_post *= pd_scale;
			//score_pd_core *= pd_scale;
			score_pd_both *= pd_scale;
			//results.push_back(TestResult("NS3[" + sis + ":pd:pre]" + history_string, score_pd_pre, -math_normaldist_to_suspicion(score_pd_pre), TestResult::TYPE_BAD_S, 0.1 / std::sqrt(1.0 + si)));
			//results.push_back(TestResult("NS3[" + sis + ":pd:post]" + history_string, score_pd_post, -math_normaldist_to_suspicion(score_pd_post), TestResult::TYPE_BAD_S, 0.1 / std::sqrt(1.0 + si)));
			results.push_back(TestResult("NS3[" + sis + ":pd:both]" + history_string, score_pd_both, -math_normaldist_to_suspicion(score_pd_both), TestResult::TYPE_BAD_S, 1.0 / std::sqrt(1.0 + si)));
		}
	}
}
std::string PractRand::Tests::NearSeq3::get_name() const {
	return "NearSeq3";
}
void PractRand::Tests::NearSeq3::deinit() {
	TestBaseclass::deinit();
}





void PractRand::Tests::LHWSCN::sample_to_counts(unsigned long level, Word *sample, long *counts) {
	const unsigned long N = 1 << level;
	Uint8 __assert_wordsize__[sizeof(Word) == 8*sizeof(Uint8) ? 1 : -1]; //this code assumes that Word is 64 bits atm
	const Uint64 xc[6] = { 0x5555555555555555ull, 0x3333333333333333ull, 0x0F0F0F0F0F0F0F0Full, 0x00FF00FF00FF00FFull, 0x0000FFFF0000FFFFull, 0x00000000FFFFFFFFull };
	Uint8 hw[1 << MAX_LEVELS];
	long thw = 0;
	for (unsigned long i = 0; i < N; i++) thw += hw[i] = count_ones64(sample[i]);
	for (unsigned long i = 0; i < 6 + level; i++) counts[i] = 0;
	for (unsigned long i = 0; i < 6; i++) for (unsigned long n = 0; n < N; n++) counts[i] += count_ones64(sample[n] ^ xc[i]);
	for (unsigned long i = 0; i < level; i++) {
		long dir[2] = { 0, 0 };
		for (unsigned long n = 0; n < N; n++) dir[(n >> i) & 1] += hw[n];
		counts[i + 6] += dir[0] - dir[1];
	}
	counts[level + 6] = thw;
	return;
}
PractRand::Tests::LHWSCN::Classification PractRand::Tests::LHWSCN::counts_to_classification(unsigned long level, long counts[]){
	Classification c = 0;
	level += 6;
	for (int i = 0; i < level; i++) c |= ((counts[i] >> (sizeof(*counts) * 8 - 1)) & 1) << i;
	return c;
}
void PractRand::Tests::LHWSCN::handle_sample(unsigned long level, Word *sample) {
	LevelData *ld = &level_data[level];
	long counts[7 + MAX_LEVELS];
	sample_to_counts(level, sample, counts);
	if (counts[6 + level] < ld->threshold_low) {
		;
	}
	long N = level + 7;
	if (!ld->parity) {
		ld->parity = true;
		for (int i = 0; i < N; i++) ld->carry_data[i] = counts[i];
	}
	else {
		ld->parity = false;
	}
}
PractRand::Tests::LHWSCN::LHWSCN(unsigned long _minimum_level) : minimum_level(_minimum_level) {
	if (sizeof(Word) != 8) issue_error("LHWSCN constructor - must operate on 64 bit words");
	if (minimum_level > MAX_MINIMUM_LEVEL) issue_error("LHWSCN constructor - minimum_level too high");
	const long word_L2 = 6;//fixed
	if (minimum_level > (TestBlock::SIZE_L2 + 3 - word_L2)) issue_error("LHWSCN constructor - minimum level too large");
	for (int i = minimum_level; i < MAX_LEVELS; i++) {
		level_data[i].storage = new Word[SAMPLES_KEPT << i];
		long bits_L2 = word_L2 + i;
		if (bits_L2 < 6) issue_error("LHWSCN constructor - level too low");
	}
}
PractRand::Tests::LHWSCN::~LHWSCN() {
	for (int i = minimum_level; i < MAX_LEVELS; i++) {
		delete[] level_data[i].storage;
	}
}
void PractRand::Tests::LHWSCN::init(PractRand::RNGs::vRNG *known_good) {
	autofail = false;
	long word_L2 = 6;//fixed
	for (int i = 0; i < MAX_LEVELS; i++) {
		level_data[i].num_samples = 0;
		level_data[i].parity = false;

		long bits_L2 = word_L2 + i;
		long bits = 1 << bits_L2;
		long bits_var;
		if (bits_L2 < 6) issue_error("LHWSCN::init - level too low");
		else bits_var = ((bits_L2 & 1) ? 17 : 12) << ((i - 6) >> 1);
		//		6		7		8		9		10		11		12		13		14		15		16		17		18		19		20
		//		12		17		24		34		48		68		96		136		192		272		384		544		768		1088	1536
		//
		//			1x: e-3			2x: e-9			3x: e-19		4x: e-34
		// for 1 kilobyte samples that should be 136
		long low = bits / 2 - 3 * bits_var - 1; if (low < 0) low = 0;
		long high = bits / 2 - 1;
		level_data[i].threshold_high = high;
		level_data[i].threshold_low = low;
		level_data[i].count_at_weight.resize(high + 1 - low, 0);
		level_data[i].outliers.clear();
	}
}
void PractRand::Tests::LHWSCN::test_blocks(TestBlock *data, int numblocks) {
	if (!numblocks) return;
	Word *cur = &data->as64[0];
	Word *end = cur + numblocks * static_cast<long long>(TestBlock::SIZE / sizeof(Word));
	long step = 1 << minimum_level;
	for (; cur < end; cur += step) {
		handle_sample(minimum_level, cur);

	}
	blocks_tested += numblocks;
}

















void PractRand::Tests::Coup16::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	for (int i = 0; i < S; i++) flags[i] = 0;
	blocks_tested = 0;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::Coup16::get_name() const {
	return "Coup16";
}
void PractRand::Tests::Coup16::test_blocks(TestBlock *data, int numblocks) {
	while (numblocks) {
		int blocks_to_use = 128 - (blocks_tested & 127);
		if (blocks_to_use > numblocks) blocks_to_use = numblocks;
		int max = blocks_to_use * (TestBlock::SIZE / sizeof(Uint16));
		for (int i = 0; i < max; i++) {
			Uint16 word = data[0].as16[i];
			flags[word >> 5] |= 1 << (word & 31);
		}
		blocks_tested += blocks_to_use;
		data += blocks_to_use;
		numblocks -= blocks_to_use;

		if (!(blocks_tested & 127)) {
			int sum = 0;
			for (int i = 0; i < S; i++) {
				sum += count_ones32(flags[i]);
				flags[i] = 0;
			}
			counts.increment(sum - 1);
		}
	}
}
void PractRand::Tests::Coup16::get_results(std::vector<TestResult> &results) {
	//static Uint64 print_at = 128 << 13;
	//if (blocks_tested < print_at) return;
	//print_at <<= 1;
	const Uint64 *count = counts.get_array();
	const double expected_mean = 41426.652943388356 - 1 + 0.185;// plus or minus about 0.001?... actually the value I got empirically was 0.184979, but it 0.185 was so close and so much prettier
	const double expected_deviation = 79.81665;// plus or minus about 0.001? - these valuse were obtained from a 1280 TB test run
	double eebar = 0.001;

	double total_error = 0;
	double weighted_error = 0;
	double inv_dev = 1.0 / expected_deviation;
	double total = blocks_tested / 128;
	if (total < 1) return;
	for (int i = 0; i < 65536; i++) {
		if (i > 40000 && i < 43000) continue;
		if (count[i]) {
			results.push_back(TestResult(get_name() + ":!", 1, 1, TestResult::TYPE_PASSFAIL, 0.0000001));
			return;
		}
	}

	std::vector<double> probs; probs.resize(3000);
	for (int i = 0; i < 3000; i++) {
		int ei = i + 40000;
		double c = count[ei];
		double f = c / total;
		double p = math_normaldist_pdf((ei - expected_mean) * inv_dev) * inv_dev;
		probs[i] = p;
		double error = std::fabs(f - p);
		double weight = -std::log(p);
		total_error += error;
		weighted_error -= error * std::log(p);
	}
	//results.push_back(TestResult(get_name() + ":A", weighted_error, weighted_error, TestResult::TYPE_RAW_NORMAL, 0.01));
	if (total > 1) {
		double norm = rarity_test(3000, &probs[0], &counts[40000], false, true);
		results.push_back(TestResult(get_name() + ":A", norm, norm, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
	if (total > 100) {
		Uint64 counts2[3000];
		for (int i = 0; i < 3000; i++) counts2[i] = counts[i+40000];
		int cat = simplify_prob_table(3000, total * 1.0, &probs[0], &counts2[0], true, false);
		double chisqr = g_test(cat, &probs[0], &counts2[0]);
		double norm = math_chisquared_to_normal(chisqr, cat - 1);
		results.push_back(TestResult(get_name() + ":B", norm, norm, TestResult::TYPE_RAW_NORMAL, 0.01));
	}


	/*
	std::vector<double> probs; probs.resize(65536);
	double total = blocks_tested / 128;
	double sum = 0;
	double sum_sqr = 0;
	for (int i = 41426 - 750; i < 41426 + 750; i++) {
		double c = count[i];
		double f = c / total;
		sum += i * c;
		sum_sqr += i * i * c;
	}
	double mean = sum / total;//should be ((65536 * (1 - exp(-1))) - 1) = 41425.652943388356, I think
	double mean_sqr = sum_sqr / total;
	double dev = std::sqrt(mean_sqr - mean * mean);//should be 79.8, I guess?
	double total_error = 0;
	double sum_p = 0;
	for (int i = 41426 - 750; i < 41426 + 750; i++) {
		double c = count[i];
		double f = c / total;
		double p = math_normaldist_pdf((i - avg) / dev) / dev;
		if (c) std::printf("%5d: %.15f    %.15f\n", i, f, p);
		total_error += std::fabs(f - p);
		sum_p += p;
	}
	std::printf("total: %.0f    mean: %.7f   dev: %.9f   err: %.9f\n", total, avg, dev, total_error);
	std::printf("sum_p: %.9f\n", sum_p);
	std::printf("blocks tested = %.0f GB, print_at = %.0f GB\n", double(blocks_tested >> 20), double(print_at >> 20));
	std::printf("\n");//*/
}

void PractRand::Tests::Coup32::init(PractRand::RNGs::vRNG *known_good) {
	blocks_tested = 0;
	set_start_position = 0;
	filled_words = 0;
	std::memset(bitmap, 0, sizeof(bitmap));
	num_complete_sets = 0;
	sum_log2_set_size = 0;
	sum_sqr_log2_set_size = 0;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::Coup32::get_name() const {
	return "Coup32";
}
void PractRand::Tests::Coup32::handle_set_completion(Uint64 position) {
	Uint64 distance = position - set_start_position;
	set_start_position = position;
	num_complete_sets += 1;
	filled_words = 0;
	double N = std::pow(2.0, COUPON_BITS);
	double expected_distance = N * std::log(N);
	// todo: needs more scoring here
}
void PractRand::Tests::Coup32::test_blocks(TestBlock *data, int numblocks) {
	while (numblocks) {
		for (long i = 0; i < TestBlock::SIZE / sizeof(Uint32); i++) {
			Uint32 c = data->as32[i];
			std::size_t index = c >> BITMAP_WORD_BITS_L2;
			BitmapWord raised_bit = BitmapWord(1) << (c & ((BitmapWord(1) << BITMAP_WORD_BITS_L2) - 1));
			BitmapWord old_word = bitmap[index];
			BitmapWord new_word = old_word | raised_bit;
			bitmap[index] = new_word;
			if (new_word == ~BitmapWord(0) && new_word != old_word) {
				filled_words += 1;
				if (filled_words == 1ull << (32 - BITMAP_WORD_BITS_L2)) {
					Uint64 pos = blocks_tested * (TestBlock::SIZE / sizeof(Uint32)) + i;
					handle_set_completion(pos);
				}
			}
		}
		numblocks--;
		data++;
		blocks_tested++;
	}
}
void PractRand::Tests::Coup32::get_results(std::vector<TestResult> &results) {
	results.push_back(TestResult("Coup32:mean", double(set_start_position) / num_complete_sets, 0, TestResult::TYPE_UNKNOWN, 1.0));
}


void PractRand::Tests::DistFreq4::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next = blocks_per - 1;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::DistFreq4::get_name() const {
	std::ostringstream buf;
	buf << "DF4(/" << blocks_per << ")";
	return buf.str();
}
void PractRand::Tests::DistFreq4::get_results(std::vector<TestResult> &results) {
	enum { TSIZE = 1 << (SIZE1 + SIZE2) };
	int num_sweeps = blocks_tested / blocks_per;
	if (num_sweeps < TSIZE * 100) return;
	int DoF = ((1 << (SIZE1 + SIZE2)) - 1) << (POSITIONS1_L2 + POSITIONS2_L2);
	const Uint64 *counts_ = counts.get_array();
	double chisqr = g_test_flat(1 << TOTAL_INDEX_BITS, counts_);
	double norm = math_chisquared_to_normal(chisqr, DoF);
	results.push_back(TestResult(get_name() + ":all", norm, norm, TestResult::TYPE_RAW_NORMAL, 0.01));
	Uint64 counts2[TSIZE];
	double highest1 = -9999999;
	int highest_pos1 = 0;
	double highest2 = -9999999;
	int highest_pos2 = 0;
	for (int pos1 = 0; pos1 < 1 << POSITIONS1_L2; pos1++) {
		int base_index1 = pos1 << (TOTAL_INDEX_BITS - POSITIONS1_L2);
		for (int pos2 = 0; pos2 < 1 << POSITIONS2_L2; pos2++) {
			int base_index2 = base_index1 + (pos2 << SIZE2);
			int new_index = 0;
			for (int value1 = 0; value1 < (1 << SIZE1); value1++) {
				int old_index = base_index2 + (value1 << (TOTAL_INDEX_BITS - POSITIONS1_L2 - SIZE1));
				for (int value2 = 0; value2 < (1 << SIZE2); value2++) {
					counts2[new_index++] = counts_[old_index++];
				}
			}
			double chisqr2 = g_test_flat(TSIZE, &counts2[0]);
			double norm2 = math_chisquared_to_normal(chisqr2, TSIZE - 1);
			if (norm2 > highest1) {
				highest2 = highest1;
				highest_pos2 = highest_pos1;
				highest1 = norm2;
				highest_pos1 = (pos1 << POSITIONS2_L2) + pos2;
			}
			else if (norm2 > highest2) {
				highest2 = norm2;
				highest_pos2 = (pos1 << POSITIONS2_L2) + pos2;
			}
		}
	}
	if (true) {
		std::ostringstream buf;
		buf << get_name() << ":h1(" << std::hex << highest_pos1 << ")";
		results.push_back(TestResult(buf.str(), highest1, highest1, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
	if (true) {
		std::ostringstream buf;
		buf << get_name() << ":h2(" << std::hex << highest_pos2 << ")";
		results.push_back(TestResult(buf.str(), highest2, highest2, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
}
void PractRand::Tests::DistFreq4::test_blocks(TestBlock *data, int numblocks) {
	blocks_tested += numblocks;
	while (numblocks > blocks_till_next) {
		data += blocks_till_next;
		numblocks -= blocks_till_next;
		for (int pos1 = 0; pos1 < 1 << POSITIONS1_L2; pos1++) {
			int bits_used = pos1 * ALIGNMENT1;
			Uint32 first = data[0].as32[bits_used >> 5] >> (bits_used & 31);
			Uint32 base_index = (pos1 << (TOTAL_INDEX_BITS - POSITIONS1_L2)) + ((first & ((1 << SIZE1) - 1)) << (TOTAL_INDEX_BITS - POSITIONS1_L2 - SIZE1));
			bits_used += SIZE1;
			if (ALIGNMENT1 % ALIGNMENT2 || SIZE1 % ALIGNMENT2) { bits_used = bits_used + ALIGNMENT2 - 1; bits_used &= 65535 ^ (ALIGNMENT2 - 1); }
			enum { ALIGNMENTS_PER_WORD = 32 / ALIGNMENT2 };
			int end_index = base_index + (1 << (TOTAL_INDEX_BITS - POSITIONS1_L2 - SIZE1));
			if (bits_used & 31) {//partial word
				Uint32 second = data[0].as32[bits_used >> 5];
				second >>= bits_used & 31;
				while (bits_used & 31) {
					counts.increment(base_index + (second & ((1 << SIZE2) - 1)));
					base_index += 1 << SIZE2;
					second >>= ALIGNMENT2;
					bits_used += ALIGNMENT2;
				}
			}
			int words_used = bits_used >> 5;
			while (base_index <= end_index - (ALIGNMENTS_PER_WORD << SIZE2)) {
				Uint32 second = data[0].as32[words_used];
				for (int i = 0; i < ALIGNMENTS_PER_WORD; i++) {
					counts.increment(base_index + (second & ((1 << SIZE2) - 1)));
					second >>= ALIGNMENT2;
					base_index += 1 << SIZE2;
				}
				words_used++;
			}
			if (base_index < end_index) {//partial word
				Uint32 second = data[0].as32[bits_used >> 5];
				while (base_index < end_index) {
					counts.increment(base_index + (second & ((1 << SIZE2) - 1)));
					second >>= ALIGNMENT2;
					base_index += 1 << SIZE2;
				}
			}
		}
		data += 1;
		numblocks -= 1;
		blocks_till_next = blocks_per - 1;
	}
	blocks_till_next -= numblocks;
}







void PractRand::Tests::TripleFreq::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next_pass = blocks_per_pass - 1;
	passes_till_next_region = PASSES_PER_REGION;
	regions_tested = 0;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::TripleFreq::get_name() const {
	std::ostringstream buf;
	buf << "TFreq(" << passes_at_once;
	if (blocks_per_pass != 1) buf << "/" << blocks_per_pass;
	buf << ")";
	return buf.str();
}
void PractRand::Tests::TripleFreq::get_results(std::vector<TestResult> &results) {
	const Uint64 *counts_ = counts.get_array();
	enum { 
		SECTOR_SIZE = 1 << SIZE3,
		PATTERN_SIZE = 1 << (SIZE3 + SIZE2 + SIZE1),
		REGION_SIZE = 1 << REGION_INDEX_BITS, 
		TOTAL_SIZE = 1 << TOTAL_INDEX_BITS
	};
	/*
					pass through to g_test		collapse		repeat across				notes
		subtest 1	SIZE3						--				SIZE1, SIZE2, POS2, POS3	minimum count on a per-sector basis
		subtest 2	SIZE1,SIZE2,SIZE3			--				POS2,POS3					--
		subtest 3	SIZE1,SIZE2,SIZE3,POS2,POS3	--				--							requires large dataset?

		subtest 1 can be done in place easily
		subtests 2 & 3 require reordering the bits.  It could be done in the desired order at test time, but that could reduce memory locality.
			maybe it's worthwhile anyway?
				NO.  tested it, big impact on performance
			so... we just reorder them here... slowing down results reporting and using twice as much memory during results reporting
	*/
	if (!regions_tested) return;
	std::vector<Uint64> counts2;//re-ordered for better testability
	counts2.resize(TOTAL_SIZE);
	Uint64 *counts2_ = &counts2[0];
	int num_regions = regions_tested;
	if (num_regions > NUMBER_OF_REGIONS) num_regions = NUMBER_OF_REGIONS;
	int max = num_regions * REGION_SIZE;
	for (int sector_base = 0; sector_base < max; sector_base += 1 << SIZE3) {
		int new_sector_base = (sector_base & (REGION_SIZE - (1 << (SIZE3 + POSITIONS3_L2)))) >> POSITIONS3_L2;//window 1 & 2
		new_sector_base |= sector_base & (TOTAL_SIZE - REGION_SIZE);//position 2
		new_sector_base |= (sector_base & (((1 << POSITIONS3_L2) - 1) << SIZE3)) << (REGION_INDEX_BITS - POSITIONS3_L2 - SIZE3);//position 3
		//for (int i = 0; i < (1 << SIZE3); i++) counts2_[new_sector_base + i] = counts_[sector_base];//window 3
		std::memcpy(counts2_ + new_sector_base, counts_ + sector_base, sizeof(Uint64) << SIZE3);
	}
	double worst_sector_n = 0;
	int worst_sector_index = -1;
	double worst_pattern_n = 0;
	int worst_pattern_index = -1;
	for (int x = 0; x < NUMBER_OF_REGIONS; x++) {
		if (x >= regions_tested) break;
		//subtest 1
		for (int y = 0; y < REGION_SIZE >> SIZE3; y++) {
			Uint64 sum = 0;
			Uint32 sector_index = y + (x << (REGION_INDEX_BITS - SIZE3));
			const Uint64 *sector = counts2_ + (sector_index << SIZE3);
			for (int z = 0; z < (1 << SIZE3); z++) sum += sector[z];
			if (sum < (10 << SIZE3)) continue;
			double chisquared = g_test_flat(1 << SIZE3, sector);
			double n = math_chisquared_to_normal(chisquared, (1 << SIZE3) - 1);
			if (n > worst_sector_n) {
				worst_sector_n = n;
				worst_sector_index = sector_index;
			}
		}
		if (regions_tested < NUMBER_OF_REGIONS * 16) continue;
		//subtest 2
		for (int y = 0; y < (1 << POSITIONS3_L2); y++) {
			int pattern_index = y + (x << POSITIONS3_L2);
			double chisquared = g_test_flat(PATTERN_SIZE, counts2_ + (pattern_index * PATTERN_SIZE));
			double n = math_chisquared_to_normal(chisquared, PATTERN_SIZE - 1);
			if (n > worst_pattern_n) {
				worst_pattern_n = n;
				worst_pattern_index = pattern_index;
			}
		}
	}
	if (worst_sector_index != -1) {
		std::ostringstream buf;
		buf << get_name() << ":sec(" << std::hex << worst_sector_index << ")";
		results.push_back(TestResult(buf.str(), worst_sector_n, worst_sector_n, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
	if (worst_pattern_index != -1) {
		std::ostringstream buf;
		buf << get_name() << ":pat(" << std::hex << worst_pattern_index << ")";
		results.push_back(TestResult(buf.str(), worst_pattern_n, worst_pattern_n, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
	if (regions_tested >= NUMBER_OF_REGIONS * 64) {
		//subtest 3
		double chisquared = g_test_flat(TOTAL_SIZE, counts2_);
		//double n = math_chisquared_to_normal(chisquared, TOTAL_SIZE - 1);
		double n = math_chisquared_to_normal(chisquared, TOTAL_SIZE - NUMBER_OF_REGIONS);
		std::ostringstream buf;
		buf << get_name() << ":all";
		results.push_back(TestResult(buf.str(), n, n, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
}
static Uint64 read_64_misaligned(const Uint64 *source, int bit_pos) {
	int i = bit_pos >> 6;
	Uint64 rv = source[i];
	int b = bit_pos & 63;
	if (b) {
		rv >>= b;
		rv |= source[i + 1] << (64 - b);
	}
	return rv;
}
void PractRand::Tests::TripleFreq::test_blocks(TestBlock *data, int numblocks) {
	blocks_tested += numblocks;
	while (numblocks > blocks_till_next_pass) {
		data += blocks_till_next_pass;
		numblocks -= blocks_till_next_pass;
		for (int pos1 = 0; pos1 < passes_at_once; pos1++) {
			Uint64 window1 = read_64_misaligned(data[0].as64, pos1 * BASE_ALIGNMENT);
			
			int p2i = regions_tested & ((1 << POSITIONS2_L2) - 1);
			Uint64 window2 = read_64_misaligned(data[0].as64, pos1 * BASE_ALIGNMENT + SIZE1 + p2i * WINDOW_ALIGNMENT);
			
			Uint32 base_index = p2i;
			base_index <<= SIZE1; base_index |= window1 & ((1 << SIZE1) - 1);
			base_index <<= SIZE2; base_index |= window2 & ((1 << SIZE2) - 1);
			base_index <<= SIZE3 + POSITIONS3_L2;
			int position = pos1 * BASE_ALIGNMENT + SIZE1 + p2i * WINDOW_ALIGNMENT + SIZE2;
			enum { N = (64 + WINDOW_ALIGNMENT - SIZE3) / WINDOW_ALIGNMENT };
			for (int p3in = 0; p3in < (1 << POSITIONS3_L2) / N; p3in++) {
				Uint64 window3 = read_64_misaligned(data[0].as64, position);
				for (int x = 0; x < N; x++) {
					counts.increment(base_index + (window3 & ((1 << SIZE3) - 1)));
					base_index += 1 << SIZE3;
					window3 >>= WINDOW_ALIGNMENT;
				}
				position += WINDOW_ALIGNMENT * N;
			}
			if ((1 << POSITIONS3_L2) % N) {//the left-overs
				Uint64 window3 = read_64_misaligned(data[0].as64, position);
				for (int x = 0; x < ((1 << POSITIONS3_L2) % N); x++) {
					counts.increment(base_index + (window3 & ((1 << SIZE3) - 1)));
					base_index += 1 << SIZE3;
					window3 >>= WINDOW_ALIGNMENT;
				}
			}
			if (--passes_till_next_region <= 0) {
				regions_tested++;
				passes_till_next_region = PASSES_PER_REGION;
			}
		}

		data += 1;
		numblocks -= 1;
		blocks_till_next_pass = blocks_per_pass - 1;
	}
	blocks_till_next_pass -= numblocks;
}


void PractRand::Tests::TripleMirrorFreq::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next_pass = blocks_per_pass - 1;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::TripleMirrorFreq::get_name() const {
	std::ostringstream buf;
	buf << "TMF(" << passes_at_once;
	if (blocks_per_pass != 1) buf << "/" << blocks_per_pass;
	buf << ")";
	return buf.str();
}
void PractRand::Tests::TripleMirrorFreq::get_results(std::vector<TestResult> &results) {
	const Uint64 *counts_ = counts.get_array();
	int repeat_blocks = get_blocks_to_repeat();
	if (blocks_tested < repeat_blocks) return;
	Sint64 passes = ((blocks_tested - repeat_blocks) / blocks_per_pass) * passes_at_once;
	double E = passes * std::pow(0.5, SIZE1 + SIZE2 + SIZE3);
	if (E < 10) return;
	int worst_position = -1;
	double worst_n = 0, overall_cs1 = 0, overall_n2 = 0;
	int cat = 1 << (SIZE1 + SIZE2 + SIZE3);
	for (int pos = 0; pos < (1 << POSITIONS_L2); pos++) {
		double chisquared = g_test_flat(cat, &counts_[pos << (SIZE1 + SIZE2 + SIZE3)]);
		double n = math_chisquared_to_normal(chisquared, cat);
		if (std::fabs(n) >= std::fabs(worst_n)) {
			worst_n = n;
			worst_position = pos;
		}
		overall_cs1 += chisquared;
		overall_n2 += n * n;
	}
	overall_n2 = (std::sqrt(overall_n2 / (1 << POSITIONS_L2)) - 1.46) * std::pow(2.0, 0.5 * POSITIONS_L2);
	std::ostringstream buf;
	buf << get_name() << ":w(" << worst_position << ")";
	results.push_back(TestResult(buf.str(), worst_n, worst_n, TestResult::TYPE_RAW_NORMAL, 0.01));
	if (E < 80) return;
	// none of these are working, and I don't know why
	double overall_n1 = math_chisquared_to_normal(overall_cs1, (cat - 1) << POSITIONS_L2);
	buf.str("");
	buf << get_name() << ":all1";
	results.push_back(TestResult(buf.str(), overall_n1, overall_n1, TestResult::TYPE_RAW_NORMAL, 0.01));
	buf.str("");
	buf << get_name() << ":all2";
	results.push_back(TestResult(buf.str(), overall_n2, overall_n2, TestResult::TYPE_RAW_NORMAL, 0.01));
}
int PractRand::Tests::TripleMirrorFreq::get_blocks_to_repeat() const {
	enum { POSITIONS = 1 << POSITIONS_L2 };
	//int bytes_needed = ((1 << POSITIONS_L2) + passes_at_once) * SAMPLE_ALIGN + TestBlock::SIZE * 2 * BLOCK_STEP;
	int bytes_needed = TestBlock::SIZE * 2 * BLOCK_STEP + (POSITIONS << POSITION_ALIGN_L2) + (passes_at_once << BASE_ALIGN_L2);// if BLOCK_STEP and 0-based positions are used
	return (bytes_needed + TestBlock::SIZE - 1) >> TestBlock::SIZE_L2;
}
void PractRand::Tests::TripleMirrorFreq::test_blocks(TestBlock *data, int numblocks) {
	while (blocks_tested < get_blocks_to_repeat()) {
		if (!numblocks) return;
		data += 1;
		numblocks -= 1;
		blocks_tested += 1;
	}
	blocks_tested += numblocks;
	while (numblocks > blocks_till_next_pass) {
		data += blocks_till_next_pass;
		numblocks -= blocks_till_next_pass;
		enum { BASE_ALIGN = 1 << BASE_ALIGN_L2, POSITION_ALIGN = 1 << POSITION_ALIGN_L2, POSITIONS = 1 << POSITIONS_L2 };
		for (long pos_code = 0; pos_code < POSITIONS; pos_code ++) {//0-based if BLOCK_STEP is used
			unsigned long base_index = pos_code << (SIZE1 + SIZE2 + SIZE3);//0-based if BLOCK_STEP is used
			long ofs = BLOCK_STEP * TestBlock::SIZE - pos_code * POSITION_ALIGN;
			for (long base_pos = -passes_at_once * BASE_ALIGN; base_pos < 0; base_pos += BASE_ALIGN) {
				unsigned long index = base_index;
				//Uint8 a = data[-BLOCK_STEP].as8[base_pos + 0], b = data[0].as8[base_pos + pos], c = data[-2 * BLOCK_STEP].as8[base_pos - pos];
				Uint8 a = data[0].as8[base_pos - 0], b = data[0].as8[base_pos - ofs], c = data[0].as8[base_pos - ofs - ofs];// if BLOCK_STEP is used
				//b -= a; a -= c; // nope.  these hurt more than they help, typically
				index |= ((unsigned long)(a & ((1 << SIZE1) - 1))) << (SIZE2 + SIZE3);
				index |= ((unsigned long)(b & ((1 << SIZE2) - 1))) << (SIZE3);
				index |= ((unsigned long)(c & ((1 << SIZE3) - 1))) << (0);
				counts.increment(index);
			}
		}

		data += 1;
		numblocks -= 1;
		blocks_till_next_pass = blocks_per_pass - 1;
	}
	blocks_till_next_pass -= numblocks;
}





void PractRand::Tests::TripleMirrorFreqN::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next_pass = 0;// blocks_per_pass - 1;
	for (int i = 0; i < MAX_LEVELS; i++) level_state[i] = 0;
	for (int i = 0; i < MAX_LEVELS; i++) level_polarity[i] = 0;
	if ((ALIGN << POSITIONS_L2) > 64) issue_error("TripleMirrorFreqN::init - bad configuration");
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::TripleMirrorFreqN::get_name() const {
	std::ostringstream buf;
	buf << "TMFn(";
	if (blocks_per_pass != 1) buf << minimum_level;
	buf << ")";
	return buf.str();
}
void PractRand::Tests::TripleMirrorFreqN::get_results(std::vector<TestResult> &results) {
	const Uint64 *counts_ = counts.get_array();
	for (int level = minimum_level; level < MAX_LEVELS; level++) {
		if ((blocks_tested >> level) < 12 << (SIZE1 + SIZE2 + SIZE3)) return;

		double all_cs = g_test_flat(1 << TOTAL_INDEX_BITS, &counts_[level << TOTAL_INDEX_BITS]);
		double all_n = math_chisquared_to_normal(all_cs, (1 << TOTAL_INDEX_BITS) - (1 << POSITIONS_L2));

		std::ostringstream buf;
		std::string level_name;
		buf << "TMFn(" << minimum_level << "+" << (level - minimum_level) << ")";
		level_name = buf.str();
		buf.str("");
		buf << level_name << ":wl";
		results.push_back(TestResult(buf.str(), all_n, all_n, TestResult::TYPE_RAW_NORMAL, 0.01));
	}
/*	Sint64 passes = ((blocks_tested - repeat_blocks) / blocks_per_pass) * passes_at_once;
	double E = passes * std::pow(0.5, SIZE1 + SIZE2 + SIZE3);
	if (E < 10) return;
	int worst_position = -1;
	double worst_n = 0, overall_cs1 = 0, overall_n2 = 0;
	int cat = 1 << (SIZE1 + SIZE2 + SIZE3);
	for (int pos = 0; pos < (1 << POSITIONS_L2); pos++) {
		double chisquared = g_test_flat(cat, &counts_[pos << (SIZE1 + SIZE2 + SIZE3)]);
		double n = math_chisquared_to_normal(chisquared, cat);
		if (std::fabs(n) >= std::fabs(worst_n)) {
			worst_n = n;
			worst_position = pos;
		}
		overall_cs1 += chisquared;
		overall_n2 += n * n;
	}
	overall_n2 = (std::sqrt(overall_n2 / (1 << POSITIONS_L2)) - 1.46) * std::pow(2.0, 0.5 * POSITIONS_L2);
	std::ostringstream buf;
	buf << get_name() << ":w(" << worst_position << ")";
	results.push_back(TestResult(buf.str(), worst_n, worst_n, TestResult::TYPE_RAW_NORMAL, 0.01));
	if (E < 80) return;
	// none of these are working, and I don't know why
	double overall_n1 = math_chisquared_to_normal(overall_cs1, (cat - 1) << POSITIONS_L2);
	buf.str("");
	buf << get_name() << ":all1";
	results.push_back(TestResult(buf.str(), overall_n1, overall_n1, TestResult::TYPE_RAW_NORMAL, 0.01));
	buf.str("");
	buf << get_name() << ":all2";
	results.push_back(TestResult(buf.str(), overall_n2, overall_n2, TestResult::TYPE_RAW_NORMAL, 0.01));
	*/
}
void PractRand::Tests::TripleMirrorFreqN::test_blocks(TestBlock *data, int numblocks) {
	while (numblocks) {
		int skip_blocks = blocks_till_next_pass;
		if (skip_blocks > numblocks) skip_blocks = numblocks;
		if (skip_blocks) {
			numblocks -= skip_blocks;
			data += skip_blocks;
			blocks_tested += skip_blocks;
			blocks_till_next_pass -= skip_blocks;
		}
		if (!numblocks) return;

		for (int level = minimum_level; level < MAX_LEVELS; level++) {
			int ostate = level_state[level]++;
			if (ostate == 0) {
				saved_blocks[level] = data[0].as64[0];
			}
			else if (ostate == 1) {
				saved_blocks[level + MAX_LEVELS] = data[0].as64[0];
			}
			else {
				Uint64 old0, old1;
				if (ostate == 2) {
					old0 = saved_blocks[level];
					old1 = saved_blocks[level + MAX_LEVELS];
					saved_blocks[level] = data[0].as64[0];
				}
				else if (ostate == 3) {
					old0 = saved_blocks[level + MAX_LEVELS];
					old1 = saved_blocks[level];
					saved_blocks[level + MAX_LEVELS] = data[0].as64[0];
					level_state[level] = 2;
				}
				else {
					issue_error();
					old0 = 0;
					old1 = 0;
				}
				int base_index = level << TOTAL_INDEX_BITS;
				base_index |= old0 & ((1 << SIZE1) - 1);
				for (int position = 0; position < (1 << POSITIONS_L2); position++) {
					int index = base_index;
					index |= position << (SIZE1 + SIZE2 + SIZE3);
					int shift1 = position << ALIGN_L2;
					index |= ((old1 >> shift1) & ((1 << SIZE2) - 1)) << SIZE1;
					// this next part REQUIRES that SIZE3 be less than or equal to ALIGN, and that ALIGN be a power of 2
					int shift2 = shift1 << 1;
					index |= ((data[0].as64[shift2 >> 6] >> (shift2 & 63)) & ((1 << SIZE3) - 1)) << (SIZE1+SIZE2);
					counts.increment(index);
				}
			}
			level_polarity[level] ^= 1;
			if (!level_polarity[level]) break;
		}
		blocks_till_next_pass = blocks_per_pass;
	}
}






void PractRand::Tests::TripleMirrorCoup::init(PractRand::RNGs::vRNG *known_good) {
	counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next_pass = blocks_per_pass - 1;
	TestBaseclass::init(known_good);
	if (SIZE1 + SIZE2 + SIZE3 < 6) issue_error("TrippleMirrorCoup - we have a problem");
	for (int i = 0; i < 1 << (TOTAL_INDEX_BITS >> 6); i++) coup_masks[i] = 0;
	for (int i = 0; i < COUP_BUCKETS; i++) coup_counts[i] = 0;
	pass_number = 0;
	for (int i = 0; i < (1 << POSITIONS_L2); i++) coup_last[i] = 0;
	coup_collected = 0;
}
std::string PractRand::Tests::TripleMirrorCoup::get_name() const {
	std::ostringstream buf;
	buf << "TMC(" << passes_at_once;
	if (blocks_per_pass != 1) buf << "/" << blocks_per_pass;
	buf << ")";
	return buf.str();
}
void PractRand::Tests::TripleMirrorCoup::get_results(std::vector<TestResult> &results) {
	const Uint64 *counts_ = counts.get_array();
	int repeat_blocks = get_blocks_to_repeat();
	if (blocks_tested < repeat_blocks) return;
	Sint64 passes = ((blocks_tested - repeat_blocks) / blocks_per_pass) * passes_at_once;
	double E = passes * std::pow(0.5, SIZE1 + SIZE2 + SIZE3);
	if (E < 10) return;
	int worst_position = -1;
	double worst_n = 0, overall_cs1 = 0, overall_n2 = 0;
	int cat = 1 << (SIZE1 + SIZE2 + SIZE3);
	for (int pos = 0; pos < (1 << POSITIONS_L2); pos++) {
		double chisquared = g_test_flat(cat, &counts_[pos << (SIZE1 + SIZE2 + SIZE3)]);
		double n = math_chisquared_to_normal(chisquared, cat);
		if (std::fabs(n) >= std::fabs(worst_n)) {
			worst_n = n;
			worst_position = pos;
		}
		overall_cs1 += chisquared;
		overall_n2 += n * n;
	}
	overall_n2 = (std::sqrt(overall_n2 / (1 << POSITIONS_L2)) - 1.46) * std::pow(2.0, 0.5 * POSITIONS_L2);
	std::ostringstream buf;
	buf << get_name() << ":w(" << worst_position << ")";
	results.push_back(TestResult(buf.str(), worst_n, worst_n, TestResult::TYPE_RAW_NORMAL, 0.01));
	if (E < 80) return;
	// none of these are working, and I don't know why
	double overall_n1 = math_chisquared_to_normal(overall_cs1, (cat - 1) << POSITIONS_L2);
	buf.str("");
	buf << get_name() << ":all1";
	results.push_back(TestResult(buf.str(), overall_n1, overall_n1, TestResult::TYPE_RAW_NORMAL, 0.01));
	buf.str("");
	buf << get_name() << ":all2";
	results.push_back(TestResult(buf.str(), overall_n2, overall_n2, TestResult::TYPE_RAW_NORMAL, 0.01));
}
int PractRand::Tests::TripleMirrorCoup::get_blocks_to_repeat() const {
	enum { POSITIONS = 1 << POSITIONS_L2 };
	//int bytes_needed = ((1 << POSITIONS_L2) + passes_at_once) * SAMPLE_ALIGN + TestBlock::SIZE * 2 * BLOCK_STEP;
	int bytes_needed = TestBlock::SIZE * 2 * BLOCK_STEP + (POSITIONS << POSITION_ALIGN_L2) + (passes_at_once << BASE_ALIGN_L2);// if BLOCK_STEP and 0-based positions are used
	return (bytes_needed + TestBlock::SIZE - 1) >> TestBlock::SIZE_L2;
}
void PractRand::Tests::TripleMirrorCoup::test_blocks(TestBlock *data, int numblocks) {
	while (blocks_tested < get_blocks_to_repeat()) {
		if (!numblocks) return;
		data += 1;
		numblocks -= 1;
		blocks_tested += 1;
	}
	blocks_tested += numblocks;
	while (numblocks > blocks_till_next_pass) {
		data += blocks_till_next_pass;
		numblocks -= blocks_till_next_pass;
		pass_number++;
		enum { BASE_ALIGN = 1 << BASE_ALIGN_L2, POSITION_ALIGN = 1 << POSITION_ALIGN_L2, POSITIONS = 1 << POSITIONS_L2 };
		for (long pos_code = 0; pos_code < POSITIONS; pos_code++) {//0-based if BLOCK_STEP is used
			unsigned long base_index = pos_code << (SIZE1 + SIZE2 + SIZE3);//0-based if BLOCK_STEP is used
			long ofs = BLOCK_STEP * TestBlock::SIZE - pos_code * POSITION_ALIGN;
			for (long base_pos = -passes_at_once * BASE_ALIGN; base_pos < 0; base_pos += BASE_ALIGN) {
				unsigned long index = base_index;
				//Uint8 a = data[-BLOCK_STEP].as8[base_pos + 0], b = data[0].as8[base_pos + pos], c = data[-2 * BLOCK_STEP].as8[base_pos - pos];
				Uint8 a = data->as8[base_pos - 0], b = data->as8[base_pos - ofs], c = data->as8[base_pos - ofs - ofs];// if BLOCK_STEP is used
				//b -= a; a -= c; // nope.  these hurt more than they help, typically
				index |= ((unsigned long)(a & ((1 << SIZE1) - 1))) << (SIZE2 + SIZE3);
				index |= ((unsigned long)(b & ((1 << SIZE2) - 1))) << (SIZE3);
				index |= ((unsigned long)(c & ((1 << SIZE3) - 1))) << (0);
				counts.increment(index);
				if (0 == ~(coup_masks[index >> 6] |= (index & 63))) {
					//we *might* have completed a coupon set here
					bool completion = true;
					unsigned long region_end = (base_index + (1 << (SIZE1 + SIZE2 + SIZE3))) >> 6;
					for (unsigned long x = base_index >> 6; x < region_end; x++) if (~coup_masks[x]) completion = false;
					if (completion) {
						for (unsigned long x = base_index >> 6; x < region_end; x++) coup_masks[x] = 0;
					}
				}
			}
		}

		data += 1;
		numblocks -= 1;
		blocks_till_next_pass = blocks_per_pass - 1;
	}
	blocks_till_next_pass -= numblocks;
}


PractRand::Tests::SimpleClassifyWord::SimpleClassifyWord() {
	if (TOTAL_BITS > 30) issue_error("SimpleClassifyWord - TOTAL_BITS way too high");
}
void PractRand::Tests::SimpleClassifyWord::init(PractRand::RNGs::vRNG *known_good) {
	TestBaseclass::init(known_good);
	for (int i = 0; i < 8; i++) root_probs[i] = 0;
	double temp_probs[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < 256; i++) {
		int nib0 = count_ones8(i & 15);
		int nib1 = count_ones8(i >> 4);
		if (nib0 > 3) nib0 = 3;
		if (nib1 > 3) nib1 = 3;
		lookup_half[i] = nib0 + (nib1 << 2);
		temp_probs[nib0] += 1.0 / 256;
	}
	for (int i = 0; i < 256; i++) {
		int nib0 = (i >> 0) & 3;
		int nib1 = (i >> 2) & 3;
		int nib2 = (i >> 4) & 3;
		int nib3 = (i >> 6) & 3;
		if (!nib0 || !nib1 || !nib2 || !nib3) {//is there at least one nibble competing?
			if (!nib0 && (!nib1 || (nib0 + nib1) <= (nib2 + nib3))) second_lookup[i] = 0;
			else if (!nib1 && (nib0 + nib1) <= (nib2 + nib3)) second_lookup[i] = 1;
			else if (!nib2) second_lookup[i] = 2;
			else second_lookup[i] = 3;
		}
		else {//no nibbles competing
			if (nib0 + nib1 <= 2 && (nib0 + nib1 <= nib2 + nib3)) second_lookup[i] = 4;
			else if (nib2 + nib3 <= 2 && (nib0 + nib1 > nib2 + nib3)) second_lookup[i] = 5;
			else if (nib0 + nib1 + nib2 + nib3 < 8) second_lookup[i] = 6;
			else second_lookup[i] = 7;
		}
		root_probs[second_lookup[i]] += temp_probs[nib0] * temp_probs[nib1] * temp_probs[nib2] * temp_probs[nib3];
	}
	counts.reset_counts();
	current = 0;
	warmup = OVERLAP_WINDOW - 1;
	for (int i = 0; i < 1 << CODE_BITS; i++) {
		Uint32 r = i & 1;
		for (int x = 1; x < CODE_BITS; x++) r |= ((i >> x) & 1) << (x * OVERLAP_WINDOW);
		reorder_code[i] = r;
	}
	reorder_mask = (1ull << TOTAL_BITS) - 1 - reorder_code[(1 << CODE_BITS) - 1];
	for (int i = 0; i < 256; i++) {
		reordered_second_lookup[i] = reorder_code[second_lookup[i]];
	}
}
std::string PractRand::Tests::SimpleClassifyWord::get_name() const { return "SCW"; }
void PractRand::Tests::SimpleClassifyWord::get_results(std::vector<TestResult> &results) {
	double expected_total = blocks_tested * (TestBlock::SIZE / 2) - (OVERLAP_WINDOW - 1);
	std::vector<double> probs;
	probs.resize(1 << TOTAL_BITS);
	std::vector<Uint64> tmp_counts;
	tmp_counts.resize(1 << TOTAL_BITS);
	counts.flush();
	Uint64 tsize = 1ull << TOTAL_BITS;
	Uint64 observed_total = 0;
	for (int i = 0; i < tsize; i++) {
		double p = 1.0;
		enum { SINGLE_CODE_MASK = (1 << CODE_BITS) - 1 };
		int reordered_index = 0;
		for (int r = 0; r < OVERLAP_WINDOW; r++) {
			p *= root_probs[(i >> (r * CODE_BITS)) & SINGLE_CODE_MASK];
			reordered_index |= reorder_code[(i >> (r * CODE_BITS)) & SINGLE_CODE_MASK] << r;
		}
		probs[i] = p;
		if (reordered_index >= tsize) issue_error("SimpleClassifyWord::get_results - reordered index too large");
		tmp_counts[i] = counts[reordered_index];
		observed_total += tmp_counts[i];
	}
	Uint64 reduced_size1 = tsize;
	//to do: shrink the overlap window here if the test was too short

	Uint64 reduced_size2 = simplify_prob_table(reduced_size1, blocks_tested * (TestBlock::SIZE / 2) / 25000.0,
		&probs[0], &tmp_counts[0], true, true);
	double r = g_test(reduced_size2, &probs[0], &tmp_counts[0]);
	double n = math_chisquared_to_normal(r, reduced_size2 - 1);
	results.push_back(TestResult(get_name(), n, 0, TestResult::TYPE_RAW_NORMAL, 0.1));
}
Uint8 PractRand::Tests::SimpleClassifyWord::word_to_code(Uint16 word) {
	Uint8 byte0 = word & 255;
	Uint8 byte1 = word >> 8;
	Uint8 intermediate0 = lookup_half[byte0];
	Uint8 intermediate1 = lookup_half[byte1];
	Uint8 merged = intermediate0 + (intermediate1 << 4);
	Uint8 code = second_lookup[merged];
	return code;
	//return second_lookup[lookup_half[word & 255] | (lookup_half[word >> 8] << 4)];
}
Uint32 PractRand::Tests::SimpleClassifyWord::word_to_reordered_code(Uint16 word) {
	Uint8 byte0 = word & 255;
	Uint8 byte1 = word >> 8;
	Uint8 intermediate0 = lookup_half[byte0];
	Uint8 intermediate1 = lookup_half[byte1];
	Uint8 merged = intermediate0 + (intermediate1 << 4);
	if (true) {
		Uint8 code = second_lookup[merged];
		Uint32 reordered = reorder_code[code];
		return reordered;
	}
	else {
		Uint32 reordered = reordered_second_lookup[merged];
		return reordered;
	}
}
void PractRand::Tests::SimpleClassifyWord::test_blocks(TestBlock *data, int numblocks) {
	Uint16 *pos = data->as16;
	Uint16 *end = pos + numblocks * (TestBlock::SIZE / 2);
	if (pos == end) return;

	blocks_tested += numblocks;

	while (warmup) {
		Uint16 word = *(pos++);
		Uint32 reordered_code = word_to_reordered_code(word);
		current <<= 1;
		current &= reorder_mask;
		current |= reordered_code;
		warmup--;
	}
	while (pos < end) {
		Uint16 word = *(pos++);
		Uint32 reordered_code = word_to_reordered_code(word);
		current <<= 1;
		current &= reorder_mask;
		current |= reordered_code;
		counts.increment(current);
	}
}

static const double _lperm16_ratio[] = {// 32768 / 2 array elements ; from 1 petabyte of hc256 (seed=0x5de6fa19): 
	0.50003673, 0.54964384, 0.59699571, 0.62975592, 0.61430835, 0.64362826, 0.65492433, 0.66203379, 0.55929909, 0.61887826, 0.67023886, 0.71175169, 0.69277723, 0.72800420, 0.74475900, 0.75474994, 0.76672578, 0.81303166, 0.84287851, 0.87317403, 0.86025155, 0.88318222, 0.89561472, 0.90238026, 0.82005262, 0.86693340, 0.89244648, 0.92198983, 0.91027262, 0.93026958, 0.94323935, 0.94931757, 0.78575513, 0.83364011, 0.86254753, 0.89400431, 0.88097009, 0.90371114, 0.91755396, 0.92457768, 0.89999686, 0.93014751, 0.94344012, 0.96154733, 0.95486285, 0.96594631, 0.97374015, 0.97706186, 0.91420083, 0.94357850, 0.95526494, 0.97237443, 0.96647448, 0.97602075, 0.98361379, 0.98632569, 0.96150139, 0.97731262, 0.98208122, 0.99031430, 0.98783217, 0.99168431, 0.99507128, 0.99603861, //    0-   63
	0.83334451, 0.85309763, 0.92223448, 0.94000808, 0.92857594, 0.96666018, 0.97140476, 0.98718621, 0.89198597, 0.90842706, 0.95456703, 0.96755369, 0.95951481, 0.98281718, 0.98609171, 0.99434811, 0.94793381, 0.95611191, 0.97533466, 0.98252472, 0.97818226, 0.98991274, 0.99193717, 0.99650169, 0.97617128, 0.98064967, 0.98946558, 0.99317484, 0.99104753, 0.99618083, 0.99720771, 0.99890210, 0.96286684, 0.96942135, 0.98343535, 0.98894767, 0.98571369, 0.99386986, 0.99538793, 0.99818112, 0.98811929, 0.99023229, 0.99412764, 0.99616173, 0.99501388, 0.99770212, 0.99832969, 0.99929862, 0.99449678, 0.99558824, 0.99741608, 0.99844609, 0.99788984, 0.99910269, 0.99941783, 0.99978430, 0.99850946, 0.99877696, 0.99920395, 0.99951407, 0.99935050, 0.99969648, 0.99980677, 0.99992227, //   64-  127
	0.45042098, 0.50003672, 0.54737718, 0.58179149, 0.56556635, 0.59635681, 0.60867794, 0.61641398, 0.50966471, 0.57066189, 0.62322925, 0.66805727, 0.64755445, 0.68557419, 0.70438123, 0.71557236, 0.71710123, 0.76876770, 0.80197115, 0.83785319, 0.82255770, 0.84969410, 0.86509920, 0.87347104, 0.77651517, 0.83066211, 0.86009198, 0.89656594, 0.88210330, 0.90679705, 0.92362408, 0.93150642, 0.73832448, 0.79238375, 0.82502396, 0.86291681, 0.84722780, 0.87459503, 0.89205645, 0.90092582, 0.86570934, 0.90329440, 0.91989018, 0.94413015, 0.93518031, 0.95003958, 0.96109857, 0.96578759, 0.88339076, 0.92075231, 0.93561855, 0.95910640, 0.95100429, 0.96412217, 0.97519880, 0.97915283, 0.94249492, 0.96464491, 0.97129191, 0.98391897, 0.98013678, 0.98601989, 0.99156399, 0.99314594, //  128-  191
	0.79276320, 0.81518659, 0.89367194, 0.91608406, 0.90167802, 0.94967031, 0.95631669, 0.97852577, 0.86146024, 0.88087069, 0.93532314, 0.95253078, 0.94189031, 0.97272656, 0.97762118, 0.98995057, 0.92702145, 0.93748321, 0.96209300, 0.97233166, 0.96615185, 0.98288991, 0.98609592, 0.99338188, 0.96498593, 0.97102239, 0.98290345, 0.98853307, 0.98530636, 0.99309076, 0.99484817, 0.99776037, 0.94709375, 0.95568089, 0.97402457, 0.98211402, 0.97737347, 0.98931702, 0.99182946, 0.99644976, 0.98104113, 0.98411856, 0.98977308, 0.99307369, 0.99121757, 0.99557046, 0.99671460, 0.99849351, 0.99100250, 0.99263455, 0.99537231, 0.99710165, 0.99616817, 0.99820900, 0.99881379, 0.99951429, 0.99725976, 0.99770484, 0.99841379, 0.99898948, 0.99868626, 0.99932811, 0.99956211, 0.99980616, //  192-  255
	0.40304063, 0.45267368, 0.50003576, 0.53599494, 0.51904020, 0.55123842, 0.56452091, 0.57286482, 0.46230575, 0.52463506, 0.57836422, 0.62634158, 0.60440312, 0.64510124, 0.66582849, 0.67817056, 0.66976199, 0.72648651, 0.76295480, 0.80414311, 0.78655104, 0.81771384, 0.83598522, 0.84588280, 0.73495657, 0.79603408, 0.82921688, 0.87229334, 0.85519382, 0.88438483, 0.90488819, 0.91449282, 0.69303328, 0.75301886, 0.78921061, 0.83323327, 0.81499712, 0.84679596, 0.86771576, 0.87833637, 0.83292221, 0.87765661, 0.89734543, 0.92750557, 0.91640460, 0.93485273, 0.94901762, 0.95505593, 0.85395152, 0.89892883, 0.91688554, 0.94646055, 0.93624098, 0.95275747, 0.96715966, 0.97230962, 0.92422758, 0.95256590, 0.96103708, 0.97781066, 0.97278392, 0.98060568, 0.98821666, 0.99039291, //  256-  319
	0.75403442, 0.77900413, 0.86641155, 0.89325209, 0.87599554, 0.93345641, 0.94192817, 0.97025088, 0.83231783, 0.85456515, 0.91695844, 0.93819089, 0.92506786, 0.96308781, 0.96954191, 0.98575743, 0.90706720, 0.91970264, 0.94952202, 0.96262373, 0.95467999, 0.97617119, 0.98053072, 0.99039294, 0.95430821, 0.96183089, 0.97664165, 0.98410139, 0.97983014, 0.99013953, 0.99260049, 0.99667512, 0.93204526, 0.94256767, 0.96503594, 0.97558024, 0.96941486, 0.98497577, 0.98843809, 0.99478830, 0.97430570, 0.97829031, 0.98561871, 0.99011938, 0.98758983, 0.99355520, 0.99516865, 0.99772259, 0.98766700, 0.98981346, 0.99341868, 0.99582068, 0.99452555, 0.99735740, 0.99823653, 0.99925596, 0.99606766, 0.99668136, 0.99765865, 0.99848625, 0.99805100, 0.99897592, 0.99932841, 0.99969590, //  320-  383
	0.37030652, 0.41825197, 0.46403611, 0.50002646, 0.48305459, 0.51524535, 0.52892453, 0.53752037, 0.42758381, 0.48883944, 0.54163624, 0.59072484, 0.56827864, 0.60987759, 0.63181851, 0.64488201, 0.62803966, 0.68681098, 0.72457648, 0.76921998, 0.75016376, 0.78399892, 0.80461558, 0.81583128, 0.69562838, 0.76059481, 0.79592742, 0.84443978, 0.82519889, 0.85806393, 0.88229280, 0.89366389, 0.65220955, 0.71497830, 0.75283372, 0.80133635, 0.78125027, 0.81630419, 0.84039802, 0.85263787, 0.79691747, 0.84719461, 0.86936161, 0.90538984, 0.89205445, 0.91414862, 0.93206114, 0.93965673, 0.82060421, 0.87204773, 0.89252256, 0.92875970, 0.91626046, 0.93647780, 0.95523589, 0.96192658, 0.89955844, 0.93455980, 0.94505560, 0.96748951, 0.96076267, 0.97121987, 0.98213693, 0.98525455, //  384-  447
	0.71673803, 0.74300232, 0.83493262, 0.86532502, 0.84576999, 0.91088595, 0.92142264, 0.95642755, 0.80145542, 0.82565652, 0.89350756, 0.91866560, 0.90310858, 0.94819350, 0.95663592, 0.97786777, 0.88234891, 0.89689069, 0.93116266, 0.94766061, 0.93770668, 0.96462977, 0.97065637, 0.98431142, 0.93943360, 0.94853341, 0.96644634, 0.97636675, 0.97067878, 0.98442113, 0.98807317, 0.99413458, 0.91250998, 0.92491900, 0.95138599, 0.96498377, 0.95701717, 0.97712112, 0.98206576, 0.99118605, 0.96359093, 0.96869494, 0.97808593, 0.98441546, 0.98085451, 0.98923070, 0.99178339, 0.99577354, 0.98205190, 0.98488448, 0.98965028, 0.99316800, 0.99127068, 0.99542070, 0.99687162, 0.99855539, 0.99364678, 0.99452850, 0.99593502, 0.99726362, 0.99656540, 0.99804813, 0.99868558, 0.99935095, //  448-  511
	0.38572943, 0.43446728, 0.48099915, 0.51699067, 0.50002299, 0.53220868, 0.54570505, 0.55417261, 0.44394173, 0.50572131, 0.55894460, 0.60751323, 0.58530234, 0.62647175, 0.64785069, 0.66058114, 0.64771086, 0.70553310, 0.74269033, 0.78569147, 0.76732509, 0.79992129, 0.81938967, 0.82998123, 0.71419868, 0.77730110, 0.81160729, 0.85758145, 0.83934198, 0.87046692, 0.89295100, 0.90348852, 0.67143595, 0.73290501, 0.77001292, 0.81637572, 0.79716470, 0.83070359, 0.85327705, 0.86475327, 0.81388365, 0.86155256, 0.88257831, 0.91581639, 0.90354614, 0.92389683, 0.94005501, 0.94691838, 0.83633088, 0.88472854, 0.90401520, 0.93710413, 0.92568098, 0.94416013, 0.96086255, 0.96681913, 0.91122531, 0.94302207, 0.95256055, 0.97236156, 0.96642767, 0.97565225, 0.98500258, 0.98767158, //  512-  575
	0.73431485, 0.75996453, 0.84975234, 0.87848649, 0.86003261, 0.92156653, 0.93110184, 0.96296032, 0.81601393, 0.83929301, 0.90456654, 0.92787267, 0.91346980, 0.95521036, 0.96272500, 0.98157524, 0.89400315, 0.90764235, 0.93979446, 0.95471946, 0.94570975, 0.97008600, 0.97531527, 0.98719092, 0.94645275, 0.95480222, 0.97125218, 0.98001642, 0.97498904, 0.98711457, 0.99021176, 0.99533535, 0.92172023, 0.93323906, 0.95781839, 0.96997680, 0.96285722, 0.98082388, 0.98507027, 0.99288693, 0.96863727, 0.97321258, 0.98163209, 0.98710684, 0.98402676, 0.99126110, 0.99338511, 0.99668809, 0.98469956, 0.98721088, 0.99142953, 0.99441652, 0.99280545, 0.99633442, 0.99751486, 0.99888532, 0.99478732, 0.99554366, 0.99674928, 0.99784104, 0.99726496, 0.99848725, 0.99898855, 0.99951322, //  576-  639
	0.35643015, 0.40370039, 0.44882427, 0.48478228, 0.46782277, 0.50004642, 0.51387695, 0.52257776, 0.41288568, 0.47368356, 0.52609645, 0.57566330, 0.55300492, 0.59499453, 0.61742522, 0.63078668, 0.61044582, 0.67003571, 0.70833263, 0.75446768, 0.73475947, 0.76971366, 0.79134267, 0.80308654, 0.67896778, 0.74562785, 0.78183123, 0.83266424, 0.81250571, 0.84688711, 0.87274171, 0.88485780, 0.63493489, 0.69888204, 0.73750656, 0.78783682, 0.76696659, 0.80338680, 0.82883513, 0.84175968, 0.78154611, 0.83428340, 0.85744520, 0.89603023, 0.88182726, 0.90540180, 0.92488069, 0.93314936, 0.80647655, 0.86065314, 0.88222982, 0.92128882, 0.90782194, 0.92960144, 0.95019242, 0.95753348, 0.88907485, 0.92693947, 0.93832048, 0.96312339, 0.95567129, 0.96727025, 0.97956538, 0.98308752, //  640-  703
	0.70096751, 0.72776066, 0.82159780, 0.85350886, 0.83299145, 0.90137804, 0.91276674, 0.95061105, 0.78841134, 0.81341522, 0.88355642, 0.91041637, 0.89380211, 0.94191498, 0.95118874, 0.97453567, 0.87187764, 0.88723117, 0.92344609, 0.94132874, 0.93052579, 0.95979608, 0.96650018, 0.98171785, 0.93314195, 0.94290749, 0.96210696, 0.97309345, 0.96680863, 0.98198993, 0.98615909, 0.99306556, 0.90426527, 0.91745089, 0.94561218, 0.96050689, 0.95178281, 0.97376994, 0.97936828, 0.98966268, 0.95904828, 0.96464602, 0.97488517, 0.98200737, 0.97800442, 0.98739338, 0.99034248, 0.99494282, 0.97967405, 0.98280115, 0.98805373, 0.99204694, 0.98989506, 0.99460295, 0.99629385, 0.99825949, 0.99262287, 0.99361856, 0.99520357, 0.99674694, 0.99593579, 0.99765887, 0.99841413, 0.99920516, //  704-  767
	0.34511680, 0.39136780, 0.43552998, 0.47111678, 0.45433378, 0.48617329, 0.50001457, 0.50871540, 0.40036635, 0.46019997, 0.51175889, 0.56115297, 0.53856060, 0.58043377, 0.60312116, 0.61663974, 0.59375246, 0.65322998, 0.69140574, 0.73836493, 0.71833966, 0.75386896, 0.77628490, 0.78847225, 0.66212797, 0.72940973, 0.76598511, 0.81855394, 0.79770326, 0.83330300, 0.86067318, 0.87350164, 0.61818456, 0.68229542, 0.72099221, 0.77258866, 0.75121476, 0.78851345, 0.81514950, 0.82868980, 0.76466527, 0.81884976, 0.84270073, 0.88353645, 0.86845662, 0.89347136, 0.91472544, 0.92374668, 0.79019450, 0.84647531, 0.86886494, 0.91078346, 0.89631649, 0.91969947, 0.94261596, 0.95078177, 0.87523651, 0.91592507, 0.92813949, 0.95592287, 0.94757666, 0.96056464, 0.97498450, 0.97910689, //  768-  831
	0.68470366, 0.71164558, 0.80593306, 0.83897285, 0.81773938, 0.88849425, 0.90067493, 0.94135978, 0.77371351, 0.79922095, 0.87074936, 0.89911242, 0.88157304, 0.93238325, 0.94269887, 0.96859605, 0.85868368, 0.87469716, 0.91242256, 0.93180418, 0.92010955, 0.95177034, 0.95942822, 0.97677777, 0.92423560, 0.93467140, 0.95521458, 0.96751697, 0.96046355, 0.97749207, 0.98245973, 0.99069765, 0.89333792, 0.90724701, 0.93694821, 0.95335260, 0.94374360, 0.96797783, 0.97450732, 0.98652262, 0.95196791, 0.95808307, 0.96934524, 0.97753288, 0.97290662, 0.98373716, 0.98736852, 0.99301464, 0.97569039, 0.97919079, 0.98506847, 0.98979788, 0.98724609, 0.99282558, 0.99498739, 0.99749582, 0.99058874, 0.99175729, 0.99361962, 0.99554445, 0.99452883, 0.99668046, 0.99770429, 0.99877681, //  832-  895
	0.33801418, 0.38363673, 0.42717028, 0.46251909, 0.44585950, 0.47746331, 0.49131619, 0.50002068, 0.39250926, 0.45171358, 0.50274864, 0.55204526, 0.52949316, 0.57129158, 0.59414244, 0.60774624, 0.58325938, 0.64266121, 0.68082954, 0.72824884, 0.70800784, 0.74390603, 0.76682237, 0.77927636, 0.65157473, 0.71925437, 0.75605495, 0.80968397, 0.78840978, 0.82473189, 0.85308658, 0.86637527, 0.60766792, 0.67189292, 0.71066419, 0.76300993, 0.74132756, 0.77916478, 0.80655605, 0.82046855, 0.75402900, 0.80913546, 0.83344550, 0.87569528, 0.86006720, 0.88594590, 0.90833726, 0.91784917, 0.77996867, 0.83757840, 0.86047536, 0.90419290, 0.88910715, 0.91349545, 0.93785363, 0.94654038, 0.86655853, 0.90899328, 0.92178188, 0.95141003, 0.94249577, 0.95635638, 0.97210827, 0.97660280, //  896-  959
	0.67448263, 0.70151768, 0.79611388, 0.82981498, 0.80814373, 0.88032075, 0.89307768, 0.93552453, 0.76450006, 0.79030501, 0.86269209, 0.89202256, 0.87388579, 0.92641431, 0.93735449, 0.96486974, 0.85042071, 0.86681445, 0.90552481, 0.92583262, 0.91356801, 0.94674725, 0.95499221, 0.97366936, 0.91864293, 0.92949890, 0.95086056, 0.96400824, 0.95647422, 0.97465852, 0.98013296, 0.98920597, 0.88647588, 0.90084004, 0.93150136, 0.94886616, 0.93870208, 0.96434007, 0.97145362, 0.98454842, 0.94751451, 0.95396881, 0.96584185, 0.97472455, 0.96970848, 0.98144511, 0.98549524, 0.99179668, 0.97318491, 0.97692192, 0.98319196, 0.98838619, 0.98558462, 0.99171006, 0.99416642, 0.99701509, 0.98931312, 0.99058744, 0.99262099, 0.99478748, 0.99364590, 0.99606709, 0.99725995, 0.99850780, //  960- 1023
	0.44077059, 0.49036840, 0.53771801, 0.57245771, 0.55609709, 0.58716295, 0.59967818, 0.60754762, 0.50003822, 0.56128380, 0.61408937, 0.65956261, 0.63876344, 0.67730835, 0.69652069, 0.70795922, 0.70743273, 0.76015873, 0.79405686, 0.83099384, 0.81520966, 0.84318351, 0.85917358, 0.86785772, 0.76803809, 0.82361338, 0.85382852, 0.89161982, 0.87661370, 0.90223304, 0.91980681, 0.92803059, 0.72909890, 0.78436845, 0.81772309, 0.85687165, 0.84064015, 0.86893398, 0.88709492, 0.89632232, 0.85900042, 0.89807931, 0.91531309, 0.94075794, 0.93133896, 0.94692329, 0.95863826, 0.96360101, 0.87738866, 0.91631329, 0.93180596, 0.95653023, 0.94799063, 0.96181236, 0.97356190, 0.97775351, 0.93873257, 0.96216816, 0.96921498, 0.98267845, 0.97862876, 0.98492836, 0.99088007, 0.99258496, // 1024- 1087
	0.78487835, 0.80782233, 0.88814255, 0.91143783, 0.89644349, 0.94636603, 0.95338292, 0.97683125, 0.85552382, 0.87551218, 0.93157831, 0.94961081, 0.93846730, 0.97076123, 0.97597706, 0.98909046, 0.92297873, 0.93385928, 0.95953165, 0.97036134, 0.96382406, 0.98150546, 0.98496710, 0.99277218, 0.96281270, 0.96914897, 0.98162882, 0.98762822, 0.98419190, 0.99248804, 0.99439241, 0.99753994, 0.94402540, 0.95301235, 0.97219120, 0.98077895, 0.97575338, 0.98843569, 0.99113590, 0.99610945, 0.97967065, 0.98292977, 0.98893019, 0.99247376, 0.99047870, 0.99515786, 0.99639858, 0.99833583, 0.99032365, 0.99205916, 0.99497504, 0.99684065, 0.99583482, 0.99803670, 0.99869534, 0.99946204, 0.99701672, 0.99749567, 0.99825876, 0.99888585, 0.99855644, 0.99925590, 0.99951390, 0.99978412, // 1088- 1151
	0.38116189, 0.42938627, 0.47541008, 0.51119367, 0.49432563, 0.52633571, 0.53984462, 0.54832904, 0.43875030, 0.50002093, 0.55281647, 0.60134463, 0.57915353, 0.62029613, 0.64179875, 0.65460919, 0.64033817, 0.69812468, 0.73526821, 0.77876079, 0.76021471, 0.79312270, 0.81303089, 0.82385145, 0.70680221, 0.77031563, 0.80483339, 0.85172160, 0.83313867, 0.86488794, 0.88809084, 0.89896668, 0.66405621, 0.72568682, 0.76288202, 0.80991543, 0.79041356, 0.82441353, 0.84759908, 0.85937812, 0.80650421, 0.85504584, 0.87643851, 0.91085567, 0.89813868, 0.91921714, 0.93616819, 0.94335523, 0.82933948, 0.87887899, 0.89858428, 0.93304897, 0.92115434, 0.94039003, 0.95806563, 0.96437039, 0.90561286, 0.93883174, 0.94878065, 0.96983271, 0.96352126, 0.97334167, 0.98348643, 0.98638480, // 1152- 1215
	0.72724089, 0.75299101, 0.84311968, 0.87244496, 0.85357657, 0.91641485, 0.92636155, 0.95961186, 0.80980649, 0.83334734, 0.89938718, 0.92343988, 0.90856320, 0.95168281, 0.95963118, 0.97960696, 0.88861673, 0.90259310, 0.93553015, 0.95116181, 0.94174077, 0.96725874, 0.97286697, 0.98560606, 0.94304014, 0.95171346, 0.96878185, 0.97811056, 0.97276434, 0.98566332, 0.98905070, 0.99465703, 0.91738153, 0.92925463, 0.95460365, 0.96742548, 0.95990668, 0.97885332, 0.98346595, 0.99194666, 0.96606789, 0.97089160, 0.97976410, 0.98566654, 0.98233567, 0.99014046, 0.99249916, 0.99616961, 0.98332673, 0.98599518, 0.99047093, 0.99373455, 0.99197562, 0.99582326, 0.99715340, 0.99869566, 0.99416563, 0.99498730, 0.99629425, 0.99751450, 0.99687159, 0.99823553, 0.99881353, 0.99941805, // 1216- 1279
	0.32980702, 0.37681417, 0.42169809, 0.45841589, 0.44109171, 0.47391334, 0.48828566, 0.49730214, 0.38595759, 0.44720922, 0.50003281, 0.55116870, 0.52778123, 0.57113127, 0.59463816, 0.60865489, 0.58247061, 0.64465912, 0.68461079, 0.73377714, 0.71278197, 0.74994952, 0.77326504, 0.78594582, 0.65401141, 0.72441558, 0.76268247, 0.81733980, 0.79563963, 0.83270824, 0.86075415, 0.87391377, 0.60805563, 0.67511301, 0.71558177, 0.76945144, 0.74713691, 0.78605550, 0.81356202, 0.82754093, 0.76116489, 0.81789178, 0.84296698, 0.88509856, 0.86952695, 0.89529350, 0.91677836, 0.92589970, 0.78799051, 0.84663781, 0.86998134, 0.91280322, 0.89802685, 0.92192847, 0.94470239, 0.95282587, 0.87698152, 0.91871807, 0.93121381, 0.95876885, 0.95049543, 0.96338626, 0.97710582, 0.98104132, // 1280- 1343
	0.67757221, 0.70576761, 0.80440917, 0.83884241, 0.81667588, 0.89055555, 0.90307699, 0.94471501, 0.77039676, 0.79700357, 0.87163505, 0.90088537, 0.88280139, 0.93522937, 0.94553598, 0.97143765, 0.85897834, 0.87563442, 0.91487407, 0.93461325, 0.92269505, 0.95499665, 0.96244552, 0.97943195, 0.92600979, 0.93668114, 0.95770585, 0.96990007, 0.96291495, 0.97976825, 0.98444565, 0.99218466, 0.89441181, 0.90878472, 0.93943441, 0.95591403, 0.94624028, 0.97061043, 0.97684905, 0.98835877, 0.95435173, 0.96051930, 0.97186124, 0.97979521, 0.97531683, 0.98581366, 0.98913346, 0.99429460, 0.97729862, 0.98076686, 0.98659142, 0.99105625, 0.98864796, 0.99391695, 0.99582411, 0.99803340, 0.99170764, 0.99282636, 0.99459807, 0.99633342, 0.99541983, 0.99735419, 0.99820870, 0.99910314, // 1344- 1407
	0.28827429, 0.33196811, 0.37370413, 0.40931273, 0.39252379, 0.42437393, 0.43888500, 0.44798115, 0.34047276, 0.39869495, 0.44886116, 0.50002210, 0.47662467, 0.51997989, 0.54459666, 0.55925343, 0.52321634, 0.58579850, 0.62602021, 0.67854500, 0.65613193, 0.69588928, 0.72213002, 0.73641709, 0.59519303, 0.66875217, 0.70871710, 0.77003426, 0.74572273, 0.78725007, 0.82081168, 0.83653520, 0.54889886, 0.61747206, 0.65884310, 0.71760853, 0.69325065, 0.73571968, 0.76758037, 0.78375946, 0.70315838, 0.76607610, 0.79386615, 0.84420957, 0.82559850, 0.85641475, 0.88400214, 0.89569777, 0.73281708, 0.79945207, 0.82595723, 0.87884829, 0.86059155, 0.89013303, 0.92057682, 0.93143533, 0.83165290, 0.88325822, 0.89877566, 0.93627142, 0.92502421, 0.94253647, 0.96301944, 0.96887709, // 1408- 1471
	0.62107508, 0.65007683, 0.75157935, 0.79030616, 0.76542814, 0.84840651, 0.86390274, 0.91552178, 0.72044428, 0.74907176, 0.82936073, 0.86415225, 0.84264212, 0.90496358, 0.91866911, 0.95312886, 0.81528645, 0.83438852, 0.87936814, 0.90426485, 0.88925129, 0.92993250, 0.94047124, 0.96442952, 0.89727164, 0.91030903, 0.93597779, 0.95258012, 0.94306052, 0.96601792, 0.97323600, 0.98517370, 0.85863094, 0.87559735, 0.91182525, 0.93342756, 0.92076431, 0.95267936, 0.96193332, 0.97894781, 0.93195343, 0.93996052, 0.95467526, 0.96612979, 0.95967315, 0.97478947, 0.98021139, 0.98863274, 0.96489360, 0.96959511, 0.97748992, 0.98429009, 0.98062139, 0.98864958, 0.99197341, 0.99583391, 0.98558171, 0.98724676, 0.98989167, 0.99280397, 0.99126979, 0.99452501, 0.99616813, 0.99789040, // 1472- 1535
	0.30726858, 0.35248193, 0.39563183, 0.43175036, 0.41472989, 0.44703971, 0.46146202, 0.47053211, 0.36126077, 0.42088122, 0.47224970, 0.52340656, 0.50001270, 0.54336177, 0.56748054, 0.58183478, 0.55030981, 0.61271890, 0.65284937, 0.70378991, 0.68203267, 0.72057114, 0.74551683, 0.75906140, 0.62208378, 0.69420302, 0.73339567, 0.79166528, 0.76855422, 0.80801473, 0.83907393, 0.85362816, 0.57594600, 0.64381408, 0.68480357, 0.74130118, 0.71788099, 0.75873783, 0.78860367, 0.80377969, 0.72969995, 0.78981666, 0.81631600, 0.86290178, 0.84566882, 0.87422462, 0.89899381, 0.90950729, 0.75803967, 0.82101292, 0.84608488, 0.89438251, 0.87771180, 0.90467789, 0.93160914, 0.94121656, 0.85240600, 0.89944837, 0.91357688, 0.94654739, 0.93667273, 0.95205156, 0.96946351, 0.97443872, // 1536- 1599
	0.64690781, 0.67553107, 0.77572317, 0.81249920, 0.78885035, 0.86768156, 0.88180348, 0.92886988, 0.74328060, 0.77099326, 0.84869164, 0.88094621, 0.86099538, 0.91880345, 0.93095785, 0.96151265, 0.83528906, 0.85324334, 0.89560449, 0.91814166, 0.90454096, 0.94137286, 0.95051504, 0.97129459, 0.91040989, 0.92237210, 0.94591997, 0.96049661, 0.95214102, 0.97232436, 0.97836086, 0.98837869, 0.87499021, 0.89077124, 0.92444508, 0.94371047, 0.93241978, 0.96087311, 0.96875716, 0.98325238, 0.94219610, 0.94935516, 0.96254620, 0.97237065, 0.96682121, 0.97983249, 0.98429003, 0.99122591, 0.97056605, 0.97470242, 0.98165189, 0.98738520, 0.98429318, 0.99105640, 0.99373401, 0.99683978, 0.98838187, 0.98979625, 0.99204646, 0.99441798, 0.99316718, 0.99581962, 0.99710144, 0.99844425, // 1600- 1663
	0.27206214, 0.31446954, 0.35497880, 0.39014040, 0.37357639, 0.40504177, 0.41959903, 0.42875137, 0.32270506, 0.37974206, 0.42890868, 0.48003658, 0.45665549, 0.50002433, 0.52506120, 0.53996533, 0.50001276, 0.56279730, 0.60319991, 0.65697663, 0.63402844, 0.67476260, 0.70217012, 0.71708705, 0.57221330, 0.64701289, 0.68769135, 0.75156769, 0.72624084, 0.76950901, 0.80521758, 0.82194327, 0.52580623, 0.59495532, 0.63671689, 0.69734500, 0.67221781, 0.71605006, 0.74962707, 0.76668157, 0.68049734, 0.74583614, 0.77470343, 0.82825555, 0.80841207, 0.84125605, 0.87120818, 0.88388765, 0.71131963, 0.78103022, 0.80879984, 0.86560455, 0.84598289, 0.87768598, 0.91115491, 0.92308441, 0.81401201, 0.86941071, 0.88603257, 0.92749070, 0.91506661, 0.93441214, 0.95751822, 0.96412544, // 1664- 1727
	0.59901736, 0.62833570, 0.73096876, 0.77136020, 0.74543245, 0.83189766, 0.84856002, 0.90416769, 0.70093510, 0.73036491, 0.81288824, 0.84981020, 0.82695496, 0.89316734, 0.90819935, 0.94601876, 0.79822160, 0.81828423, 0.86548142, 0.89242202, 0.87616503, 0.92017316, 0.93190318, 0.95854105, 0.88605599, 0.90001741, 0.92748711, 0.94581833, 0.93531159, 0.96066746, 0.96885971, 0.98244124, 0.84466364, 0.86264214, 0.90104091, 0.92464670, 0.91080671, 0.94568992, 0.95610369, 0.97528346, 0.92320649, 0.93192893, 0.94795127, 0.96078865, 0.95355934, 0.97052327, 0.97673594, 0.98642534, 0.96004985, 0.96523371, 0.97394249, 0.98165296, 0.97748876, 0.98659481, 0.99047337, 0.99497185, 0.98319194, 0.98506784, 0.98805610, 0.99142677, 0.98964861, 0.99341850, 0.99537181, 0.99741819, // 1728- 1791
	0.25527671, 0.29565838, 0.33420737, 0.36821231, 0.35218500, 0.38261181, 0.39690807, 0.40589632, 0.30351323, 0.35822824, 0.40538809, 0.45543050, 0.43255013, 0.47497241, 0.50001419, 0.51492747, 0.47236557, 0.53371259, 0.57316546, 0.62708419, 0.60407839, 0.64489012, 0.67317165, 0.68855155, 0.54293450, 0.61720911, 0.65758972, 0.72322044, 0.69719711, 0.74162351, 0.77968357, 0.79752081, 0.49755208, 0.56559534, 0.60665811, 0.66812045, 0.64265167, 0.68710069, 0.72223051, 0.74008522, 0.64879530, 0.71520719, 0.74449485, 0.80111746, 0.78019000, 0.81488582, 0.84790043, 0.86191772, 0.68008580, 0.75187545, 0.78045262, 0.84171259, 0.82057324, 0.85475105, 0.89273407, 0.90627511, 0.78439616, 0.84402967, 0.86192557, 0.90897057, 0.89486486, 0.91681495, 0.94476261, 0.95275533, // 1792- 1855
	0.57009657, 0.59906340, 0.70044947, 0.74176542, 0.71521185, 0.80377447, 0.82165191, 0.88133113, 0.67276482, 0.70244395, 0.78567283, 0.82467538, 0.80055558, 0.87042651, 0.88731251, 0.92973035, 0.77075043, 0.79151676, 0.84045241, 0.86972588, 0.85206014, 0.89984323, 0.91343412, 0.94424658, 0.86543681, 0.88044566, 0.90995687, 0.93086237, 0.91888521, 0.94780212, 0.95790028, 0.97464557, 0.82079621, 0.83978351, 0.88030062, 0.90661992, 0.89120046, 0.93007217, 0.94255807, 0.96551926, 0.90549560, 0.91514821, 0.93291677, 0.94801107, 0.93950271, 0.95944540, 0.96737732, 0.97971157, 0.94936497, 0.95528738, 0.96523417, 0.97470418, 0.96959553, 0.98076888, 0.98599463, 0.99205958, 0.97692092, 0.97919052, 0.98280536, 0.98721005, 0.98488825, 0.98981203, 0.99263351, 0.99558788, // 1856- 1919
	0.24528362, 0.28446086, 0.32185488, 0.35514631, 0.33946117, 0.36924940, 0.38339546, 0.39228361, 0.29207947, 0.34542314, 0.39139706, 0.44077839, 0.41819279, 0.46005289, 0.48510113, 0.50001553, 0.45587048, 0.51639882, 0.55531926, 0.60930216, 0.58625084, 0.62713020, 0.65591408, 0.67155023, 0.52546169, 0.59947209, 0.63969130, 0.70633205, 0.67991472, 0.72503397, 0.76448802, 0.78297770, 0.48072988, 0.54809844, 0.58875230, 0.65071792, 0.62504980, 0.66983463, 0.70591767, 0.72425160, 0.62986682, 0.69697745, 0.72654040, 0.78498724, 0.76336376, 0.79914517, 0.83403522, 0.84881752, 0.66142930, 0.73452278, 0.76356353, 0.82748104, 0.80544151, 0.84111151, 0.88176669, 0.89626712, 0.76679249, 0.82893116, 0.84755867, 0.89794180, 0.88282046, 0.90633048, 0.93716681, 0.94597908, // 1920- 1983
	0.55287841, 0.58163142, 0.68225554, 0.72414471, 0.69721063, 0.78701157, 0.80563361, 0.86776440, 0.65598041, 0.68580935, 0.76946807, 0.80969505, 0.78481979, 0.85691104, 0.87487363, 0.92004802, 0.75439983, 0.77558057, 0.82554742, 0.85621742, 0.83769573, 0.88772606, 0.90242877, 0.93572558, 0.85315954, 0.86877979, 0.89952575, 0.92196362, 0.90910501, 0.94013917, 0.95137381, 0.97000232, 0.80659591, 0.82616710, 0.86795750, 0.89588436, 0.87952134, 0.92077348, 0.93448925, 0.95972035, 0.89495111, 0.90515893, 0.92396336, 0.94040222, 0.93112102, 0.95286189, 0.96180271, 0.97571359, 0.94300073, 0.94936380, 0.96005140, 0.97056517, 0.96489545, 0.97730403, 0.98332882, 0.99032389, 0.97318641, 0.97568944, 0.97967506, 0.98469751, 0.98205197, 0.98766568, 0.99100252, 0.99449622, // 1984- 2047
	0.23333298, 0.28296070, 0.33033737, 0.37195052, 0.35227657, 0.38956437, 0.40627133, 0.41673216, 0.29258706, 0.35966437, 0.41750345, 0.47686042, 0.44975668, 0.49999622, 0.52767514, 0.54416931, 0.49993396, 0.57491375, 0.62294217, 0.68334403, 0.65754246, 0.70321552, 0.73162768, 0.74699715, 0.58617905, 0.67188132, 0.71850305, 0.78535319, 0.75888814, 0.80403723, 0.83776815, 0.85357369, 0.53086006, 0.61196509, 0.66083720, 0.72685664, 0.69958831, 0.74724409, 0.78047871, 0.79738905, 0.71529029, 0.78570846, 0.81675412, 0.86802275, 0.84907959, 0.88039710, 0.90580269, 0.91657266, 0.74866653, 0.82087849, 0.84962791, 0.90106991, 0.88334154, 0.91201351, 0.93838594, 0.94778968, 0.85881375, 0.90909923, 0.92426661, 0.95591286, 0.94642606, 0.96122820, 0.97621154, 0.98050778, // 2048- 2111
	0.61519108, 0.64924026, 0.76863454, 0.81133994, 0.78398228, 0.87535466, 0.89041551, 0.94040641, 0.72791544, 0.76030545, 0.85114674, 0.88678769, 0.86480215, 0.92854519, 0.94057179, 0.97075354, 0.83553970, 0.85595777, 0.90409116, 0.92775208, 0.91357246, 0.95208976, 0.96054510, 0.97974293, 0.91606420, 0.92890869, 0.95427484, 0.96824081, 0.96017285, 0.97955254, 0.98454573, 0.99277351, 0.87808997, 0.89556502, 0.93285928, 0.95218722, 0.94086522, 0.96942942, 0.97624285, 0.98886104, 0.95011476, 0.95735497, 0.97073383, 0.97958130, 0.97458623, 0.98622187, 0.98964323, 0.99497079, 0.97571311, 0.97971039, 0.98641534, 0.99122508, 0.98863619, 0.99430455, 0.99616863, 0.99833488, 0.99179440, 0.99301223, 0.99495298, 0.99669304, 0.99577405, 0.99770911, 0.99849173, 0.99930018, // 2112- 2175
	0.18695654, 0.23124853, 0.27357357, 0.31321979, 0.29450671, 0.32999585, 0.34681525, 0.35739109, 0.23988643, 0.30190234, 0.35535045, 0.41421600, 0.38731890, 0.43722722, 0.46631754, 0.48365126, 0.42510642, 0.49999387, 0.54822308, 0.61302000, 0.58535060, 0.63443589, 0.66699050, 0.68468011, 0.51125844, 0.60106931, 0.64992293, 0.72598105, 0.69581316, 0.74733588, 0.78843040, 0.80770197, 0.45595570, 0.53858562, 0.58854424, 0.66127625, 0.63112430, 0.68371241, 0.72293213, 0.74288108, 0.64058838, 0.71967255, 0.75463612, 0.81753143, 0.79425476, 0.83276334, 0.86623249, 0.88042693, 0.67783856, 0.76146898, 0.79478214, 0.85994920, 0.83742054, 0.87383285, 0.91001135, 0.92291181, 0.80223488, 0.86640037, 0.88575547, 0.93032480, 0.91692218, 0.93779971, 0.96082566, 0.96742705, // 2176- 2239
	0.54348163, 0.57856022, 0.70133218, 0.75007167, 0.71871676, 0.82309032, 0.84237327, 0.90644396, 0.66492335, 0.70005458, 0.79853727, 0.84184920, 0.81509048, 0.89267525, 0.90916126, 0.95059604, 0.78083050, 0.80477137, 0.86098391, 0.89176036, 0.87327338, 0.92340393, 0.93581050, 0.96379979, 0.88148811, 0.89756688, 0.92917731, 0.94887661, 0.93758093, 0.96483959, 0.97277868, 0.98596595, 0.83405753, 0.85514262, 0.90010269, 0.92623440, 0.91093519, 0.94952118, 0.96000355, 0.97932019, 0.92406933, 0.93380633, 0.95167360, 0.96487824, 0.95744984, 0.97489245, 0.98067371, 0.98964819, 0.96180382, 0.96736891, 0.97673403, 0.98429058, 0.98021400, 0.98912753, 0.99249815, 0.99639626, 0.98549187, 0.98736567, 0.99034613, 0.99338156, 0.99178463, 0.99517681, 0.99671425, 0.99832809, // 2240- 2303
	0.15713360, 0.19801901, 0.23704120, 0.27549427, 0.25740227, 0.29175329, 0.30859080, 0.31921006, 0.20595301, 0.26473925, 0.31542564, 0.37396589, 0.34713715, 0.39684199, 0.42684998, 0.44472666, 0.37689579, 0.45176885, 0.50009905, 0.56777656, 0.53897453, 0.59032554, 0.62544228, 0.64460019, 0.46310287, 0.55564370, 0.60586237, 0.68781364, 0.65526955, 0.71081614, 0.75674734, 0.77821247, 0.40770957, 0.49141266, 0.54196999, 0.61906558, 0.58708625, 0.64288126, 0.68591192, 0.70779716, 0.59224414, 0.67716313, 0.71473459, 0.78501423, 0.75900603, 0.80207385, 0.84075091, 0.85719675, 0.63247301, 0.72327561, 0.75961250, 0.83352145, 0.80794668, 0.84931492, 0.89176360, 0.90690079, 0.76595797, 0.83881096, 0.86085802, 0.91390296, 0.89804270, 0.92266509, 0.95092146, 0.95901191, // 2304- 2367
	0.49740548, 0.53309529, 0.65802996, 0.71060131, 0.67684059, 0.78953497, 0.81142730, 0.88471089, 0.62443924, 0.66131884, 0.76474478, 0.81299896, 0.78313504, 0.86951798, 0.88893992, 0.93772542, 0.74569997, 0.77180181, 0.83345695, 0.86856910, 0.84733583, 0.90489338, 0.91987252, 0.95352809, 0.85925671, 0.87743035, 0.91304390, 0.93641985, 0.92307508, 0.95532692, 0.96522471, 0.98156476, 0.80575923, 0.82912721, 0.87906148, 0.90953361, 0.89167996, 0.93666784, 0.94954143, 0.97324051, 0.90736618, 0.91865667, 0.93941254, 0.95543950, 0.94642818, 0.96764552, 0.97490073, 0.98624026, 0.95285059, 0.95945444, 0.97050395, 0.97983432, 0.97479794, 0.98580106, 0.99013720, 0.99516126, 0.98143465, 0.98372569, 0.98739576, 0.99125708, 0.98921712, 0.99354527, 0.99557520, 0.99770184, // 2368- 2431
	0.12684809, 0.16217287, 0.19589835, 0.23077242, 0.21432509, 0.24558834, 0.26167399, 0.27180747, 0.16902288, 0.22126385, 0.26625059, 0.32147906, 0.29624233, 0.34308330, 0.37291726, 0.39070985, 0.31671118, 0.38696672, 0.43214805, 0.50002329, 0.47106931, 0.52247151, 0.55983603, 0.58012293, 0.39751386, 0.48796751, 0.53714979, 0.62361616, 0.58932064, 0.64786020, 0.69976944, 0.72409359, 0.34556217, 0.42553311, 0.47376036, 0.55276154, 0.52005683, 0.57714874, 0.62427764, 0.64822880, 0.51874298, 0.60661378, 0.64525529, 0.72375427, 0.69471589, 0.74276293, 0.78941200, 0.80916197, 0.56015115, 0.65652925, 0.69491929, 0.78059435, 0.75106071, 0.79885481, 0.85235271, 0.87142302, 0.69814104, 0.78225687, 0.80759394, 0.87483994, 0.85468994, 0.88607811, 0.92540831, 0.93664917, // 2432- 2495
	0.43280954, 0.46719001, 0.58761258, 0.64269490, 0.60728466, 0.72535232, 0.75067909, 0.83512577, 0.56098917, 0.59835893, 0.70314827, 0.75685919, 0.72365435, 0.81985837, 0.84392862, 0.90440010, 0.68335124, 0.71137871, 0.77738727, 0.81890492, 0.79386810, 0.86153821, 0.88106223, 0.92531794, 0.81427757, 0.83507966, 0.87599360, 0.90574749, 0.88870536, 0.92984166, 0.94401664, 0.96753567, 0.75256754, 0.77852017, 0.83394494, 0.87129459, 0.84940962, 0.90459741, 0.92230358, 0.95489950, 0.86968655, 0.88338624, 0.90859357, 0.93011754, 0.91796199, 0.94639005, 0.95743046, 0.97459286, 0.93111712, 0.93949622, 0.95355479, 0.96682882, 0.95966775, 0.97532579, 0.98233722, 0.99047643, 0.96971671, 0.97291080, 0.97799171, 0.98402578, 0.98085368, 0.98758812, 0.99121123, 0.99501892, // 2496- 2559
	0.13977901, 0.17746653, 0.21343598, 0.24987898, 0.23271136, 0.26526607, 0.28170664, 0.29200401, 0.18482103, 0.23980607, 0.28722806, 0.34389010, 0.31801631, 0.36601472, 0.39595155, 0.41376081, 0.34243088, 0.41467497, 0.46111267, 0.52900138, 0.50003292, 0.55137228, 0.58783911, 0.60764917, 0.42543508, 0.51680628, 0.56644607, 0.65103813, 0.61750041, 0.67472648, 0.72407709, 0.74720871, 0.37208374, 0.45366917, 0.50293858, 0.58108834, 0.54868385, 0.60519762, 0.65059840, 0.67365238, 0.55033300, 0.63666198, 0.67495468, 0.74987155, 0.72217556, 0.76811302, 0.81132170, 0.82968917, 0.59088978, 0.68505491, 0.72246659, 0.80317703, 0.77535653, 0.82038449, 0.86916716, 0.88658011, 0.72683593, 0.80634495, 0.83026391, 0.89150318, 0.87321173, 0.90168834, 0.93630429, 0.94619468, // 2560- 2623
	0.46035085, 0.49531937, 0.61768481, 0.67173258, 0.63700803, 0.75264481, 0.77668885, 0.85626606, 0.58806572, 0.62522444, 0.72942918, 0.78077808, 0.74904845, 0.84106592, 0.86314250, 0.91861815, 0.70995286, 0.73716763, 0.80125925, 0.84013971, 0.81671268, 0.88002259, 0.89759743, 0.93739515, 0.83350887, 0.85314197, 0.89185524, 0.91885424, 0.90336868, 0.94070627, 0.95308287, 0.97352070, 0.77524440, 0.80010021, 0.85318067, 0.88764232, 0.86745738, 0.91835111, 0.93393243, 0.96268700, 0.88576888, 0.89843095, 0.92171302, 0.94091685, 0.93011434, 0.95547702, 0.96491986, 0.97957699, 0.94040784, 0.94800918, 0.96078186, 0.97238197, 0.96612005, 0.97979822, 0.98566583, 0.99246779, 0.97471805, 0.97752862, 0.98201888, 0.98711017, 0.98442877, 0.99012423, 0.99307503, 0.99615939, // 2624- 2687
	0.11683770, 0.15031692, 0.18225923, 0.21601932, 0.20011276, 0.23032261, 0.24615451, 0.25609364, 0.15682286, 0.20687498, 0.25005024, 0.30414054, 0.27941704, 0.32530932, 0.35513027, 0.37289033, 0.29679921, 0.36560595, 0.40977061, 0.47753422, 0.44870387, 0.50011651, 0.53814290, 0.55882525, 0.37590701, 0.46573986, 0.51451038, 0.60241134, 0.56750006, 0.62703839, 0.68095381, 0.70621354, 0.32505075, 0.40378428, 0.45122937, 0.53084298, 0.49792543, 0.55540121, 0.60392138, 0.62856749, 0.49447326, 0.58308964, 0.62231429, 0.70346424, 0.67348836, 0.72325571, 0.77247296, 0.79335393, 0.53629178, 0.63447318, 0.67350096, 0.76314937, 0.73224263, 0.78217273, 0.83933283, 0.85970474, 0.67540641, 0.76361427, 0.79005702, 0.86194408, 0.84043540, 0.87395694, 0.91693973, 0.92922893, // 2688- 2751
	0.41147341, 0.44543611, 0.56440517, 0.62030850, 0.58428310, 0.70404888, 0.73059792, 0.81843990, 0.54003031, 0.57758095, 0.68280453, 0.73831150, 0.70397423, 0.80351997, 0.82912458, 0.89349087, 0.66281732, 0.69137704, 0.75900762, 0.80244425, 0.77613837, 0.84716720, 0.86821992, 0.91590266, 0.79943748, 0.82111300, 0.86381312, 0.89561125, 0.87736524, 0.92143872, 0.93703890, 0.96290182, 0.73497820, 0.76179884, 0.81900601, 0.85867655, 0.83544111, 0.89396273, 0.91330633, 0.94879821, 0.85728523, 0.87169529, 0.89842242, 0.92172403, 0.90860560, 0.93944083, 0.95167300, 0.97075542, 0.92394847, 0.93290474, 0.94799052, 0.96252644, 0.95466859, 0.97186802, 0.97975132, 0.98893037, 0.96584943, 0.96932833, 0.97490561, 0.98163631, 0.97808749, 0.98561699, 0.98978241, 0.99413111, // 2752- 2815
	0.10440294, 0.13491717, 0.16403697, 0.19540449, 0.18061392, 0.20868445, 0.22373597, 0.23320807, 0.14085048, 0.18698699, 0.22675168, 0.27788994, 0.25449895, 0.29786751, 0.32684574, 0.34410595, 0.26845208, 0.33300660, 0.37447526, 0.44023957, 0.41218213, 0.46185499, 0.50000667, 0.52071472, 0.34271381, 0.42877198, 0.47555210, 0.56311425, 0.52839481, 0.58766580, 0.64366648, 0.66991255, 0.29499140, 0.36957497, 0.41460523, 0.49269006, 0.46033833, 0.51674138, 0.56617075, 0.59128551, 0.45418962, 0.54115135, 0.57953184, 0.66270472, 0.63197827, 0.68288356, 0.73567909, 0.75806939, 0.49512734, 0.59300600, 0.63196301, 0.72570015, 0.69334336, 0.74568911, 0.80895583, 0.83152561, 0.63181868, 0.72332382, 0.75078706, 0.83013446, 0.80630203, 0.84328021, 0.89401398, 0.90849900, // 2816- 2879
	0.37764783, 0.41000754, 0.52331609, 0.57856045, 0.54300414, 0.66149329, 0.68897176, 0.78060682, 0.50321678, 0.53989270, 0.64277058, 0.69964821, 0.66447285, 0.76634350, 0.79426724, 0.86443266, 0.62308283, 0.65182932, 0.71965801, 0.76544438, 0.73780149, 0.81262136, 0.83621251, 0.88956648, 0.76667213, 0.78934496, 0.83397863, 0.86930459, 0.84905085, 0.89794617, 0.91665907, 0.94769223, 0.69897418, 0.72646720, 0.78513113, 0.82810839, 0.80294892, 0.86638571, 0.88888950, 0.93027238, 0.82742773, 0.84303722, 0.87178916, 0.89844698, 0.88338148, 0.91865769, 0.93380031, 0.95737650, 0.90515459, 0.91514276, 0.93193752, 0.94935493, 0.93994879, 0.96052007, 0.97089235, 0.98293196, 0.95396504, 0.95807982, 0.96464210, 0.97320796, 0.96869447, 0.97828068, 0.98411972, 0.99022980, // 2880- 2943
	0.09762694, 0.12654140, 0.15411240, 0.18418924, 0.17002291, 0.19689171, 0.21155159, 0.22073280, 0.13216470, 0.17617176, 0.21407082, 0.26359022, 0.24096696, 0.28292395, 0.31147653, 0.32847976, 0.25298065, 0.31534618, 0.35537886, 0.41992594, 0.39237796, 0.44120435, 0.47931694, 0.50002692, 0.32468361, 0.40872200, 0.45433118, 0.54173487, 0.50713190, 0.56624916, 0.62340509, 0.65018114, 0.27862873, 0.35098363, 0.39466365, 0.47192380, 0.43989753, 0.49573145, 0.54566503, 0.57103125, 0.43223518, 0.51831557, 0.55623097, 0.64051275, 0.60937839, 0.66091377, 0.71569704, 0.73888530, 0.47282643, 0.57045971, 0.60927774, 0.70532870, 0.67219536, 0.72582756, 0.79244773, 0.81621939, 0.60794768, 0.70147764, 0.72946814, 0.81278695, 0.78779362, 0.82665675, 0.88153533, 0.89720555, // 2944- 3007
	0.35924543, 0.39072935, 0.50097283, 0.55596559, 0.52060003, 0.63830909, 0.66641954, 0.75992391, 0.48319755, 0.51941836, 0.62102576, 0.67861884, 0.64299444, 0.74621026, 0.77537577, 0.84864523, 0.60148698, 0.63033361, 0.69822277, 0.74527878, 0.71691560, 0.79376849, 0.81874232, 0.87532661, 0.74887244, 0.77210164, 0.81777311, 0.85499287, 0.83365223, 0.88516089, 0.90556398, 0.93941892, 0.67939320, 0.70725741, 0.76671360, 0.81150512, 0.78527322, 0.85134044, 0.87563126, 0.92017278, 0.81122719, 0.82742876, 0.85727158, 0.88576955, 0.86967819, 0.90735642, 0.92409247, 0.95011461, 0.89493368, 0.90548553, 0.92320544, 0.94219615, 0.93194851, 0.95435411, 0.96606848, 0.97966587, 0.94751677, 0.95196421, 0.95906313, 0.96863219, 0.96358761, 0.97429185, 0.98104025, 0.98811172, // 3008- 3071
	0.18000626, 0.22353420, 0.26501267, 0.30440018, 0.28582506, 0.32110682, 0.33788523, 0.34844283, 0.23197883, 0.29320109, 0.34601163, 0.40482474, 0.37795773, 0.42777302, 0.45710593, 0.47455633, 0.41391033, 0.48875665, 0.53699902, 0.60249090, 0.57451792, 0.62410899, 0.65729359, 0.67533103, 0.49997479, 0.59043699, 0.63963994, 0.71707519, 0.68639207, 0.73881736, 0.78102002, 0.80084895, 0.44463691, 0.52764826, 0.57762890, 0.65136250, 0.62088798, 0.67416881, 0.71429668, 0.73468136, 0.62920588, 0.70981052, 0.74527029, 0.80995869, 0.78606084, 0.82558603, 0.86030578, 0.87500483, 0.66720400, 0.75253876, 0.78653393, 0.85380492, 0.83062157, 0.86809844, 0.90574250, 0.91918265, 0.79369696, 0.86000683, 0.87994387, 0.92647800, 0.91250632, 0.93426895, 0.95851498, 0.96544877, // 3072- 3135
	0.53274182, 0.56794386, 0.69129988, 0.74083622, 0.70898790, 0.81532001, 0.83507671, 0.90144210, 0.65545520, 0.69101475, 0.79064111, 0.83514335, 0.80763725, 0.88730635, 0.90451090, 0.94758589, 0.77266760, 0.79707624, 0.85460852, 0.88635690, 0.86723273, 0.91914763, 0.93206290, 0.96143223, 0.87630573, 0.89286094, 0.92540017, 0.94597107, 0.93418107, 0.96259235, 0.97101610, 0.98495214, 0.82743580, 0.84905195, 0.89515004, 0.92234377, 0.90643261, 0.94648340, 0.95754857, 0.97789548, 0.92017763, 0.93025872, 0.94883629, 0.96271121, 0.95486979, 0.97318815, 0.97933215, 0.98886734, 0.95971464, 0.96551877, 0.97526463, 0.98324946, 0.97895527, 0.98835794, 0.99194578, 0.99610769, 0.98454319, 0.98651939, 0.98965991, 0.99288612, 0.99118638, 0.99479073, 0.99644874, 0.99818123, // 3136- 3199
	0.13307854, 0.16937461, 0.20401292, 0.23942526, 0.22272646, 0.25444104, 0.27061465, 0.28078642, 0.17641943, 0.22969801, 0.27565834, 0.33127452, 0.30582515, 0.35298268, 0.38280063, 0.40056660, 0.32814190, 0.39895910, 0.44451089, 0.51203491, 0.48320857, 0.53438617, 0.57124253, 0.59133501, 0.40953047, 0.49999823, 0.54915371, 0.63456818, 0.60068483, 0.65855934, 0.70931156, 0.73309507, 0.35724367, 0.43751541, 0.48597418, 0.56436083, 0.53188254, 0.58854945, 0.63483244, 0.65835371, 0.53172490, 0.61862155, 0.65691471, 0.73387798, 0.70545845, 0.75255350, 0.79773909, 0.81688885, 0.57268492, 0.66779713, 0.70566672, 0.78920442, 0.76039113, 0.80697344, 0.85863238, 0.87702627, 0.70915855, 0.79135643, 0.81611996, 0.88092995, 0.86146228, 0.89169651, 0.92928691, 0.94002235, // 3200- 3263
	0.44455218, 0.47902214, 0.59971239, 0.65422427, 0.61922480, 0.73578584, 0.76059164, 0.84284054, 0.57209041, 0.60923593, 0.71344278, 0.76604753, 0.73350180, 0.82779811, 0.85105642, 0.90949017, 0.69382790, 0.72142967, 0.78652069, 0.82680748, 0.80242471, 0.86833518, 0.88703455, 0.92952299, 0.82157183, 0.84186663, 0.88181666, 0.91047028, 0.89402707, 0.93369559, 0.94720189, 0.96958462, 0.76134095, 0.78679313, 0.84110948, 0.87730228, 0.85609939, 0.90955399, 0.92644285, 0.95757291, 0.87559420, 0.88888547, 0.91329955, 0.93393590, 0.92229342, 0.94951402, 0.96000539, 0.97625424, 0.93448667, 0.94255432, 0.95610731, 0.96875350, 0.96192708, 0.97686050, 0.98346441, 0.99113915, 0.97145642, 0.97450995, 0.97937431, 0.98506433, 0.98206164, 0.98843155, 0.99182846, 0.99538846, // 3264- 3327
	0.10755300, 0.13994613, 0.17083812, 0.20409984, 0.18842452, 0.21821521, 0.23404920, 0.24399994, 0.14626312, 0.19517451, 0.23739508, 0.29126887, 0.26664461, 0.31236113, 0.34243761, 0.36034218, 0.28145320, 0.35005315, 0.39426470, 0.46296016, 0.43358150, 0.48552613, 0.52446343, 0.54564269, 0.36041849, 0.45085055, 0.49998820, 0.58975276, 0.55412581, 0.61495624, 0.67032233, 0.69629434, 0.30967196, 0.38852558, 0.43619943, 0.51702222, 0.48353685, 0.54194520, 0.59165410, 0.61688783, 0.47874082, 0.56912090, 0.60883316, 0.69254134, 0.66168013, 0.71295128, 0.76374570, 0.78531066, 0.52130355, 0.62181018, 0.66171951, 0.75410893, 0.72220788, 0.77379212, 0.83301584, 0.85413863, 0.66323117, 0.75414478, 0.78126570, 0.85612012, 0.83367023, 0.86852926, 0.91340571, 0.92620733, // 3328- 3391
	0.39661325, 0.43074515, 0.54999371, 0.60709596, 0.57045985, 0.69260605, 0.72009073, 0.81106882, 0.52674107, 0.56477267, 0.67138830, 0.72851337, 0.69317913, 0.79550557, 0.82202487, 0.88879417, 0.65093574, 0.68030381, 0.74944915, 0.79441763, 0.76724178, 0.84083782, 0.86261959, 0.91217329, 0.79181186, 0.81416111, 0.85814543, 0.89117036, 0.87221994, 0.91792285, 0.93424023, 0.96124890, 0.72540266, 0.75293838, 0.81169647, 0.85280479, 0.82874132, 0.88941702, 0.90954799, 0.94652515, 0.85134226, 0.86640020, 0.89398759, 0.91830687, 0.90460074, 0.93675879, 0.94947941, 0.96944586, 0.92075466, 0.93006762, 0.94567929, 0.96086608, 0.95268318, 0.97060582, 0.97885054, 0.98843638, 0.96434426, 0.96797720, 0.97377238, 0.98081769, 0.97711105, 0.98497877, 0.98931680, 0.99386887, // 3392- 3455
	0.07802708, 0.10344572, 0.12772850, 0.15556190, 0.14243059, 0.16735532, 0.18146204, 0.19032068, 0.10839674, 0.14828489, 0.18266605, 0.22997044, 0.20834847, 0.24844257, 0.27680034, 0.29368763, 0.21464975, 0.27401895, 0.31223595, 0.37643127, 0.34902441, 0.39763404, 0.43691818, 0.45826619, 0.28293535, 0.36542547, 0.41024709, 0.50001089, 0.46442076, 0.52517084, 0.58585521, 0.61429272, 0.23907241, 0.30893813, 0.35111085, 0.42915527, 0.39683457, 0.45327142, 0.50550701, 0.53206247, 0.38545145, 0.47274994, 0.51123699, 0.60017336, 0.56734089, 0.62175887, 0.68126333, 0.70650279, 0.42649226, 0.52702415, 0.56708629, 0.67000679, 0.63449827, 0.69193171, 0.76546904, 0.79169567, 0.56382722, 0.66366740, 0.69359840, 0.78571935, 0.75807196, 0.80108900, 0.86319188, 0.88095431, // 3456- 3519
	0.31794535, 0.34854129, 0.45559525, 0.51187549, 0.47571678, 0.59620133, 0.62619336, 0.72604784, 0.44215759, 0.47857189, 0.58063627, 0.64136258, 0.60380845, 0.71263892, 0.74465562, 0.82512546, 0.56073379, 0.59060292, 0.66098050, 0.71173609, 0.68112864, 0.76399731, 0.79175920, 0.85484103, 0.71852962, 0.74333987, 0.79213603, 0.83334068, 0.80972724, 0.86672252, 0.89000658, 0.92856937, 0.64415029, 0.67350446, 0.73615025, 0.78512362, 0.75643440, 0.82877325, 0.85608418, 0.90643588, 0.78528695, 0.80294031, 0.83544112, 0.86744951, 0.84939698, 0.89169699, 0.91093028, 0.94088192, 0.87951396, 0.89119810, 0.91080589, 0.93241467, 0.92075628, 0.94625643, 0.95990679, 0.97575138, 0.93868767, 0.94374064, 0.95177470, 0.96285790, 0.95701246, 0.96941360, 0.97737313, 0.98571537, // 3520- 3583
	0.08975132, 0.11792305, 0.14480421, 0.17481756, 0.16066204, 0.18755065, 0.20231828, 0.21161354, 0.12340218, 0.16687944, 0.20436503, 0.25427967, 0.23146094, 0.27375253, 0.30282802, 0.32012194, 0.24117306, 0.30419403, 0.34469109, 0.41071446, 0.38254036, 0.43251130, 0.47163045, 0.49293916, 0.31366936, 0.39932161, 0.44588244, 0.53558597, 0.50002147, 0.56076363, 0.61936021, 0.64680759, 0.26705924, 0.34050311, 0.38484923, 0.46401173, 0.43119558, 0.48841684, 0.53966906, 0.56569425, 0.42246134, 0.51090599, 0.54997436, 0.63679474, 0.60469547, 0.65785966, 0.71397871, 0.73776900, 0.46406961, 0.56460628, 0.60463968, 0.70337343, 0.66929496, 0.72437309, 0.79225439, 0.81645816, 0.60327312, 0.69956854, 0.72848257, 0.81363657, 0.78809031, 0.82783942, 0.88310175, 0.89890295, // 3584- 3647
	0.34912323, 0.38111417, 0.49309400, 0.54964704, 0.51327223, 0.63450852, 0.66341391, 0.75981899, 0.47571562, 0.51274417, 0.61661397, 0.67590919, 0.63924511, 0.74542933, 0.77531887, 0.85034138, 0.59651192, 0.62618456, 0.69604338, 0.74454635, 0.71528777, 0.79441299, 0.81987115, 0.87750581, 0.74758933, 0.77140930, 0.81832792, 0.85627880, 0.83450117, 0.88704185, 0.90753836, 0.94150434, 0.67637289, 0.70501793, 0.76609487, 0.81197866, 0.78510431, 0.85283480, 0.87729856, 0.92232735, 0.81148685, 0.82810922, 0.85866509, 0.88760026, 0.87132510, 0.90953361, 0.92624262, 0.95219803, 0.89587519, 0.90661248, 0.92464823, 0.94370472, 0.93342052, 0.95590721, 0.96741986, 0.98077741, 0.94886295, 0.95334844, 0.96049193, 0.96997822, 0.96498374, 0.97557130, 0.98210733, 0.98894570, // 3648- 3711
	0.06974016, 0.09320575, 0.11561485, 0.14195809, 0.12956669, 0.15308806, 0.16671208, 0.17528667, 0.09779373, 0.13514079, 0.16734255, 0.21279101, 0.19200595, 0.23052094, 0.25839701, 0.27500132, 0.19588728, 0.25270043, 0.28916704, 0.35217332, 0.32530530, 0.37292893, 0.41233776, 0.43374181, 0.26122349, 0.34151658, 0.38511079, 0.47482843, 0.43926111, 0.50002901, 0.56215360, 0.59131537, 0.21927124, 0.28659230, 0.32724500, 0.40451626, 0.37250324, 0.42837375, 0.48135767, 0.50826920, 0.35926709, 0.44572559, 0.48392124, 0.57425423, 0.54090659, 0.59611898, 0.65811666, 0.68440522, 0.39985215, 0.50046112, 0.54061011, 0.64638888, 0.60985044, 0.66899274, 0.74652482, 0.77419288, 0.53577640, 0.63830591, 0.66910779, 0.76599547, 0.73691026, 0.78214320, 0.84913883, 0.86826148, // 3712- 3775
	0.29589583, 0.32549165, 0.42916792, 0.48511381, 0.44913427, 0.56915503, 0.59982556, 0.70209841, 0.41845327, 0.45438280, 0.55520619, 0.61689113, 0.57871605, 0.68935668, 0.72293230, 0.80732995, 0.53538580, 0.56541641, 0.63615533, 0.68850187, 0.65688323, 0.74252936, 0.77193516, 0.83858764, 0.69795257, 0.72343962, 0.77362517, 0.81713089, 0.79219748, 0.85237793, 0.87759582, 0.91938610, 0.62134330, 0.65121038, 0.71497044, 0.76607609, 0.73617056, 0.81171608, 0.84112417, 0.89521827, 0.76673351, 0.78515008, 0.81899833, 0.85319412, 0.83392947, 0.87908475, 0.90009645, 0.93286439, 0.86794860, 0.88029397, 0.90103521, 0.92443318, 0.91180247, 0.93941019, 0.95458609, 0.97218715, 0.93149549, 0.93692887, 0.94559059, 0.95781797, 0.95138079, 0.96505667, 0.97401760, 0.98343028, // 3776- 3839
	0.05676625, 0.07639133, 0.09513291, 0.11770566, 0.10706508, 0.12727329, 0.13934482, 0.14692595, 0.08019332, 0.11192317, 0.13923999, 0.17920962, 0.16093983, 0.19480260, 0.22032665, 0.23553040, 0.16224501, 0.21158697, 0.24330035, 0.30025592, 0.27593129, 0.31905734, 0.35635732, 0.37661035, 0.21897522, 0.29071733, 0.32970757, 0.41415323, 0.38066899, 0.43786178, 0.50001251, 0.52913126, 0.18250472, 0.24177734, 0.27751898, 0.34886553, 0.31932792, 0.37088254, 0.42251623, 0.44875338, 0.30410749, 0.38355168, 0.41863395, 0.50648290, 0.47402699, 0.52780265, 0.59201002, 0.61925697, 0.34151750, 0.43595073, 0.47354419, 0.57991169, 0.54320298, 0.60256791, 0.68707754, 0.71721230, 0.46633141, 0.56765166, 0.59804582, 0.70162655, 0.67052860, 0.71888125, 0.79740192, 0.81983874, // 3840- 3903
	0.25091675, 0.27719223, 0.36917157, 0.42120378, 0.38773275, 0.49931422, 0.52965468, 0.63090155, 0.36386756, 0.39699255, 0.48985130, 0.55016664, 0.51288509, 0.62091883, 0.65622602, 0.74505752, 0.47170660, 0.50024785, 0.56752529, 0.62025632, 0.58842106, 0.67458566, 0.70660669, 0.77930696, 0.63723506, 0.66294642, 0.71356999, 0.76071616, 0.73369944, 0.79892314, 0.82905161, 0.87895865, 0.55919393, 0.58855678, 0.65123295, 0.70502393, 0.67350777, 0.75295381, 0.78679364, 0.84908381, 0.70726649, 0.72646425, 0.76180661, 0.80013166, 0.77851712, 0.82914171, 0.85513554, 0.89556356, 0.82616319, 0.83977429, 0.86263326, 0.89076326, 0.87558901, 0.90878285, 0.92925117, 0.95301425, 0.90083694, 0.90724190, 0.91745096, 0.93323375, 0.92491587, 0.94256455, 0.95567426, 0.96942394, // 3904- 3967
	0.05069201, 0.06850626, 0.08550800, 0.10634075, 0.09652092, 0.11516358, 0.12650671, 0.13363655, 0.07196441, 0.10104658, 0.12610371, 0.16348169, 0.14638809, 0.17807425, 0.20249629, 0.21703344, 0.14645469, 0.19230632, 0.22177158, 0.27593549, 0.25277693, 0.29379686, 0.33011226, 0.34984853, 0.19918035, 0.26690128, 0.30374173, 0.38573505, 0.35320225, 0.40874316, 0.47088602, 0.50001782, 0.16528671, 0.22077276, 0.25427696, 0.32279828, 0.29439697, 0.34392713, 0.39494408, 0.42085945, 0.27829700, 0.35442679, 0.38811956, 0.47472575, 0.44272701, 0.49574849, 0.56105355, 0.58873938, 0.31416299, 0.40570696, 0.44215195, 0.54875458, 0.51196917, 0.57143146, 0.65922097, 0.69053557, 0.43383851, 0.53447699, 0.56477834, 0.67145878, 0.63942565, 0.68922263, 0.77317384, 0.79715720, // 3968- 4031
	0.22982988, 0.25454820, 0.34104081, 0.39127467, 0.35898958, 0.46655738, 0.49677694, 0.59755308, 0.33829603, 0.37011105, 0.45927210, 0.51888706, 0.48203321, 0.58881980, 0.62498979, 0.71598206, 0.44183084, 0.46967936, 0.53529844, 0.58823192, 0.55627730, 0.64275972, 0.67605366, 0.75148837, 0.60878148, 0.63459746, 0.68542102, 0.73429410, 0.70628595, 0.77390513, 0.80631767, 0.86002310, 0.53007655, 0.55919253, 0.62134125, 0.67638294, 0.64413573, 0.72542469, 0.76134190, 0.82744540, 0.67941303, 0.69899158, 0.73498874, 0.77526295, 0.75255769, 0.80575316, 0.83405796, 0.87811295, 0.80658077, 0.82079918, 0.84463954, 0.87498249, 0.85861833, 0.89441809, 0.91737469, 0.94402560, 0.88646958, 0.89333069, 0.90425836, 0.92171146, 0.91250913, 0.93203369, 0.94709148, 0.96285765, // 4032- 4095
	0.21427695, 0.26172920, 0.30693298, 0.34780255, 0.32855834, 0.36510929, 0.38184716, 0.39236779, 0.27094217, 0.33595257, 0.39197391, 0.45111295, 0.42409116, 0.47422783, 0.50247365, 0.51930123, 0.46929729, 0.54417420, 0.59229433, 0.65444814, 0.62791949, 0.67492776, 0.70506211, 0.72139417, 0.55538261, 0.64281047, 0.69041916, 0.76095507, 0.73296467, 0.78076787, 0.81751185, 0.83475455, // 4120- 4127
	0.50002099, 0.58186686, 0.63121167, 0.69990326, 0.67143326, 0.72113781, 0.75684210, 0.77500732, 0.68458920, 0.75860925, 0.79116303, 0.84721213, 0.82660245, 0.86083169, 0.88953787, 0.90171647, 0.71943150, 0.79644648, 0.82710243, 0.88420756, 0.86449214, 0.89637090, 0.92672643, 0.93757216, 0.83569243, 0.89159591, 0.90842362, 0.94543434, 0.93431721, 0.95158795, 0.96989541, 0.97511601, // 4152- 4159
	0.58574610, 0.62024330, 0.74104434, 0.78618957, 0.75712811, 0.85390560, 0.87065329, 0.92660869, 0.70202770, 0.73557389, 0.82956128, 0.86834939, 0.84436412, 0.91384168, 0.92764672, 0.96244179, 0.81310077, 0.83495921, 0.88655189, 0.91297447, 0.89701752, 0.94030331, 0.95035728, 0.97322160, 0.90184625, 0.91602258, 0.94393460, 0.96028193, 0.95091574, 0.97351832, 0.97971046, 0.98996987, // 4184- 4191
	0.86000623, 0.87895280, 0.91940529, 0.94152269, 0.92856926, 0.96125982, 0.96957688, 0.98496519, 0.93940515, 0.94770047, 0.96290286, 0.97353914, 0.96755055, 0.98159302, 0.98597109, 0.99276897, 0.96999755, 0.97463875, 0.98243380, 0.98837434, 0.98517596, 0.99218103, 0.99465999, 0.99753865, 0.98920794, 0.99069085, 0.99306715, 0.99532978, 0.99413345, 0.99667138, 0.99776173, 0.99890240, // 4216- 4223
	0.16639107, 0.20763955, 0.24699586, 0.28506087, 0.26710504, 0.30115638, 0.31772092, 0.32813102, 0.21565849, 0.27433931, 0.32490319, 0.38255958, 0.35620211, 0.40505656, 0.43444534, 0.45191833, 0.38805281, 0.46141443, 0.50856730, 0.57447224, 0.54634987, 0.59622088, 0.63045548, 0.64901198, 0.47241749, 0.56249231, 0.61144877, 0.69109541, 0.65952160, 0.71340992, 0.75823176, 0.77923222, // 4248- 4255
	0.41825641, 0.50003055, 0.54936787, 0.62426290, 0.59322758, 0.64735551, 0.68930037, 0.71061876, 0.59904094, 0.68114669, 0.71743383, 0.78556034, 0.76035059, 0.80214379, 0.84008341, 0.85617952, 0.63763541, 0.72569142, 0.76077876, 0.83295274, 0.80803861, 0.84833501, 0.89038835, 0.90538494, 0.76668942, 0.83777361, 0.85911787, 0.91147908, 0.89578562, 0.92020048, 0.94879211, 0.95696426, // 4280- 4287
	0.50570744, 0.54059230, 0.66258638, 0.71345075, 0.68072108, 0.78981658, 0.81111371, 0.88221812, 0.62950299, 0.66539267, 0.76606309, 0.81292536, 0.78395895, 0.86792519, 0.88696242, 0.93489311, 0.74763675, 0.77293428, 0.83259730, 0.86696454, 0.84623124, 0.90234966, 0.91720430, 0.95075784, 0.85862207, 0.87630085, 0.91112296, 0.93422613, 0.92096488, 0.95292608, 0.96301639, 0.97971284, // 4312- 4319
	0.80629990, 0.82904032, 0.87759610, 0.90753966, 0.88999388, 0.93426128, 0.94719697, 0.97102148, 0.90558198, 0.91666230, 0.93702882, 0.95307209, 0.94402874, 0.96523769, 0.97277594, 0.98454538, 0.95137478, 0.95789680, 0.96885672, 0.97835809, 0.97323114, 0.98444316, 0.98904790, 0.99439080, 0.98013347, 0.98245861, 0.98616620, 0.99020601, 0.98807600, 0.99259860, 0.99484754, 0.99720468, // 4344- 4351
	0.13750027, 0.17498796, 0.21078244, 0.24714869, 0.23002151, 0.26255266, 0.27902999, 0.28935714, 0.18228540, 0.23713662, 0.28446121, 0.34119082, 0.31523223, 0.36334736, 0.39338012, 0.41127023, 0.33905578, 0.41153175, 0.45801541, 0.52624706, 0.49709250, 0.54876057, 0.58543140, 0.60537924, 0.42233908, 0.51405745, 0.56391505, 0.64892581, 0.61519028, 0.67278787, 0.72246149, 0.74573045, // 4376- 4383
	0.36881692, 0.45068783, 0.50007946, 0.57858509, 0.54603720, 0.60279817, 0.64852897, 0.67175368, 0.54721655, 0.63438099, 0.67278357, 0.74836579, 0.72046785, 0.76663094, 0.81022376, 0.82869740, 0.58816833, 0.68301282, 0.72071559, 0.80204900, 0.77399937, 0.81936288, 0.86847059, 0.88598150, 0.72516403, 0.80531763, 0.82928779, 0.89096474, 0.87248483, 0.90130789, 0.93605713, 0.94599499, // 4408- 4415
	0.45745849, 0.49248827, 0.61521604, 0.66961744, 0.63467846, 0.75125009, 0.77525523, 0.85546376, 0.58570427, 0.62301716, 0.72772898, 0.77948495, 0.74750619, 0.84021936, 0.86240343, 0.91832804, 0.70812642, 0.73552615, 0.80006669, 0.83914565, 0.81559423, 0.87941177, 0.89712647, 0.93717037, 0.83251320, 0.85232485, 0.89130585, 0.91848032, 0.90290257, 0.94054176, 0.95293564, 0.97351519, // 4440- 4447
	0.77389014, 0.79892381, 0.85241056, 0.88702723, 0.86670021, 0.91794843, 0.93367886, 0.96258100, 0.88516920, 0.89793456, 0.92140446, 0.94070845, 0.92984318, 0.95530771, 0.96483650, 0.97955412, 0.94012311, 0.94779304, 0.96065502, 0.97230846, 0.96602396, 0.97978382, 0.98566424, 0.99249483, 0.97465481, 0.97748826, 0.98197870, 0.98712366, 0.98441289, 0.99013844, 0.99309074, 0.99618033, // 4472- 4479
	0.10600972, 0.13709971, 0.16679068, 0.19866502, 0.18364250, 0.21218301, 0.22743247, 0.23701869, 0.14315178, 0.19012763, 0.23059993, 0.28242050, 0.25873562, 0.30267896, 0.33189634, 0.34929129, 0.27317909, 0.33878827, 0.38092485, 0.44724606, 0.41896686, 0.46913438, 0.50733083, 0.52811228, 0.34863207, 0.43568102, 0.48295484, 0.57084251, 0.53602584, 0.59552181, 0.65115011, 0.67721942, // 4504- 4511
	0.30012264, 0.37577628, 0.42139642, 0.49999833, 0.46744350, 0.52426621, 0.57357713, 0.59864609, 0.46178463, 0.54941785, 0.58804216, 0.67100013, 0.64035474, 0.69116940, 0.74325054, 0.76536016, 0.50314315, 0.60125988, 0.64033965, 0.73339963, 0.70130754, 0.75322921, 0.81527247, 0.83739117, 0.64063780, 0.73168425, 0.75896370, 0.83679945, 0.81347280, 0.84980049, 0.89884957, 0.91287653, // 4536- 4543
	0.38379580, 0.41656129, 0.53125827, 0.58688890, 0.55110855, 0.67031733, 0.69752657, 0.78853356, 0.51027027, 0.54722302, 0.65088577, 0.70752899, 0.67250085, 0.77408964, 0.80156026, 0.87061476, 0.63103155, 0.65982264, 0.72772797, 0.77312999, 0.74571552, 0.81993704, 0.84296167, 0.89518399, 0.77340767, 0.79591918, 0.84022461, 0.87484041, 0.85500251, 0.90288091, 0.92096312, 0.95091051, // 4568- 4575
	0.70628756, 0.73369105, 0.79219350, 0.83452417, 0.80974067, 0.87220636, 0.89405349, 0.93418332, 0.83366899, 0.84904866, 0.87736460, 0.90338092, 0.88871857, 0.92305716, 0.93758377, 0.96020429, 0.90909900, 0.91888012, 0.93531221, 0.95213620, 0.94306206, 0.96291830, 0.97276161, 0.98419164, 0.95646775, 0.96045857, 0.96681462, 0.97499409, 0.97067637, 0.97982367, 0.98530465, 0.99104487, // 4600- 4607
	0.11903848, 0.15279213, 0.18500248, 0.21877450, 0.20286096, 0.23304122, 0.24880831, 0.25870495, 0.15936393, 0.20961123, 0.25291005, 0.30676805, 0.28214525, 0.32780645, 0.35737510, 0.37497252, 0.30050689, 0.36889171, 0.41292834, 0.47994756, 0.45134437, 0.50212347, 0.53968403, 0.56011578, 0.37915699, 0.46815319, 0.51651964, 0.60322353, 0.56880427, 0.62750410, 0.68069409, 0.70560999, // 4632- 4639
	0.32859829, 0.40679210, 0.45395964, 0.53257128, 0.50002534, 0.55681714, 0.60464331, 0.62892511, 0.49722879, 0.58461821, 0.62314032, 0.70306107, 0.67353281, 0.72247311, 0.77099305, 0.79160732, 0.53837090, 0.63513696, 0.67363338, 0.76183864, 0.73141072, 0.78061230, 0.83731046, 0.85752501, 0.67562588, 0.76218770, 0.78815324, 0.85924448, 0.83792594, 0.87109023, 0.91426935, 0.92660323, // 4664- 4671
	0.41429836, 0.44800833, 0.56606999, 0.62113986, 0.58575044, 0.70376773, 0.72970907, 0.81622513, 0.54154388, 0.57862460, 0.68269405, 0.73734415, 0.70356870, 0.80153187, 0.82676919, 0.89036153, 0.66297417, 0.69118532, 0.75772607, 0.80051392, 0.77468364, 0.84457039, 0.86540805, 0.91253263, 0.79790567, 0.81927813, 0.86138415, 0.89292555, 0.87484228, 0.91849554, 0.93421672, 0.96027219, // 4696- 4703
	0.73429336, 0.76072196, 0.81713342, 0.85627438, 0.83332574, 0.89119787, 0.91047447, 0.94598748, 0.85501680, 0.86930390, 0.89562379, 0.91885980, 0.90575696, 0.93643534, 0.94887230, 0.96821714, 0.92195230, 0.93085339, 0.94580684, 0.96050063, 0.95257687, 0.96990170, 0.97810617, 0.98763161, 0.96400222, 0.96750934, 0.97309889, 0.98001390, 0.97637128, 0.98410623, 0.98852905, 0.99317239, // 4728- 4735
	0.09629197, 0.12541490, 0.15320163, 0.18371561, 0.16933422, 0.19664140, 0.21150989, 0.22085977, 0.13107856, 0.17559763, 0.21394198, 0.26431168, 0.24128114, 0.28395781, 0.31293083, 0.33017583, 0.25277615, 0.31637456, 0.35716040, 0.42289402, 0.39482960, 0.44463006, 0.48330373, 0.50428049, 0.32591589, 0.41149306, 0.45803606, 0.54673578, 0.51158709, 0.57167968, 0.62914814, 0.65610314, // 4760- 4767
	0.27893263, 0.35265296, 0.39717488, 0.47576795, 0.44322265, 0.50003226, 0.55045835, 0.57607362, 0.43538636, 0.52313816, 0.56194628, 0.64709847, 0.61561881, 0.66784926, 0.72257931, 0.74578785, 0.47665930, 0.57605225, 0.61553032, 0.71222742, 0.67885127, 0.73278628, 0.79886197, 0.82239161, 0.61469177, 0.70898901, 0.73726598, 0.82006019, 0.79527296, 0.83395856, 0.88739153, 0.90264042, // 4792- 4799
	0.36104426, 0.39312087, 0.50542120, 0.56128572, 0.52538263, 0.64534336, 0.67351352, 0.76780555, 0.48701222, 0.52383468, 0.62715846, 0.68537935, 0.64937952, 0.75375966, 0.78280501, 0.85582604, 0.60726179, 0.63647474, 0.70540568, 0.75278537, 0.72420336, 0.80154399, 0.82630265, 0.88218217, 0.75517041, 0.77851831, 0.82445503, 0.86137231, 0.84023721, 0.89131723, 0.91109806, 0.94394998, // 4824- 4831
	0.68542137, 0.71356225, 0.77363921, 0.81828842, 0.79213815, 0.85813447, 0.88179912, 0.92541511, 0.81776562, 0.83396903, 0.86376059, 0.89183609, 0.87599436, 0.91303750, 0.92918041, 0.95422446, 0.89952713, 0.90995713, 0.92747922, 0.94591355, 0.93597906, 0.95771629, 0.96878161, 0.98162809, 0.95086262, 0.95519794, 0.96211648, 0.97124456, 0.96643234, 0.97664754, 0.98290994, 0.98946572, // 4856- 4863
	0.08246536, 0.10795717, 0.13231432, 0.15962089, 0.14673861, 0.17118719, 0.18486828, 0.19346203, 0.11292289, 0.15241169, 0.18645227, 0.23243898, 0.21140949, 0.25038709, 0.27778609, 0.29410166, 0.21951086, 0.27709961, 0.31410452, 0.37573446, 0.34943157, 0.39609574, 0.43384665, 0.45436359, 0.28570994, 0.36517001, 0.40838672, 0.49451488, 0.46034300, 0.51868773, 0.57750750, 0.60507997, // 4888- 4895
	0.24317230, 0.31072262, 0.35148021, 0.42644655, 0.39537826, 0.44956674, 0.50001198, 0.52564243, 0.38506836, 0.46858853, 0.50544653, 0.59076726, 0.55921166, 0.61147953, 0.66951055, 0.69413451, 0.42439217, 0.52066031, 0.55899854, 0.65856157, 0.62419468, 0.67978433, 0.75282858, 0.77888555, 0.55558305, 0.65160389, 0.68042925, 0.77093357, 0.74376359, 0.78605005, 0.84951916, 0.86767096, // 4920- 4927
	0.31927822, 0.34883745, 0.45227932, 0.50609256, 0.47149484, 0.58679635, 0.61559386, 0.71162783, 0.43903948, 0.47407753, 0.57233930, 0.63078855, 0.59466244, 0.69937430, 0.73071018, 0.80944341, 0.55335510, 0.58201355, 0.64957149, 0.69857769, 0.66900076, 0.74904647, 0.77652700, 0.83882354, 0.70705442, 0.73113685, 0.77852637, 0.81929522, 0.79593088, 0.85233674, 0.87630547, 0.91602293, // 4952- 4959
	0.63459263, 0.66295003, 0.72346125, 0.77142523, 0.74332692, 0.81414829, 0.84186147, 0.89286807, 0.77208615, 0.78935320, 0.82109458, 0.85316562, 0.83509669, 0.87741530, 0.89756774, 0.92893021, 0.86877299, 0.88043226, 0.90001066, 0.92236811, 0.91030832, 0.93669221, 0.95170887, 0.96914875, 0.92949152, 0.93466919, 0.94290284, 0.95480069, 0.94853134, 0.96182991, 0.97102081, 0.98064935, // 4984- 4991
	0.07542771, 0.09909100, 0.12167318, 0.14737204, 0.13526381, 0.15826024, 0.17132938, 0.17954845, 0.10369375, 0.14064059, 0.17248054, 0.21624836, 0.19623132, 0.23334150, 0.25992686, 0.27576293, 0.20259370, 0.25715863, 0.29222357, 0.35180041, 0.32636155, 0.37147376, 0.40873237, 0.42900514, 0.26534593, 0.34165615, 0.38313744, 0.46796041, 0.43431190, 0.49177480, 0.55125409, 0.57914420, // 5016- 5023
	0.22500435, 0.28939585, 0.32826683, 0.40138115, 0.37109491, 0.42393137, 0.47437732, 0.50000383, 0.35945263, 0.44078598, 0.47684105, 0.56207788, 0.53052949, 0.58286392, 0.64254297, 0.66789454, 0.39779444, 0.49250431, 0.53025935, 0.63129348, 0.59639881, 0.65281210, 0.72943369, 0.75676329, 0.52567039, 0.62249916, 0.65155570, 0.74594303, 0.71762354, 0.76165912, 0.83028949, 0.84988344, // 5048- 5055
	0.29805308, 0.32632820, 0.42529512, 0.47801091, 0.44415125, 0.55707339, 0.58616151, 0.68299695, 0.41466179, 0.44879611, 0.54452774, 0.60306994, 0.56687902, 0.67177678, 0.70424269, 0.78588528, 0.52599165, 0.55432205, 0.62115298, 0.67102498, 0.64096566, 0.72239188, 0.75123358, 0.81673293, 0.68260116, 0.70705497, 0.75517447, 0.79791642, 0.77342960, 0.83252618, 0.85862740, 0.90185416, // 5080- 5087
	0.60877756, 0.63723280, 0.69799721, 0.74759029, 0.71853135, 0.79182453, 0.82157823, 0.87633191, 0.74888212, 0.76668528, 0.79941723, 0.83350242, 0.81428494, 0.85929634, 0.88151176, 0.91603648, 0.85315249, 0.86543070, 0.88605050, 0.91040637, 0.89726659, 0.92599684, 0.94303223, 0.96280937, 0.91863594, 0.92423368, 0.93313463, 0.94644422, 0.93943259, 0.95430388, 0.96498387, 0.97616518, // 5112- 5119
	0.10004651, 0.13433231, 0.16722257, 0.20314835, 0.18615823, 0.21827537, 0.23532245, 0.24604606, 0.14103153, 0.19352974, 0.23878243, 0.29679286, 0.27033798, 0.31939118, 0.35121164, 0.37011350, 0.28472216, 0.35946328, 0.40744530, 0.48140122, 0.44988746, 0.50546468, 0.54582093, 0.56773660, 0.37077642, 0.46833966, 0.52139817, 0.61470596, 0.57753196, 0.64082126, 0.69590080, 0.72161791, // 5144- 5151
	0.31546422, 0.40106156, 0.45270855, 0.53824624, 0.50284969, 0.56455889, 0.61497338, 0.64062178, 0.49989106, 0.59598505, 0.63802022, 0.72283597, 0.69164359, 0.74330656, 0.79207111, 0.81266730, 0.54516871, 0.65015175, 0.69211733, 0.78283474, 0.75134549, 0.80216727, 0.85682650, 0.87627025, 0.69595660, 0.78609414, 0.81316066, 0.88225317, 0.86162168, 0.89375693, 0.93196511, 0.94288135, // 5176- 5183
	0.40914651, 0.44599432, 0.57544608, 0.63510082, 0.59643630, 0.72517266, 0.75252806, 0.84217562, 0.54682200, 0.58705628, 0.69995991, 0.75771017, 0.72202617, 0.82529924, 0.85029385, 0.91294495, 0.67840978, 0.70863788, 0.78007489, 0.82419439, 0.79759510, 0.86940489, 0.88927179, 0.93411910, 0.81675251, 0.83879640, 0.88229807, 0.91258068, 0.89517370, 0.93717090, 0.95080554, 0.97321439, // 5208- 5215
	0.75148391, 0.77932711, 0.83859334, 0.87746679, 0.85474690, 0.91221083, 0.92954140, 0.96143874, 0.87530276, 0.88961031, 0.91598438, 0.93733546, 0.92528676, 0.95373649, 0.96379439, 0.97976196, 0.93568887, 0.94421109, 0.95855278, 0.97128437, 0.96442094, 0.97945183, 0.98561772, 0.99278390, 0.97365222, 0.97677854, 0.98170699, 0.98718535, 0.98430846, 0.99039459, 0.99338665, 0.99650190, // 5240- 5247
	0.06987354, 0.09674699, 0.12235176, 0.15282618, 0.13846434, 0.16573954, 0.18117896, 0.19089304, 0.10194242, 0.14499206, 0.18206286, 0.23389725, 0.21020554, 0.25418648, 0.28480919, 0.30309362, 0.21417915, 0.28032502, 0.32270980, 0.39347894, 0.36328479, 0.41691229, 0.45892435, 0.48170732, 0.29029616, 0.38144199, 0.43102037, 0.52729687, 0.48912172, 0.55424784, 0.61645613, 0.64555482, // 5272- 5279
	0.24144588, 0.31889387, 0.36565894, 0.45058017, 0.41542155, 0.47684817, 0.53144089, 0.55920655, 0.40408660, 0.50013160, 0.54232715, 0.63565175, 0.60124520, 0.65839608, 0.71738508, 0.74232369, 0.44941056, 0.55789235, 0.60122968, 0.70620661, 0.66996517, 0.72847439, 0.79830633, 0.82319393, 0.59964949, 0.70334153, 0.73434173, 0.82293937, 0.79642558, 0.83770919, 0.89229037, 0.90789306, // 5304- 5311
	0.32842888, 0.36228054, 0.48064443, 0.54172622, 0.50245900, 0.63336162, 0.66459675, 0.76790407, 0.46310714, 0.50261841, 0.61346453, 0.67704469, 0.63770773, 0.75166829, 0.78282102, 0.86114964, 0.59166636, 0.62361779, 0.69902669, 0.75056511, 0.71945530, 0.80370773, 0.82993288, 0.88915802, 0.75125478, 0.77649611, 0.82621720, 0.86537762, 0.84295445, 0.89712155, 0.91716448, 0.95036042, // 5336- 5343
	0.67603377, 0.70659882, 0.77189425, 0.81988079, 0.79178999, 0.86261805, 0.88708612, 0.93198336, 0.81871243, 0.83622406, 0.86827555, 0.89765105, 0.88107132, 0.91981914, 0.93578698, 0.96054670, 0.90241389, 0.91342002, 0.93190827, 0.95050767, 0.94046996, 0.96245406, 0.97286155, 0.98495809, 0.95497949, 0.95942235, 0.96646983, 0.97530015, 0.97065482, 0.98051444, 0.98608855, 0.99193822, // 5368- 5375
	0.05658387, 0.08014826, 0.10264639, 0.13063655, 0.11747096, 0.14251160, 0.15726757, 0.16655429, 0.08475238, 0.12353618, 0.15708535, 0.20616652, 0.18370148, 0.22536313, 0.25552548, 0.27346593, 0.18315913, 0.24536764, 0.28530706, 0.35478913, 0.32497107, 0.37767602, 0.42052230, 0.44380886, 0.25469097, 0.34311059, 0.39109075, 0.48877665, 0.45000618, 0.51601694, 0.58138034, 0.61188275, // 5400- 5407
	0.20886771, 0.28262571, 0.32741516, 0.41198361, 0.37696265, 0.43825288, 0.49457763, 0.52326661, 0.36152173, 0.45761844, 0.49974614, 0.59739248, 0.56122024, 0.62098064, 0.68444031, 0.71124015, 0.40708011, 0.51716099, 0.56119354, 0.67234788, 0.63398765, 0.69601757, 0.77248147, 0.79973999, 0.55710138, 0.66706720, 0.69933804, 0.79675009, 0.76765041, 0.81296742, 0.87472733, 0.89246941, // 5432- 5439
	0.29278978, 0.32523115, 0.43887798, 0.50031401, 0.46097270, 0.59307251, 0.62574718, 0.73462962, 0.42622164, 0.46541914, 0.57536221, 0.64144219, 0.60054947, 0.71901709, 0.75307794, 0.83863967, 0.55349308, 0.58620738, 0.66306011, 0.71817288, 0.68479660, 0.77473175, 0.80371602, 0.86954518, 0.72234810, 0.74904007, 0.80147802, 0.84455378, 0.81991916, 0.87946245, 0.90229795, 0.94036690, // 5464- 5471
	0.64273467, 0.67453909, 0.74236950, 0.79438670, 0.76397691, 0.84065215, 0.86828121, 0.91910287, 0.79381160, 0.81266180, 0.84722638, 0.87997283, 0.86156422, 0.90486454, 0.92340970, 0.95210891, 0.88771712, 0.89982035, 0.92012350, 0.94136090, 0.92990094, 0.95494902, 0.96723378, 0.98152680, 0.94672631, 0.95176590, 0.95977204, 0.97006365, 0.96463283, 0.97619369, 0.98288177, 0.98992995, // 5496- 5503
	0.03845933, 0.05587656, 0.07247755, 0.09464371, 0.08419856, 0.10398454, 0.11647636, 0.12432675, 0.05925580, 0.08915292, 0.11493799, 0.15581028, 0.13709316, 0.17177603, 0.19890009, 0.21504301, 0.13207389, 0.18251927, 0.21493386, 0.27627117, 0.25014251, 0.29645374, 0.33732550, 0.35950517, 0.19009936, 0.26611823, 0.30748571, 0.39987847, 0.36322008, 0.42572721, 0.49355960, 0.52529681, // 5528- 5535
	0.15275301, 0.21446019, 0.25162085, 0.32901183, 0.29693947, 0.35287667, 0.40926516, 0.43791877, 0.27719894, 0.36425660, 0.40267255, 0.49995797, 0.46410251, 0.52359360, 0.59361164, 0.62326944, 0.31830198, 0.42210514, 0.46347252, 0.58040281, 0.54000597, 0.60522721, 0.69505048, 0.72708438, 0.45494434, 0.56754134, 0.60116478, 0.71347996, 0.67986142, 0.73204779, 0.81280116, 0.83583745, // 5560- 5567
	0.22485983, 0.25238705, 0.34871172, 0.40583701, 0.36905537, 0.49165637, 0.52541479, 0.63793741, 0.34585551, 0.38149538, 0.48130261, 0.54777617, 0.50666454, 0.62572757, 0.66411990, 0.76072280, 0.46140792, 0.49275512, 0.56668682, 0.62492391, 0.58969178, 0.68475092, 0.71939496, 0.79765259, 0.64089969, 0.66896382, 0.72413680, 0.77464730, 0.74571522, 0.81560532, 0.84619096, 0.89699983, // 5592- 5599
	0.55629843, 0.58840890, 0.65690854, 0.71528112, 0.68109464, 0.76724038, 0.80245263, 0.86720677, 0.71690521, 0.73776166, 0.77617457, 0.81666261, 0.79384554, 0.84737440, 0.87320225, 0.91353741, 0.83768455, 0.85205105, 0.87617973, 0.90452934, 0.88924278, 0.92271465, 0.94172742, 0.96380616, 0.91354776, 0.92008835, 0.93054964, 0.94570138, 0.93770765, 0.95465658, 0.96616377, 0.97817069, // 5624- 5631
	0.04515107, 0.06482831, 0.08363463, 0.10793060, 0.09649251, 0.11820649, 0.13157087, 0.13992574, 0.06863694, 0.10188113, 0.13047167, 0.17442694, 0.15432667, 0.19156925, 0.21982120, 0.23663362, 0.15096538, 0.20573833, 0.24097092, 0.30527891, 0.27780438, 0.32650212, 0.36809969, 0.39066579, 0.21390111, 0.29460085, 0.33835372, 0.43268363, 0.39531537, 0.45914717, 0.52602640, 0.55734591, // 5656- 5663
	0.17336031, 0.23962951, 0.27960918, 0.35968473, 0.32649293, 0.38439832, 0.44081964, 0.46947718, 0.30840829, 0.39867456, 0.43854219, 0.53591893, 0.50003366, 0.55961114, 0.62716748, 0.65588100, 0.35106365, 0.45734205, 0.49964557, 0.61430576, 0.57470572, 0.63887089, 0.72368682, 0.75393328, 0.49297359, 0.60423386, 0.63771956, 0.74424926, 0.71225314, 0.76201333, 0.83568559, 0.85676895, // 5688- 5695
	0.24998948, 0.27930055, 0.38191237, 0.44089534, 0.40306312, 0.52906134, 0.56252523, 0.67360832, 0.37554460, 0.41251186, 0.51613500, 0.58238963, 0.54136678, 0.66021481, 0.69706893, 0.78954937, 0.49539103, 0.52727136, 0.60226226, 0.65926974, 0.62493832, 0.71816035, 0.75061342, 0.82421352, 0.67101730, 0.69855256, 0.75277658, 0.80049655, 0.77314771, 0.83916193, 0.86696579, 0.91299808, // 5720- 5727
	0.58821875, 0.62023555, 0.68849491, 0.74450701, 0.71169519, 0.79442067, 0.82683980, 0.88636174, 0.74525819, 0.76542932, 0.80249008, 0.84015962, 0.81892674, 0.86862932, 0.89172705, 0.92776244, 0.85618480, 0.86971590, 0.89239236, 0.91812483, 0.90426820, 0.93461716, 0.95115407, 0.97034324, 0.92581248, 0.93181187, 0.94131939, 0.95469493, 0.94764935, 0.96260562, 0.97233238, 0.98251141, // 5752- 5759
	0.03406239, 0.04994905, 0.06516113, 0.08583571, 0.07610321, 0.09461596, 0.10655647, 0.11406796, 0.05305644, 0.08080100, 0.10469466, 0.14359386, 0.12581832, 0.15874489, 0.18514261, 0.20084251, 0.11965277, 0.16721440, 0.19791095, 0.25721369, 0.23183205, 0.27676011, 0.31709417, 0.33903282, 0.17440136, 0.24748913, 0.28714330, 0.37825466, 0.34210758, 0.40383900, 0.47222277, 0.50425869, // 5784- 5791
	0.13918283, 0.19789315, 0.23326731, 0.30889679, 0.27749257, 0.33215352, 0.38852151, 0.41719469, 0.25663960, 0.34177225, 0.37923662, 0.47638823, 0.44048563, 0.50007943, 0.57152389, 0.60181014, 0.29668665, 0.39901907, 0.43977338, 0.55802794, 0.51730472, 0.58307763, 0.67625559, 0.70939679, 0.43036550, 0.54332102, 0.57730667, 0.69323762, 0.65857155, 0.71263096, 0.79773336, 0.82204919, // 5816- 5823
	0.20837084, 0.23466375, 0.32667771, 0.38282138, 0.34676076, 0.46682283, 0.50086688, 0.61423340, 0.32635114, 0.36105048, 0.45847973, 0.52492920, 0.48386399, 0.60284624, 0.64251832, 0.74187384, 0.43893776, 0.47001408, 0.54323019, 0.60218528, 0.56672912, 0.66281877, 0.69900414, 0.78046084, 0.62113874, 0.64951923, 0.70539796, 0.75769808, 0.72767014, 0.80002725, 0.83255464, 0.88649216, // 5848- 5855
	0.53524587, 0.56750927, 0.63613031, 0.69607780, 0.66096799, 0.74954093, 0.78648953, 0.85458553, 0.69817369, 0.71955541, 0.75891844, 0.80129567, 0.77741607, 0.83334504, 0.86109775, 0.90418777, 0.82552114, 0.84043537, 0.86546427, 0.89559681, 0.87935474, 0.91484939, 0.93553507, 0.95952165, 0.90548814, 0.91240077, 0.92341392, 0.93977952, 0.93116695, 0.94947135, 0.96209662, 0.97534002, // 5880- 5887
	0.02626150, 0.03890460, 0.05097719, 0.06793935, 0.05994624, 0.07513230, 0.08528852, 0.09166928, 0.04135957, 0.06384122, 0.08321874, 0.11601443, 0.10101212, 0.12879193, 0.15211611, 0.16598615, 0.09424063, 0.13380277, 0.15927723, 0.21061444, 0.18866072, 0.22755652, 0.26435944, 0.28433138, 0.13971160, 0.20228576, 0.23627677, 0.31875056, 0.28607143, 0.34190565, 0.40800632, 0.43897210, // 5912- 5919
	0.11046066, 0.15993726, 0.18976065, 0.25677818, 0.22902804, 0.27745432, 0.33051753, 0.35747611, 0.20788563, 0.28265740, 0.31555769, 0.40645309, 0.37286393, 0.42840244, 0.50001382, 0.53036104, 0.24319018, 0.33541602, 0.37214268, 0.48655137, 0.44703731, 0.51089757, 0.60901643, 0.64401799, 0.36059698, 0.46799881, 0.50025188, 0.61934766, 0.58350524, 0.63920481, 0.73552523, 0.76306509, // 5944- 5951
	0.16975955, 0.19212627, 0.27044323, 0.32041627, 0.28828901, 0.39543212, 0.42763170, 0.53493946, 0.27464917, 0.30547706, 0.39194472, 0.45440369, 0.41578507, 0.52776419, 0.56780091, 0.66861965, 0.37475850, 0.40314705, 0.47005468, 0.52725423, 0.49273305, 0.58622239, 0.62365019, 0.70864671, 0.55433265, 0.58198987, 0.63645772, 0.69117615, 0.65981962, 0.73554475, 0.77293378, 0.83495013, // 5976- 5983
	0.46968576, 0.50023558, 0.56544997, 0.62616933, 0.59061475, 0.68032431, 0.72140965, 0.79707411, 0.63031200, 0.65181656, 0.69139720, 0.73714215, 0.71137582, 0.77181272, 0.80475751, 0.85602537, 0.77556345, 0.79150153, 0.81827066, 0.85324294, 0.83436625, 0.87562588, 0.90257772, 0.93385525, 0.86680232, 0.87467640, 0.88722919, 0.90763484, 0.89687906, 0.91969454, 0.93746990, 0.95609326, // 6008- 6015
	0.02295212, 0.03420571, 0.04494551, 0.06033319, 0.05309309, 0.06685006, 0.07626511, 0.08216769, 0.03639575, 0.05664905, 0.07410774, 0.10430204, 0.09050008, 0.11610007, 0.13809719, 0.15120365, 0.08350503, 0.11960785, 0.14285986, 0.19083671, 0.17034675, 0.20668822, 0.24194657, 0.26113523, 0.12499657, 0.18307354, 0.21472973, 0.29350711, 0.26225776, 0.31562969, 0.38075723, 0.41126770, // 6040- 6047
	0.09829595, 0.14385524, 0.17131851, 0.23465950, 0.20840672, 0.25421927, 0.30586757, 0.33214314, 0.18728502, 0.25759267, 0.28875402, 0.37673013, 0.34417629, 0.39810469, 0.46967642, 0.50001297, 0.22044160, 0.30843331, 0.34338366, 0.45622963, 0.41727451, 0.48023603, 0.58048810, 0.61625277, 0.33096716, 0.43590762, 0.46724771, 0.58788688, 0.55173082, 0.60793091, 0.70910811, 0.73803829, // 6072- 6079
	0.15336543, 0.17408709, 0.24656275, 0.29397206, 0.26345117, 0.36509353, 0.39653996, 0.50145596, 0.25269923, 0.28186698, 0.36373287, 0.42447297, 0.38690457, 0.49584964, 0.53610553, 0.63759550, 0.34747499, 0.37474818, 0.43903914, 0.49545423, 0.46136263, 0.55354034, 0.59173156, 0.67819487, 0.52595363, 0.55335436, 0.60723484, 0.66296964, 0.63100655, 0.70810106, 0.74762851, 0.81306761, // 6104- 6111
	0.44181457, 0.47169512, 0.53544403, 0.59652086, 0.56075209, 0.65096855, 0.69382501, 0.77264935, 0.60148531, 0.62306734, 0.66271056, 0.70993203, 0.68333547, 0.74565871, 0.78085426, 0.83551746, 0.75438632, 0.77074287, 0.79822141, 0.83527529, 0.81529218, 0.85899045, 0.88860275, 0.92296142, 0.85038690, 0.85868412, 0.87188009, 0.89398944, 0.88233075, 0.90705804, 0.92702695, 0.94793245, // 6136- 6143
	0.08578533, 0.11662490, 0.14599485, 0.17941717, 0.16365393, 0.19354005, 0.20981730, 0.22004549, 0.12259705, 0.17065812, 0.21207704, 0.26721212, 0.24201828, 0.28866940, 0.31996257, 0.33853501, 0.25157422, 0.32220737, 0.36771891, 0.43995140, 0.40906514, 0.46370933, 0.50489631, 0.52719376, 0.33281403, 0.42741840, 0.47877513, 0.57354834, 0.53590482, 0.60002731, 0.65847530, 0.68580857, // 6168- 6175
	0.28060609, 0.36235116, 0.41169081, 0.49694461, 0.46159323, 0.52322049, 0.57564505, 0.60224412, 0.45460888, 0.55076399, 0.59297163, 0.68179629, 0.64891423, 0.70332817, 0.75682747, 0.77953905, 0.50016087, 0.60660090, 0.64893095, 0.74673317, 0.71301350, 0.76748584, 0.82921705, 0.85127521, 0.65083930, 0.74723754, 0.77601689, 0.85431193, 0.83089376, 0.86735713, 0.91326469, 0.92639012, // 6200- 6207
	0.37104144, 0.40649353, 0.53057594, 0.59103662, 0.55224565, 0.68206781, 0.71089713, 0.80754697, 0.50738734, 0.54726642, 0.65914156, 0.71964754, 0.68227785, 0.79068252, 0.81856146, 0.88866981, 0.63745065, 0.66862147, 0.74188063, 0.78955449, 0.76085355, 0.83860435, 0.86132757, 0.91300415, 0.78585633, 0.80942979, 0.85580899, 0.89034238, 0.87059536, 0.91834398, 0.93493127, 0.96243227, // 6232- 6239
	0.71594722, 0.74502381, 0.80717959, 0.85035718, 0.82510544, 0.88889259, 0.90950363, 0.94761100, 0.84858079, 0.86446341, 0.89355160, 0.91861879, 0.90444582, 0.93772675, 0.95061840, 0.97071694, 0.92002280, 0.92970630, 0.94600742, 0.96150318, 0.95312279, 0.97142689, 0.97961022, 0.98909523, 0.96486017, 0.96860303, 0.97454710, 0.98158442, 0.97787866, 0.98573186, 0.98994682, 0.99435101, // 6264- 6271
	0.05642486, 0.07925694, 0.10106622, 0.12797088, 0.11528307, 0.13934886, 0.15354482, 0.16244310, 0.08369534, 0.12112912, 0.15339667, 0.20057669, 0.17900019, 0.21900736, 0.24814481, 0.26549558, 0.17918784, 0.23846669, 0.27678950, 0.34347898, 0.31497538, 0.36549783, 0.40699081, 0.42959315, 0.24746763, 0.33220123, 0.37827642, 0.47295049, 0.43543035, 0.49952702, 0.56405235, 0.59429630, // 6296- 6303
	0.20353635, 0.27429891, 0.31702524, 0.39875584, 0.36491350, 0.42394378, 0.47934125, 0.50749039, 0.34992529, 0.44210024, 0.48282676, 0.57785904, 0.54270848, 0.60097000, 0.66460105, 0.69162159, 0.39339125, 0.50002929, 0.54238207, 0.65256849, 0.61454464, 0.67600931, 0.75434200, 0.78226996, 0.53813784, 0.64550073, 0.67783489, 0.77632164, 0.74671638, 0.79275222, 0.85843732, 0.87721756, // 6328- 6335
	0.28422226, 0.31533613, 0.42433145, 0.48381027, 0.44551561, 0.57305092, 0.60526771, 0.71259818, 0.41310693, 0.45096142, 0.55717328, 0.62194879, 0.58186803, 0.69796017, 0.73228476, 0.81855526, 0.53613097, 0.56783259, 0.64257433, 0.69701797, 0.66414786, 0.75308943, 0.78284067, 0.85025144, 0.70423035, 0.73070194, 0.78278133, 0.82676736, 0.80155075, 0.86240373, 0.88698372, 0.92767714, // 6360- 6367
	0.62499178, 0.65624364, 0.72298849, 0.77533358, 0.74460298, 0.82205406, 0.85104735, 0.90447341, 0.77536502, 0.79428914, 0.82907921, 0.86315573, 0.84393256, 0.88891631, 0.90915162, 0.94053501, 0.87486350, 0.88729819, 0.90818460, 0.93094869, 0.91866586, 0.94552619, 0.95961759, 0.97598135, 0.93735456, 0.94268942, 0.95119942, 0.96271885, 0.95663455, 0.96953053, 0.97762446, 0.98609021, // 6392- 6399
	0.04473912, 0.06438685, 0.08315749, 0.10747895, 0.09602685, 0.11776901, 0.13110579, 0.13954716, 0.06820247, 0.10141472, 0.13002856, 0.17406078, 0.15393615, 0.19122972, 0.21955185, 0.23643825, 0.15032029, 0.20527414, 0.24052314, 0.30506392, 0.27754338, 0.32645158, 0.36810372, 0.39071738, 0.21347539, 0.29429713, 0.33822242, 0.43294360, 0.39538558, 0.45948839, 0.52646708, 0.55787686, // 6424- 6431
	0.17293711, 0.23928049, 0.27935320, 0.35963266, 0.32638187, 0.38441153, 0.44101473, 0.46977313, 0.30815266, 0.39891435, 0.43891730, 0.53647495, 0.50038088, 0.56018570, 0.62791742, 0.65668643, 0.35100207, 0.45752877, 0.50000908, 0.61507274, 0.57532967, 0.63961367, 0.72453210, 0.75480564, 0.49314449, 0.60505684, 0.63838392, 0.74522125, 0.71325443, 0.76310091, 0.83661196, 0.85765499, // 6456- 6463
	0.24961727, 0.27906882, 0.38202874, 0.44107575, 0.40314413, 0.52980832, 0.56309206, 0.67516548, 0.37558559, 0.41261646, 0.51656551, 0.58306995, 0.54188448, 0.66114205, 0.69794326, 0.79066188, 0.49585110, 0.52771237, 0.60301997, 0.66012132, 0.62566008, 0.71909599, 0.75160773, 0.82537834, 0.67179068, 0.69935226, 0.75367851, 0.80144203, 0.77409886, 0.84019280, 0.86789626, 0.91385877, // 6488- 6495
	0.58878031, 0.62087387, 0.68941636, 0.74542324, 0.71259109, 0.79544040, 0.82777592, 0.88727500, 0.74619417, 0.76635385, 0.80352379, 0.84105864, 0.81980003, 0.86963435, 0.89266705, 0.92856368, 0.85687778, 0.87040923, 0.89314752, 0.91879405, 0.90495233, 0.93522344, 0.95167794, 0.97074820, 0.92640356, 0.93238637, 0.94191855, 0.95521144, 0.94820384, 0.96305063, 0.97271928, 0.98280911, // 6520- 6527
	0.02763403, 0.04089360, 0.05355046, 0.07124790, 0.06290725, 0.07873724, 0.08923385, 0.09582388, 0.04347821, 0.06695780, 0.08719534, 0.12115983, 0.10563611, 0.13441622, 0.15830536, 0.17252735, 0.09893833, 0.14005722, 0.16654167, 0.21940140, 0.19680985, 0.23687975, 0.27431305, 0.29466445, 0.14620176, 0.21079785, 0.24590007, 0.33002260, 0.29666962, 0.35360785, 0.42011040, 0.45126302, // 6552- 6559
	0.11580811, 0.16706385, 0.19799907, 0.26661178, 0.23817710, 0.28777053, 0.34145885, 0.36872946, 0.21716395, 0.29385408, 0.32762623, 0.41963892, 0.38565238, 0.44206302, 0.51347683, 0.54380174, 0.25328351, 0.34749178, 0.38491848, 0.50001275, 0.46028747, 0.52453265, 0.62159426, 0.65618583, 0.37376470, 0.48221134, 0.51482576, 0.63317504, 0.59759348, 0.65290210, 0.74705047, 0.77395002, // 6584- 6591
	0.17701069, 0.20013707, 0.28114577, 0.33222239, 0.29937723, 0.40903115, 0.44150705, 0.54993349, 0.28439759, 0.31595984, 0.40455299, 0.46774060, 0.42868104, 0.54191803, 0.58184977, 0.68229921, 0.38694105, 0.41579745, 0.48385054, 0.54138488, 0.50671086, 0.60058428, 0.63774524, 0.72200116, 0.56686824, 0.59465326, 0.64932874, 0.70356085, 0.67248396, 0.74746473, 0.78394467, 0.84437220, // 6616- 6623
	0.48202473, 0.51287429, 0.57874710, 0.63926040, 0.60380156, 0.69317793, 0.73347774, 0.80765132, 0.64298865, 0.66446338, 0.70396272, 0.74904322, 0.72363033, 0.78315887, 0.81510531, 0.86472133, 0.78481013, 0.80054268, 0.82695821, 0.86098494, 0.84262221, 0.88278816, 0.90855319, 0.93845591, 0.87387593, 0.88156839, 0.89379894, 0.91344562, 0.90309867, 0.92504553, 0.94187797, 0.95950641, // 6648- 6655
	0.03354048, 0.04899988, 0.06376051, 0.08374505, 0.07432955, 0.09217481, 0.10369317, 0.11091372, 0.05201019, 0.07884759, 0.10198500, 0.13941277, 0.12229773, 0.15402570, 0.17944853, 0.19458523, 0.11671305, 0.16254391, 0.19206982, 0.24895048, 0.22465956, 0.26776147, 0.30669427, 0.32784530, 0.16940770, 0.23965808, 0.27776632, 0.36554237, 0.33073073, 0.39014875, 0.45682825, 0.48806383, // 6680- 6687
	0.13551541, 0.19198876, 0.22602942, 0.29871729, 0.26861706, 0.32111625, 0.37582311, 0.40359287, 0.24866453, 0.33008884, 0.36603556, 0.45998747, 0.42520866, 0.48287035, 0.55297726, 0.58277497, 0.28699675, 0.38547236, 0.42461729, 0.53971292, 0.50002617, 0.56426945, 0.65710399, 0.69022091, 0.41502358, 0.52456835, 0.55746743, 0.67184614, 0.63755304, 0.69087399, 0.77796580, 0.80283903, // 6712- 6719
	0.20205456, 0.22734924, 0.31593843, 0.36980442, 0.33515263, 0.45056574, 0.48342774, 0.59315879, 0.31586516, 0.34932998, 0.44318648, 0.50753289, 0.46772715, 0.58307545, 0.62194240, 0.71969241, 0.42450590, 0.45444182, 0.52504049, 0.58240439, 0.54772641, 0.64150541, 0.67707485, 0.75750402, 0.60305613, 0.63078605, 0.68534593, 0.73733194, 0.70753940, 0.77948001, 0.81292686, 0.86834618, // 6744- 6751
	0.51887008, 0.55015137, 0.61694852, 0.67592918, 0.64131941, 0.72846062, 0.76603288, 0.83515726, 0.67860131, 0.69961488, 0.73833354, 0.78079153, 0.75683044, 0.81298499, 0.84187447, 0.88680894, 0.80969086, 0.82465010, 0.84978571, 0.88094099, 0.86413628, 0.90087742, 0.92343679, 0.94960489, 0.89200538, 0.89909680, 0.91039127, 0.92785941, 0.91865642, 0.93820861, 0.95251395, 0.96754308, // 6776- 6783
	0.02399378, 0.03589121, 0.04724317, 0.06352499, 0.05585437, 0.07042442, 0.08030250, 0.08651440, 0.03820306, 0.05961414, 0.07808050, 0.10989201, 0.09534377, 0.12231415, 0.14525953, 0.15890760, 0.08796608, 0.12614766, 0.15068601, 0.20112331, 0.17963147, 0.21781064, 0.25430505, 0.27421091, 0.13192745, 0.19301452, 0.22627836, 0.30811253, 0.27565545, 0.33107637, 0.39744188, 0.42856995, // 6808- 6815
	0.10364484, 0.15166130, 0.18068636, 0.24677930, 0.21938957, 0.26721252, 0.32026132, 0.34721155, 0.19773493, 0.27147354, 0.30393860, 0.39479482, 0.36122035, 0.41678121, 0.48914301, 0.51973746, 0.23248704, 0.32396102, 0.36046587, 0.47547747, 0.43577953, 0.50005317, 0.59966689, 0.63521995, 0.34816772, 0.45595852, 0.48830163, 0.60927290, 0.57305352, 0.62943517, 0.72797064, 0.75614995, // 6840- 6847
	0.16151604, 0.18332346, 0.25967079, 0.30903928, 0.27721880, 0.38323857, 0.41560485, 0.52346161, 0.26495561, 0.29536638, 0.38065955, 0.44319674, 0.40454052, 0.51651632, 0.55715357, 0.65923069, 0.36375217, 0.39196158, 0.45850121, 0.51605847, 0.48132685, 0.57532690, 0.61342695, 0.70010485, 0.54450165, 0.57233982, 0.62711527, 0.68269358, 0.65084064, 0.72776928, 0.76605397, 0.82954200, // 6872- 6879
	0.45929944, 0.48985836, 0.55518509, 0.61662327, 0.58062171, 0.67145166, 0.71336036, 0.79061047, 0.62098471, 0.64278905, 0.68279026, 0.72946904, 0.70311724, 0.76472181, 0.79855368, 0.85120117, 0.76945559, 0.78565799, 0.81286070, 0.84867361, 0.82934078, 0.87162371, 0.89936737, 0.93157476, 0.86269121, 0.87073524, 0.88356320, 0.90455786, 0.89350196, 0.91695330, 0.93531289, 0.95455572, // 6904- 6911
	0.01638881, 0.02480135, 0.03284078, 0.04476795, 0.03914365, 0.04981684, 0.05739190, 0.06215201, 0.02644398, 0.04194031, 0.05529742, 0.07942808, 0.06839530, 0.08885232, 0.10727648, 0.11824490, 0.06163555, 0.08999373, 0.10823755, 0.14764856, 0.13082538, 0.16065996, 0.19105220, 0.20756051, 0.09426590, 0.14138676, 0.16701786, 0.23455060, 0.20776918, 0.25349538, 0.31293442, 0.34079003, // 6936- 6943
	0.07327445, 0.10961876, 0.13155598, 0.18473972, 0.16270070, 0.20114451, 0.24718689, 0.27058250, 0.14321697, 0.20167510, 0.22749095, 0.30497861, 0.27634969, 0.32379635, 0.39100667, 0.41951697, 0.17077312, 0.24566989, 0.27551676, 0.37844429, 0.34291453, 0.40036960, 0.50001050, 0.53554779, 0.26261479, 0.35622968, 0.38437602, 0.50016739, 0.46540578, 0.51952374, 0.62651556, 0.65711485, // 6968- 6975
	0.11799281, 0.13469592, 0.19312767, 0.23311998, 0.20739888, 0.29304563, 0.32115875, 0.41470209, 0.20183653, 0.22642351, 0.29539579, 0.34933256, 0.31598209, 0.41262946, 0.45097358, 0.54728937, 0.28188446, 0.30547372, 0.36112579, 0.41251140, 0.38150032, 0.46541733, 0.50264767, 0.58710494, 0.44879843, 0.47408141, 0.52383095, 0.57863261, 0.54722360, 0.62307480, 0.66538216, 0.73558699, // 7000- 7007
	0.37010376, 0.39698829, 0.45437156, 0.51275129, 0.47855467, 0.56477814, 0.60923073, 0.69102954, 0.51944091, 0.53989914, 0.57757903, 0.62525987, 0.59836341, 0.66134873, 0.70006723, 0.76030721, 0.68579791, 0.70243051, 0.73036042, 0.77098069, 0.74906994, 0.79699207, 0.83333693, 0.87550256, 0.79030178, 0.79921094, 0.81340214, 0.83927480, 0.82564195, 0.85456252, 0.88086319, 0.90841673, // 7032- 7039
	0.01367899, 0.02084951, 0.02769853, 0.03807648, 0.03318409, 0.04247083, 0.04922133, 0.05346175, 0.02224749, 0.03563930, 0.04718382, 0.06857076, 0.05878584, 0.07692359, 0.09373159, 0.10374075, 0.05224291, 0.07709015, 0.09308583, 0.12857613, 0.11343711, 0.14030334, 0.16848710, 0.18380814, 0.08085120, 0.12297168, 0.14586718, 0.20831864, 0.18356477, 0.22583527, 0.28280552, 0.30949191, // 7064- 7071
	0.06245371, 0.09461817, 0.11403428, 0.16262416, 0.14248689, 0.17761360, 0.22113633, 0.24325367, 0.12372010, 0.17682267, 0.20022848, 0.27293913, 0.24603790, 0.29058559, 0.35599437, 0.38375933, 0.14871425, 0.21775655, 0.24522500, 0.34380320, 0.30979843, 0.36480363, 0.46447780, 0.50000701, 0.23213750, 0.32062983, 0.34715661, 0.46127773, 0.42702299, 0.48028827, 0.59034785, 0.62178616, // 7096- 7103
	0.10247282, 0.11734791, 0.16943451, 0.20601793, 0.18250769, 0.26090009, 0.28741455, 0.37596605, 0.17932131, 0.20183414, 0.26496383, 0.31587826, 0.28439608, 0.37557857, 0.41308990, 0.50739333, 0.25270193, 0.27462886, 0.32636010, 0.37556894, 0.34588837, 0.42621290, 0.46314450, 0.54682300, 0.41464772, 0.43902661, 0.48700669, 0.54152639, 0.51027137, 0.58573722, 0.62949755, 0.70206875, // 7128- 7135
	0.33830331, 0.36387120, 0.41841421, 0.47570628, 0.44217045, 0.52676063, 0.57209856, 0.65547630, 0.48319418, 0.50320748, 0.54003826, 0.58808012, 0.56101232, 0.62443145, 0.66493351, 0.72793581, 0.65596680, 0.67275241, 0.70093006, 0.74326972, 0.72044146, 0.77039767, 0.80978921, 0.85550831, 0.76448471, 0.77370344, 0.78838237, 0.81599281, 0.80144525, 0.83230949, 0.86144556, 0.89196079, // 7160- 7167
	0.03840149, 0.05757083, 0.07573690, 0.10044248, 0.08882592, 0.11095294, 0.12476647, 0.13350355, 0.06127307, 0.09443317, 0.12290929, 0.16838768, 0.14761222, 0.18612941, 0.21565545, 0.23323432, 0.14099626, 0.19781418, 0.23435613, 0.30233952, 0.27321032, 0.32449141, 0.36834043, 0.39208462, 0.20628998, 0.29083155, 0.33675492, 0.43637296, 0.39687594, 0.46415346, 0.53369496, 0.56617164, // 7192- 7199
	0.16432387, 0.23326231, 0.27495410, 0.35940059, 0.32428115, 0.38530221, 0.44446916, 0.47446598, 0.30385990, 0.40031923, 0.44228553, 0.54520634, 0.50706014, 0.56994321, 0.63952552, 0.66912638, 0.34933642, 0.46173380, 0.50678548, 0.62627008, 0.58504440, 0.65174134, 0.73737583, 0.76789540, 0.49998289, 0.61709786, 0.65255012, 0.76097239, 0.72855047, 0.77907312, 0.85092731, 0.87142084, // 7224- 7231
	0.24433061, 0.27498004, 0.38196120, 0.44434135, 0.40468256, 0.53778246, 0.57322266, 0.69097775, 0.37598034, 0.41470060, 0.52345433, 0.59309159, 0.54998544, 0.67471099, 0.71266797, 0.80760813, 0.50152918, 0.53512882, 0.61412936, 0.67372898, 0.63779022, 0.73460416, 0.76783779, 0.84294392, 0.68306414, 0.71162250, 0.76784734, 0.81623115, 0.78848871, 0.85559694, 0.88219485, 0.92654901, // 7256- 7263
	0.59744048, 0.63089468, 0.70232099, 0.75976046, 0.72600537, 0.81101759, 0.84283318, 0.90147305, 0.75977530, 0.78046902, 0.81864356, 0.85623788, 0.83506959, 0.88459805, 0.90629434, 0.94046613, 0.86776217, 0.88135907, 0.90413458, 0.92885333, 0.91555042, 0.94474826, 0.95958763, 0.97682346, 0.93549109, 0.94134166, 0.95059335, 0.96294908, 0.95642959, 0.97022953, 0.97848588, 0.98715617, // 7288- 7295
	0.02267603, 0.03536678, 0.04747788, 0.06544722, 0.05697755, 0.07302253, 0.08409957, 0.09099679, 0.03781554, 0.06117345, 0.08130663, 0.11672465, 0.10056176, 0.13056813, 0.15596637, 0.17108360, 0.09080302, 0.13361393, 0.16115452, 0.21774046, 0.19357498, 0.23640946, 0.27666070, 0.29855165, 0.14009517, 0.20858074, 0.24585033, 0.33634228, 0.30050216, 0.36170818, 0.43240540, 0.46552195, // 7320- 7327
	0.10839936, 0.16218873, 0.19479683, 0.26834375, 0.23787228, 0.29103954, 0.34840109, 0.37754354, 0.21393771, 0.29670645, 0.33307926, 0.43252972, 0.39575208, 0.45672312, 0.53205620, 0.56414195, 0.25291404, 0.35446038, 0.39509082, 0.51784451, 0.47539168, 0.54406891, 0.64381259, 0.67941280, 0.38292680, 0.50007420, 0.53523195, 0.65933557, 0.62205527, 0.67979174, 0.77342281, 0.80016489, // 7352- 7359
	0.17312311, 0.19748820, 0.28292612, 0.33800167, 0.30248154, 0.42089206, 0.45607119, 0.57313037, 0.28744866, 0.32112728, 0.41556749, 0.48350719, 0.44146904, 0.56312292, 0.60522644, 0.71098688, 0.39649686, 0.42763040, 0.50082325, 0.56247722, 0.52538903, 0.62599809, 0.66453020, 0.75253350, 0.58615452, 0.61561610, 0.67362688, 0.72976269, 0.69754676, 0.77523151, 0.81111308, 0.87063487, // 7384- 7391
	0.49671304, 0.52966066, 0.59988386, 0.66336881, 0.62619303, 0.72003735, 0.76056585, 0.83512846, 0.66635282, 0.68894809, 0.73059657, 0.77660433, 0.75065002, 0.81152978, 0.84238860, 0.89032260, 0.80563551, 0.82164186, 0.84860491, 0.88180250, 0.86388633, 0.90305641, 0.92636005, 0.95338612, 0.89307416, 0.90065779, 0.91273655, 0.93109347, 0.92141115, 0.94193320, 0.95631508, 0.97139289, // 7416- 7423
	0.01793914, 0.02868219, 0.03897451, 0.05494180, 0.04743246, 0.06169162, 0.07188694, 0.07827487, 0.03081939, 0.05119847, 0.06883446, 0.10124773, 0.08643627, 0.11393318, 0.13806584, 0.15245680, 0.07576511, 0.11435387, 0.13909342, 0.19233338, 0.16962305, 0.20993794, 0.24923473, 0.27059132, 0.12017724, 0.18393540, 0.21860539, 0.30642308, 0.27153595, 0.33105595, 0.40200115, 0.43531198, // 7448- 7455
	0.09162813, 0.14091880, 0.17071996, 0.24103159, 0.21181594, 0.26281066, 0.31956375, 0.34846344, 0.18639072, 0.26569607, 0.30031404, 0.39864477, 0.36240600, 0.42283397, 0.49988013, 0.53266578, 0.22369208, 0.32230130, 0.36154845, 0.48532785, 0.44253103, 0.51163889, 0.61580018, 0.65287218, 0.34776631, 0.46473168, 0.49981873, 0.62870338, 0.59014132, 0.64994769, 0.75021998, 0.77877299, // 7480- 7487
	0.15176528, 0.17421097, 0.25305713, 0.30616644, 0.27208007, 0.38567882, 0.42083151, 0.53778667, 0.26087214, 0.29301528, 0.38314084, 0.45056162, 0.40891756, 0.52984247, 0.57289761, 0.68185917, 0.36495926, 0.39542087, 0.46675047, 0.52931647, 0.49142204, 0.59348386, 0.63343179, 0.72555364, 0.55703865, 0.58679424, 0.64524617, 0.70373056, 0.67024805, 0.75118484, 0.78976685, 0.85400175, // 7512- 7519
	0.46647609, 0.49928652, 0.56901665, 0.63443253, 0.59620705, 0.69294690, 0.73590227, 0.81525588, 0.63826830, 0.66141908, 0.70406086, 0.75270860, 0.72532851, 0.78962103, 0.82310925, 0.87531131, 0.78695487, 0.80374938, 0.83188596, 0.86763749, 0.84836120, 0.89056296, 0.91640179, 0.94633886, 0.88034686, 0.88845174, 0.90135071, 0.92152732, 0.91089308, 0.93347449, 0.94967014, 0.96667007, // 7544- 7551
	0.00969011, 0.01608929, 0.02217987, 0.03250582, 0.02764864, 0.03688518, 0.04408309, 0.04859456, 0.01733057, 0.03017115, 0.04123265, 0.06373566, 0.05344975, 0.07250852, 0.09104104, 0.10207120, 0.04407033, 0.06968697, 0.08613737, 0.12515835, 0.10850435, 0.13806230, 0.16990191, 0.18722190, 0.07353920, 0.11912580, 0.14388866, 0.21431003, 0.18637818, 0.23408067, 0.29840092, 0.32858241, // 7576- 7583
	0.05457035, 0.08854829, 0.10905738, 0.16319895, 0.14076859, 0.17990138, 0.22908866, 0.25409459, 0.11773812, 0.17712454, 0.20331951, 0.28645299, 0.25575574, 0.30682299, 0.38072196, 0.41207232, 0.14566671, 0.22371001, 0.25476719, 0.36685019, 0.32817640, 0.39076486, 0.49984341, 0.53874911, 0.23889384, 0.34077427, 0.37120516, 0.50000616, 0.46140019, 0.52143094, 0.63810174, 0.67144935, // 7608- 7615
	0.09768306, 0.11354322, 0.16916617, 0.21028831, 0.18384315, 0.27206462, 0.30256842, 0.40428623, 0.18249366, 0.20741258, 0.27724502, 0.33516210, 0.29933630, 0.40308138, 0.44556429, 0.55223329, 0.26345357, 0.28827662, 0.34674013, 0.40298641, 0.36907551, 0.46093615, 0.50252670, 0.59651162, 0.44409798, 0.47147001, 0.52530408, 0.58571019, 0.55111677, 0.63473156, 0.68076972, 0.75709235, // 7640- 7647
	0.35894315, 0.38775148, 0.44913534, 0.51323968, 0.47567400, 0.57044075, 0.61925296, 0.70896169, 0.52056174, 0.54300927, 0.58429371, 0.63695011, 0.60728064, 0.67689377, 0.71876852, 0.78397753, 0.69720030, 0.71517684, 0.74541471, 0.78882993, 0.76541371, 0.81663884, 0.85356869, 0.89643195, 0.80815100, 0.81770015, 0.83297195, 0.85999985, 0.84574858, 0.87594378, 0.90166492, 0.92857850, // 7672- 7679
	0.01216638, 0.01986458, 0.02720607, 0.03924556, 0.03356936, 0.04434511, 0.05242240, 0.05749593, 0.02136563, 0.03647973, 0.04951496, 0.07499787, 0.06333982, 0.08491262, 0.10514587, 0.11717906, 0.05355870, 0.08306030, 0.10202481, 0.14533016, 0.12685669, 0.15958810, 0.19371246, 0.21219891, 0.08748769, 0.13852208, 0.16629848, 0.24191552, 0.21189156, 0.26315958, 0.32948348, 0.36058962, // 7704- 7711
	0.06567806, 0.10425941, 0.12749131, 0.18655242, 0.16211140, 0.20473851, 0.25625470, 0.28238830, 0.13846506, 0.20361111, 0.23244014, 0.32024567, 0.28768646, 0.34151763, 0.41647729, 0.44831343, 0.16905598, 0.25328813, 0.28677117, 0.40239628, 0.36245582, 0.42700885, 0.53461651, 0.57298838, 0.27158431, 0.37813253, 0.40990484, 0.53862733, 0.50000853, 0.56007041, 0.67172055, 0.70364698, // 7736- 7743
	0.11391332, 0.13176103, 0.19432501, 0.23903678, 0.21028970, 0.30603083, 0.33810752, 0.44467089, 0.20602939, 0.23307450, 0.30903618, 0.36980452, 0.33226144, 0.44101486, 0.48379941, 0.59111244, 0.29388639, 0.32036821, 0.38276588, 0.44088257, 0.40574326, 0.50049067, 0.54180190, 0.63487400, 0.47804560, 0.50608869, 0.56131891, 0.62111176, 0.58682822, 0.66959197, 0.71344931, 0.78619217, // 7768- 7775
	0.39122295, 0.42116926, 0.48509591, 0.54962619, 0.51181252, 0.60712975, 0.65423012, 0.74081316, 0.55587200, 0.57855481, 0.62024654, 0.67168936, 0.64267426, 0.71066970, 0.75007884, 0.81135062, 0.72413375, 0.74175649, 0.77134398, 0.81246952, 0.79031415, 0.83883895, 0.87241782, 0.91141726, 0.82980988, 0.83892722, 0.85347076, 0.87847031, 0.86531219, 0.89316879, 0.91605825, 0.93999712, // 7800- 7807
	0.00831127, 0.01398952, 0.01936942, 0.02878492, 0.02434520, 0.03272850, 0.03944666, 0.04363619, 0.01507330, 0.02666209, 0.03664918, 0.05746836, 0.04794562, 0.06562424, 0.08319948, 0.09367963, 0.03876539, 0.06227590, 0.07729274, 0.11390668, 0.09832360, 0.12604806, 0.15667496, 0.17331359, 0.06577037, 0.10829263, 0.13145994, 0.19895224, 0.17218411, 0.21790214, 0.28113072, 0.31078676, // 7832- 7839
	0.04843907, 0.07981998, 0.09880759, 0.15020627, 0.12891756, 0.16614012, 0.21400106, 0.23832848, 0.10630347, 0.16225861, 0.18724033, 0.26782512, 0.23808306, 0.28761875, 0.36088906, 0.39200376, 0.13272176, 0.20728444, 0.23700051, 0.34714629, 0.30913068, 0.37059082, 0.48053294, 0.51971047, 0.22096580, 0.31998295, 0.34986428, 0.47853243, 0.44000428, 0.49993770, 0.61944583, 0.65347598, // 7864- 7871
	0.08866052, 0.10345022, 0.15518590, 0.19434080, 0.16917081, 0.25321451, 0.28279662, 0.38245777, 0.16944507, 0.19310269, 0.25958603, 0.31589486, 0.28113147, 0.38203377, 0.42431818, 0.53054289, 0.24652152, 0.27044794, 0.32669897, 0.38205595, 0.34855968, 0.43899561, 0.48053056, 0.57524646, 0.42529641, 0.45226347, 0.50530395, 0.56605786, 0.53127369, 0.61523127, 0.66259566, 0.74104799, // 7896- 7903
	0.34108601, 0.36911777, 0.42906274, 0.49304418, 0.45557226, 0.55001003, 0.59975968, 0.69124238, 0.50095449, 0.52326881, 0.56426317, 0.61771834, 0.58763556, 0.65799695, 0.70133473, 0.76869126, 0.68221914, 0.70042848, 0.73094664, 0.77566990, 0.75161390, 0.80432542, 0.84310848, 0.88810989, 0.79606233, 0.80594091, 0.82155335, 0.84971854, 0.83487220, 0.86637782, 0.89366260, 0.92219971, // 7928- 7935
	0.00492711, 0.00843851, 0.01178472, 0.01786480, 0.01500108, 0.02043832, 0.02502195, 0.02790020, 0.00911996, 0.01651737, 0.02288687, 0.03698194, 0.03054122, 0.04248277, 0.05523864, 0.06283857, 0.02378410, 0.03917422, 0.04906679, 0.07460166, 0.06370966, 0.08304418, 0.10599366, 0.11847402, 0.04148068, 0.07072257, 0.08663083, 0.13680951, 0.11689027, 0.15090593, 0.20261440, 0.22684163, // 7960- 7967
	0.03010415, 0.05121556, 0.06395923, 0.10115515, 0.08574649, 0.11262947, 0.15049664, 0.16972646, 0.06800087, 0.10773278, 0.12524663, 0.18721385, 0.16431654, 0.20224218, 0.26449938, 0.29089606, 0.08674936, 0.14154797, 0.16336586, 0.25297585, 0.22204655, 0.27203101, 0.37351259, 0.40967941, 0.14908767, 0.22656716, 0.24988277, 0.36193496, 0.32828385, 0.38060145, 0.50000573, 0.53412702, // 7992- 7999
	0.05792783, 0.06803913, 0.10346201, 0.13177824, 0.11357859, 0.17430838, 0.19748537, 0.27498000, 0.11735551, 0.13469058, 0.18330287, 0.22737154, 0.20013933, 0.27909388, 0.31535260, 0.40654998, 0.17406788, 0.19213273, 0.23471563, 0.27934154, 0.25238390, 0.32533654, 0.36222418, 0.44601806, 0.32633406, 0.34883020, 0.39314974, 0.44802706, 0.41657772, 0.49252873, 0.54058122, 0.62026490, // 8024- 8031
	0.25455254, 0.27718397, 0.32547879, 0.38112810, 0.34854670, 0.43076598, 0.47903051, 0.56795412, 0.39075064, 0.40999149, 0.44550802, 0.49533268, 0.46721928, 0.53314460, 0.57856880, 0.64929190, 0.58161213, 0.59905200, 0.62832153, 0.67551930, 0.65005962, 0.70573422, 0.75296517, 0.80779676, 0.70149969, 0.71161988, 0.72776164, 0.75994024, 0.74298274, 0.77898187, 0.81517710, 0.85308126, // 8056- 8063
	0.00396354, 0.00685244, 0.00961526, 0.01475173, 0.01232623, 0.01691835, 0.02089967, 0.02339175, 0.00741499, 0.01361739, 0.01896206, 0.03112246, 0.02556383, 0.03587828, 0.04725158, 0.05402494, 0.01950425, 0.03258539, 0.04098769, 0.06336420, 0.05380739, 0.07075111, 0.09150874, 0.10281581, 0.03455179, 0.05997369, 0.07381398, 0.11905031, 0.10111730, 0.13175000, 0.18017351, 0.20286510, // 8088- 8095
	0.02487197, 0.04303841, 0.05399330, 0.08713247, 0.07342267, 0.09735578, 0.13234652, 0.15012596, 0.05713785, 0.09212721, 0.10757968, 0.16419127, 0.14325933, 0.17792730, 0.23696647, 0.26201296, 0.07361401, 0.12278894, 0.14240836, 0.22605768, 0.19717298, 0.24389637, 0.34291601, 0.37823265, 0.12859325, 0.19984503, 0.22117440, 0.32859922, 0.29636287, 0.34653321, 0.46592458, 0.50001736, // 8120- 8127
	0.04913852, 0.05792493, 0.08866565, 0.11387755, 0.09768470, 0.15175845, 0.17314114, 0.24444761, 0.10246535, 0.11799392, 0.16155504, 0.20207269, 0.17699240, 0.24960282, 0.28419187, 0.37107437, 0.15335835, 0.16975075, 0.20836965, 0.24999645, 0.22490232, 0.29279315, 0.32840921, 0.40906501, 0.29804137, 0.31927755, 0.36107864, 0.41431798, 0.38378714, 0.45747193, 0.50572668, 0.58573920, // 8152- 8159
	0.22984507, 0.25091916, 0.29589295, 0.34914669, 0.31797310, 0.39662692, 0.44453470, 0.53273272, 0.35925899, 0.37764272, 0.41147213, 0.46043701, 0.43283487, 0.49737170, 0.54350954, 0.61513492, 0.55285912, 0.57008359, 0.59899932, 0.64689604, 0.62106223, 0.67755959, 0.72722521, 0.78485603, 0.67446219, 0.68468558, 0.70094015, 0.73428448, 0.71671194, 0.75399380, 0.79274474, 0.83336518, // 8184- 8191
	0.16666982, 0.20725503, 0.24600736, 0.28328621, 0.26570982, 0.29906406, 0.31531607, 0.32553320, 0.21514827, 0.27277875, 0.32242638, 0.37895259, 0.35312698, 0.40098713, 0.42992060, 0.44713646, 0.38487914, 0.45654234, 0.50259318, 0.56725067, 0.53961213, 0.58856481, 0.62238131, 0.64076868, 0.46722082, 0.55550406, 0.60341196, 0.68207581, 0.65088488, 0.70414585, 0.74910665, 0.77017667, // 8216- 8223
	0.41429176, 0.49429542, 0.54257562, 0.61624910, 0.58573360, 0.63897641, 0.68074048, 0.70196703, 0.59085098, 0.67159331, 0.70726444, 0.77514291, 0.75001655, 0.79164759, 0.83026197, 0.84664286, 0.62892715, 0.71580070, 0.75043942, 0.82300888, 0.79796316, 0.83847992, 0.88201912, 0.89753238, 0.75570813, 0.82688660, 0.84824202, 0.90235074, 0.88611110, 0.91137719, 0.94209366, 0.95086939, // 8248- 8255
	0.50001603, 0.53412685, 0.65360528, 0.70365443, 0.67144101, 0.77882075, 0.80020128, 0.87136618, 0.62178654, 0.65709364, 0.75615423, 0.80286160, 0.77396465, 0.85762968, 0.87722026, 0.92637460, 0.73801486, 0.76305367, 0.82208759, 0.85672652, 0.83583027, 0.89243467, 0.90788249, 0.94290201, 0.84988692, 0.86765558, 0.90265584, 0.92660850, 0.91286820, 0.94599873, 0.95696924, 0.97512784, // 8280- 8287
	0.79714333, 0.81982747, 0.86825196, 0.89890696, 0.88095092, 0.92618076, 0.94003378, 0.96546805, 0.89721037, 0.90849903, 0.92925866, 0.94617697, 0.93665189, 0.95901346, 0.96743026, 0.98050056, 0.94597734, 0.95275161, 0.96413076, 0.97443331, 0.96887621, 0.98104364, 0.98638433, 0.99258431, 0.97660519, 0.97910183, 0.98308038, 0.98767109, 0.98525011, 0.99038393, 0.99314578, 0.99603639, // 8312- 8319
	0.14691325, 0.18483416, 0.22103248, 0.25702575, 0.24006710, 0.27225642, 0.28838484, 0.29850686, 0.19220420, 0.24703203, 0.29428734, 0.34994229, 0.32449258, 0.37167646, 0.40095114, 0.41837692, 0.35078630, 0.42144153, 0.46691730, 0.53280761, 0.50467108, 0.55460965, 0.59003115, 0.60927190, 0.43210196, 0.52097749, 0.56928631, 0.65148494, 0.61889618, 0.67453624, 0.72282586, 0.74546069, // 8344- 8351
	0.37977152, 0.45943812, 0.50750132, 0.58347437, 0.55201140, 0.60689646, 0.65119046, 0.67369431, 0.55407220, 0.63778833, 0.67474990, 0.74763083, 0.72069923, 0.76531718, 0.80788976, 0.82593457, 0.59352044, 0.68467767, 0.72094035, 0.79987866, 0.77268626, 0.81670102, 0.86531911, 0.88266816, 0.72497419, 0.80248923, 0.82577301, 0.88646360, 0.86824960, 0.89655388, 0.93196904, 0.94208227, // 8376- 8383
	0.46589421, 0.50001322, 0.61943902, 0.67174319, 0.63811993, 0.75023099, 0.77345711, 0.85095874, 0.59033942, 0.62652199, 0.72802804, 0.77798180, 0.74705707, 0.83663628, 0.85846146, 0.91330509, 0.70913805, 0.73550962, 0.79775343, 0.83573488, 0.81280484, 0.87478360, 0.89227376, 0.93198549, 0.83027175, 0.84951419, 0.88737311, 0.91427085, 0.89886382, 0.93606851, 0.94878305, 0.96989602, // 8408- 8415
	0.77316978, 0.79739627, 0.84913062, 0.88310311, 0.86320590, 0.91339773, 0.92929288, 0.95851230, 0.88153546, 0.89401129, 0.91697044, 0.93629431, 0.92539745, 0.95093755, 0.96082498, 0.97622490, 0.93716492, 0.94476093, 0.95752082, 0.96946248, 0.96302102, 0.97711362, 0.98348416, 0.99088122, 0.97210754, 0.97498350, 0.97956518, 0.98500240, 0.98213496, 0.98821627, 0.99156135, 0.99506981, // 8440- 8447
	0.07778071, 0.10633632, 0.13360443, 0.16509654, 0.15025858, 0.17845032, 0.19407395, 0.20393212, 0.11187617, 0.15689085, 0.19565245, 0.24843890, 0.22431100, 0.26901146, 0.29955110, 0.31775059, 0.23124788, 0.29866532, 0.34206872, 0.41239245, 0.38230046, 0.43558621, 0.47670729, 0.49903366, 0.30870309, 0.40020857, 0.44987486, 0.54435773, 0.50695100, 0.57089261, 0.63088521, 0.65897140, // 8472- 8479
	0.25892855, 0.33746148, 0.38487458, 0.46875815, 0.43398284, 0.49459643, 0.54773836, 0.57469126, 0.42495235, 0.51940374, 0.56122574, 0.65135200, 0.61799737, 0.67334696, 0.72961251, 0.75346823, 0.46939910, 0.57573794, 0.61798231, 0.71893003, 0.68415221, 0.74044443, 0.80687379, 0.83059121, 0.61785336, 0.71709209, 0.74687466, 0.83082951, 0.80567211, 0.84484063, 0.89654747, 0.91134177, // 8504- 8511
	0.34647554, 0.38059650, 0.50000259, 0.56010035, 0.52147973, 0.65002815, 0.68001892, 0.77938089, 0.48026388, 0.51947651, 0.62941279, 0.69094221, 0.65289968, 0.76311639, 0.79274195, 0.86746870, 0.60801074, 0.63917057, 0.71254878, 0.76206223, 0.73212305, 0.81285878, 0.83775316, 0.89362332, 0.76167006, 0.78600059, 0.83388565, 0.87107339, 0.84978307, 0.90125182, 0.92016187, 0.95160816, // 8536- 8543
	0.68921526, 0.71884216, 0.78213866, 0.82785663, 0.80104452, 0.86862474, 0.89174595, 0.93425031, 0.82665007, 0.84331281, 0.87395389, 0.90171452, 0.88607388, 0.92267098, 0.93780462, 0.96122177, 0.90632067, 0.91679532, 0.93441323, 0.95205180, 0.94253365, 0.96336062, 0.97333058, 0.98493000, 0.95635476, 0.96055229, 0.96727389, 0.97564368, 0.97122336, 0.98062368, 0.98601622, 0.99168815, // 8568- 8575
	0.06001193, 0.08392438, 0.10678998, 0.13469417, 0.12155577, 0.14646816, 0.16106029, 0.17020420, 0.08859547, 0.12756713, 0.16120254, 0.20971589, 0.18752416, 0.22863184, 0.25823727, 0.27584663, 0.18871231, 0.24996486, 0.28945161, 0.35732373, 0.32833111, 0.37982186, 0.42147684, 0.44414595, 0.25918267, 0.34581459, 0.39284155, 0.48817864, 0.45038322, 0.51488366, 0.57882879, 0.60876005, // 8600- 8607
	0.21382242, 0.28654303, 0.33040501, 0.41313508, 0.37888074, 0.43873813, 0.49394098, 0.52201450, 0.36484857, 0.45828883, 0.49955214, 0.59421405, 0.55917957, 0.61717901, 0.67962514, 0.70608611, 0.40887771, 0.51628751, 0.55894188, 0.66776250, 0.63024324, 0.69097419, 0.76689830, 0.79399277, 0.55547736, 0.66188484, 0.69391752, 0.78973377, 0.76098797, 0.80569387, 0.86824968, 0.88610631, // 8632- 8639
	0.29636411, 0.32830525, 0.43993792, 0.50001414, 0.46146634, 0.59019261, 0.62212054, 0.72851563, 0.42701917, 0.46539282, 0.57298934, 0.63750610, 0.59758940, 0.71325100, 0.74666127, 0.83089923, 0.55175193, 0.58351634, 0.65854558, 0.71226692, 0.67979036, 0.76746402, 0.79637889, 0.86159610, 0.71759499, 0.74375848, 0.79529025, 0.83790687, 0.81344320, 0.87251936, 0.89574681, 0.93431573, // 8664- 8671
	0.63940966, 0.67051955, 0.73687660, 0.78804837, 0.75810163, 0.83378403, 0.86146893, 0.91252326, 0.78780907, 0.80629787, 0.84040963, 0.87311066, 0.85470409, 0.89795240, 0.91693052, 0.94645956, 0.88281932, 0.89485525, 0.91506079, 0.93666125, 0.92501963, 0.95049502, 0.96351477, 0.97863591, 0.94250653, 0.94758266, 0.95568383, 0.96641473, 0.96076016, 0.97275030, 0.98014066, 0.98782971, // 8696- 8703
	0.07139873, 0.09834413, 0.12401696, 0.15422837, 0.13998919, 0.16701892, 0.18227011, 0.19186819, 0.10357020, 0.14641210, 0.18334441, 0.23458699, 0.21116514, 0.25460427, 0.28481431, 0.30279876, 0.21607035, 0.28129219, 0.32324506, 0.39275176, 0.36304418, 0.41577455, 0.45699068, 0.47941656, 0.29104911, 0.38077487, 0.42963534, 0.52432881, 0.48677567, 0.55088317, 0.61227978, 0.64105578, // 8728- 8735
	0.24283924, 0.31929488, 0.36533804, 0.44888131, 0.41430476, 0.47464483, 0.52850824, 0.55588362, 0.40361181, 0.49760956, 0.53915198, 0.63092414, 0.59702866, 0.65323894, 0.71175847, 0.73653527, 0.44785685, 0.55448273, 0.59679480, 0.70067378, 0.66485048, 0.72278810, 0.79259932, 0.81748889, 0.59563559, 0.69740729, 0.72783796, 0.81616925, 0.78970431, 0.83082134, 0.88644495, 0.90233306, // 8760- 8767
	0.32859624, 0.36191998, 0.47857258, 0.53860158, 0.49996374, 0.62862713, 0.65928235, 0.76134189, 0.46124934, 0.50017065, 0.60923642, 0.67183535, 0.63315909, 0.74533365, 0.77631470, 0.85433915, 0.58790187, 0.61926899, 0.69328132, 0.74424444, 0.71351336, 0.79667758, 0.82290740, 0.88222271, 0.74594531, 0.77091103, 0.82009355, 0.85925076, 0.83678117, 0.89095239, 0.91146449, 0.94544242, // 8792- 8799
	0.67142169, 0.70160377, 0.76596602, 0.81365590, 0.78573698, 0.85616799, 0.88089824, 0.92653504, 0.81280343, 0.83010514, 0.86199484, 0.89151538, 0.87488399, 0.91387666, 0.93032568, 0.95591373, 0.89793020, 0.90896583, 0.92748928, 0.94655458, 0.93627233, 0.95876045, 0.96982327, 0.98267092, 0.95139924, 0.95592023, 0.96312779, 0.97236348, 0.96748625, 0.97780966, 0.98391739, 0.99030500, // 8824- 8831
	0.03329746, 0.05032073, 0.06655091, 0.08904730, 0.07846208, 0.09863124, 0.11150952, 0.11964584, 0.05366304, 0.08360287, 0.10943160, 0.15164281, 0.13232648, 0.16807970, 0.19625053, 0.21300869, 0.12477347, 0.17691914, 0.21037481, 0.27460831, 0.24726364, 0.29577791, 0.33847525, 0.36167496, 0.18482945, 0.26409175, 0.30733158, 0.40378242, 0.36550191, 0.43081531, 0.50072763, 0.53347011, // 8856- 8863
	0.14615600, 0.21020480, 0.24885906, 0.32984943, 0.29621963, 0.35464148, 0.41321445, 0.44291707, 0.27481703, 0.36646263, 0.40667468, 0.50835600, 0.47096062, 0.53355475, 0.60471033, 0.63495645, 0.31816676, 0.42700884, 0.47030595, 0.59104276, 0.54932375, 0.61675654, 0.70695342, 0.73916345, 0.46185045, 0.57911267, 0.61453351, 0.72795574, 0.69394453, 0.74689331, 0.82578744, 0.84826339, // 8888- 8895
	0.22120963, 0.24977203, 0.34981251, 0.40988489, 0.37132738, 0.50007822, 0.53521148, 0.65202565, 0.34713987, 0.38428378, 0.48843399, 0.55731616, 0.51477967, 0.63852500, 0.67775106, 0.77623944, 0.46724023, 0.50017816, 0.57725655, 0.63765443, 0.60118819, 0.69980729, 0.73431325, 0.81340757, 0.65151828, 0.68037377, 0.73731936, 0.78816113, 0.75902475, 0.82937021, 0.85909372, 0.90842604, // 8920- 8927
	0.56467752, 0.59797868, 0.66897283, 0.72846791, 0.69363071, 0.78157732, 0.81594593, 0.87989720, 0.72953942, 0.75076135, 0.79018802, 0.83037205, 0.80770460, 0.86102384, 0.88572208, 0.92429462, 0.84753829, 0.86190603, 0.88612045, 0.91358545, 0.89873784, 0.93123207, 0.94877568, 0.96919388, 0.92173220, 0.92810722, 0.93833351, 0.95258789, 0.94507887, 0.96101398, 0.97127631, 0.98204534, // 8952- 8959
	0.02858469, 0.04367087, 0.05804743, 0.07857807, 0.06892597, 0.08727060, 0.09932750, 0.10693152, 0.04661529, 0.07364565, 0.09693470, 0.13612204, 0.11821296, 0.15140379, 0.17835487, 0.19439797, 0.10960538, 0.15763547, 0.18865819, 0.24941835, 0.22338480, 0.26954844, 0.31100418, 0.33362336, 0.16497024, 0.23944162, 0.27997621, 0.37382304, 0.33668594, 0.40007732, 0.47034409, 0.50329570, // 8984- 8991
	0.12937031, 0.18889052, 0.22478026, 0.30244910, 0.27032080, 0.32647697, 0.38443488, 0.41386585, 0.24777927, 0.33547485, 0.37406596, 0.47479178, 0.43745610, 0.49905040, 0.57232726, 0.60348887, 0.28905527, 0.39473812, 0.43686826, 0.55851637, 0.51650392, 0.58446482, 0.67888971, 0.71260536, 0.42675916, 0.54383038, 0.57918549, 0.69733442, 0.66189268, 0.71725006, 0.80249366, 0.82689134, // 9016- 9023
	0.19982523, 0.22654849, 0.32024050, 0.37791650, 0.34082191, 0.46472126, 0.50007992, 0.61709567, 0.32061921, 0.35616903, 0.45599978, 0.52459255, 0.48220838, 0.60523861, 0.64555591, 0.74720506, 0.43590041, 0.46787663, 0.54340629, 0.60423302, 0.56750905, 0.66686800, 0.70348274, 0.78622045, 0.62246641, 0.65161315, 0.70896842, 0.76217145, 0.73162537, 0.80529943, 0.83785059, 0.89161953, // 9048- 9055
	0.53453202, 0.56762376, 0.63832414, 0.69953850, 0.66372956, 0.75409228, 0.79141901, 0.85999138, 0.70135429, 0.72329655, 0.76357618, 0.80647327, 0.78227346, 0.83885274, 0.86633516, 0.90919475, 0.82892678, 0.84400965, 0.86939070, 0.89946513, 0.88324190, 0.91868614, 0.93881280, 0.96217527, 0.90901015, 0.91593747, 0.92692969, 0.94301509, 0.93453359, 0.95256225, 0.96463628, 0.97732379, // 9080- 9087
	0.01281220, 0.02150561, 0.02978601, 0.04353931, 0.03704974, 0.04944989, 0.05866084, 0.06449388, 0.02318077, 0.04044929, 0.05530449, 0.08448894, 0.07112589, 0.09587235, 0.11870367, 0.13229393, 0.05943479, 0.09340918, 0.11544350, 0.16491740, 0.14379367, 0.18118778, 0.21949885, 0.24015425, 0.09849232, 0.15724544, 0.18889415, 0.27387020, 0.24026668, 0.29774694, 0.36907392, 0.40251806, // 9112- 9119
	0.07348639, 0.11791697, 0.14468375, 0.21149954, 0.18375488, 0.23215235, 0.28831752, 0.31697372, 0.15735365, 0.23241813, 0.26480531, 0.36214290, 0.32641684, 0.38574497, 0.46499674, 0.49863098, 0.19250053, 0.28753806, 0.32513696, 0.45006485, 0.40689835, 0.47667362, 0.58524849, 0.62413167, 0.30985069, 0.42655220, 0.46197711, 0.59550469, 0.55541471, 0.61809931, 0.72497009, 0.75572332, // 9144- 9151
	0.12856077, 0.14908474, 0.22077163, 0.27145972, 0.23884881, 0.34741155, 0.38301136, 0.50009208, 0.23206176, 0.26259597, 0.34833338, 0.41503313, 0.37376593, 0.49347711, 0.53816180, 0.65062366, 0.33084795, 0.36032464, 0.43002381, 0.49318493, 0.45504067, 0.55795972, 0.60000444, 0.69567890, 0.52555765, 0.55554019, 0.61449410, 0.67552130, 0.64067822, 0.72514330, 0.76663669, 0.83559083, // 9176- 9183
	0.43379885, 0.46634318, 0.53578272, 0.60307564, 0.56374536, 0.66310532, 0.70909635, 0.79362749, 0.60799273, 0.63171923, 0.67543260, 0.72680847, 0.69793812, 0.76561848, 0.80228592, 0.85900250, 0.76679906, 0.78435863, 0.81380932, 0.85236962, 0.83161326, 0.87702467, 0.90560216, 0.93873843, 0.86656245, 0.87524908, 0.88909987, 0.91117962, 0.89949663, 0.92425749, 0.94245392, 0.96152171, // 9208- 9215
	0.10803499, 0.13856334, 0.16768962, 0.19856120, 0.18400506, 0.21161205, 0.22630075, 0.23553314, 0.14449724, 0.19022103, 0.22960403, 0.27957589, 0.25673413, 0.29909897, 0.32726125, 0.34403117, 0.27210781, 0.33508229, 0.37559208, 0.43904778, 0.41193450, 0.46002091, 0.49682550, 0.51683252, 0.34455464, 0.42795557, 0.47328732, 0.55785991, 0.52430792, 0.58157639, 0.63614438, 0.66172152, // 9240- 9247
	0.29796021, 0.37053725, 0.41432409, 0.48974707, 0.45848515, 0.51302299, 0.56098326, 0.58535827, 0.45317873, 0.53689553, 0.57388620, 0.65415276, 0.62448049, 0.67364892, 0.72537306, 0.74736209, 0.49262682, 0.58689999, 0.62446800, 0.71563448, 0.68416445, 0.73506828, 0.79818283, 0.82069460, 0.62416361, 0.71264044, 0.73915734, 0.81750609, 0.79399696, 0.83056473, 0.88266533, 0.89754804, // 9272- 9279
	0.37824850, 0.40968678, 0.51971069, 0.57299554, 0.53876244, 0.65283346, 0.67943664, 0.76793901, 0.50000360, 0.53553946, 0.63519309, 0.69023441, 0.65621214, 0.75481739, 0.78227791, 0.85128528, 0.61624383, 0.64400680, 0.70944626, 0.75394360, 0.72708291, 0.79984033, 0.82320300, 0.87627052, 0.75676022, 0.77887383, 0.82240466, 0.85752347, 0.83740276, 0.88598049, 0.90540030, 0.93756628, // 9304- 9311
	0.69052092, 0.71721799, 0.77418110, 0.81645311, 0.79169760, 0.85413061, 0.87703883, 0.91917971, 0.81621688, 0.83152527, 0.85971329, 0.88656535, 0.87141881, 0.90691548, 0.92290078, 0.94775959, 0.89626690, 0.90627289, 0.92308661, 0.94121734, 0.93142904, 0.95282995, 0.96436440, 0.97775390, 0.94653739, 0.95078133, 0.95753339, 0.96681860, 0.96192431, 0.97230398, 0.97914909, 0.98632355, // 9336- 9343
	0.09158698, 0.11913524, 0.14544393, 0.17436492, 0.16072761, 0.18660246, 0.20079910, 0.20970809, 0.12450019, 0.16667300, 0.20299408, 0.25094568, 0.22903051, 0.26966569, 0.29757988, 0.31420667, 0.23970606, 0.29995754, 0.33871827, 0.40166787, 0.37479800, 0.42247558, 0.46014105, 0.48059478, 0.30899233, 0.39078716, 0.43525788, 0.52147430, 0.48727915, 0.54563807, 0.60302877, 0.62992152, // 9368- 9375
	0.26446670, 0.33462440, 0.37697884, 0.45280250, 0.42139014, 0.47618634, 0.52594371, 0.55121848, 0.41293951, 0.49739608, 0.53461854, 0.61852918, 0.58753702, 0.63887936, 0.69455423, 0.71813740, 0.45268719, 0.54905047, 0.58740224, 0.68403455, 0.65069452, 0.70464575, 0.77358463, 0.79817697, 0.58528173, 0.67891695, 0.70699869, 0.79260678, 0.76692068, 0.80686402, 0.86532686, 0.88201618, // 9400- 9407
	0.34291620, 0.37350767, 0.48055765, 0.53462955, 0.49985843, 0.61577340, 0.64380921, 0.73741472, 0.46447187, 0.50001244, 0.59966785, 0.65711326, 0.62159942, 0.72455222, 0.75435834, 0.82924489, 0.58050394, 0.60902999, 0.67625155, 0.72367072, 0.69509444, 0.77249752, 0.79833177, 0.85679681, 0.72944452, 0.75283084, 0.79885827, 0.83730929, 0.81526920, 0.86848058, 0.89040455, 0.92673333, // 9432- 9439
	0.65922490, 0.68707401, 0.74652163, 0.79226092, 0.76547198, 0.83301441, 0.85862327, 0.90575676, 0.79245640, 0.80896245, 0.83933566, 0.86918425, 0.85236734, 0.89178120, 0.91000175, 0.93836266, 0.88176050, 0.89273226, 0.91116052, 0.93160901, 0.92057216, 0.94470425, 0.95806287, 0.97355952, 0.93785467, 0.94260991, 0.95018938, 0.96085652, 0.95523823, 0.96716380, 0.97519603, 0.98361458, // 9464- 9471
	0.04544268, 0.06468044, 0.08305986, 0.10651079, 0.09544074, 0.11643387, 0.12926226, 0.13731461, 0.06841720, 0.10062495, 0.12838544, 0.17064422, 0.15131880, 0.18714263, 0.21433823, 0.23052718, 0.14888241, 0.20144020, 0.23521575, 0.29687908, 0.27058174, 0.31728307, 0.35724558, 0.37902510, 0.20932591, 0.28662201, 0.32858535, 0.41941553, 0.38336326, 0.44483261, 0.51013838, 0.54071444, // 9496- 9503
	0.17042219, 0.23394568, 0.27224856, 0.34919162, 0.31730480, 0.37286929, 0.42766019, 0.45549478, 0.30012213, 0.38648272, 0.42465297, 0.51869032, 0.48396679, 0.54140541, 0.60809913, 0.63630943, 0.34070917, 0.44288886, 0.48349818, 0.59548819, 0.55682912, 0.61930174, 0.70464065, 0.73504964, 0.47670266, 0.58442575, 0.61673375, 0.72276989, 0.69097907, 0.74041048, 0.81670215, 0.83846914, // 9528- 9535
	0.24389021, 0.27203429, 0.37060679, 0.42707697, 0.39081232, 0.51152410, 0.54398487, 0.65175461, 0.36480602, 0.40032701, 0.50002278, 0.56420551, 0.52451856, 0.63961165, 0.67600881, 0.76754169, 0.48026991, 0.51086740, 0.58320258, 0.63882719, 0.60522091, 0.69621538, 0.72861362, 0.80232692, 0.65279872, 0.67975105, 0.73281300, 0.78063785, 0.75323563, 0.81931298, 0.84832063, 0.89636723, // 9560- 9567
	0.57147501, 0.60255569, 0.66896787, 0.72440114, 0.69192921, 0.77383540, 0.80697717, 0.86811432, 0.72581015, 0.74567886, 0.78229231, 0.82039010, 0.79888554, 0.84930749, 0.87384383, 0.91204063, 0.84107223, 0.85473571, 0.87768842, 0.90466374, 0.89010454, 0.92193202, 0.94038478, 0.96179906, 0.91350764, 0.91970384, 0.92958916, 0.94416157, 0.93647635, 0.95272753, 0.96410387, 0.97601000, // 9592- 9599
	0.03246251, 0.04747997, 0.06183517, 0.08132937, 0.07214430, 0.08960034, 0.10090739, 0.10800105, 0.05040822, 0.07656031, 0.09912137, 0.13586456, 0.11906083, 0.15021642, 0.17534064, 0.19031140, 0.11321161, 0.15815928, 0.18705048, 0.24319046, 0.21923320, 0.26174052, 0.30040352, 0.32138543, 0.16489275, 0.23398699, 0.27151956, 0.35867623, 0.32410974, 0.38311723, 0.44986330, 0.48112905, // 9624- 9631
	0.13167187, 0.18711552, 0.22052168, 0.29247647, 0.26267550, 0.31467961, 0.36923580, 0.39695048, 0.24250271, 0.32296715, 0.35853149, 0.45229460, 0.41765844, 0.47501087, 0.54561915, 0.57553479, 0.28034467, 0.37807325, 0.41697906, 0.53229798, 0.49248440, 0.55684427, 0.65068190, 0.68416261, 0.40698820, 0.51651814, 0.54949299, 0.66487521, 0.63020410, 0.68407761, 0.77265601, 0.79793571, // 9656- 9663
	0.19716629, 0.22202797, 0.30911581, 0.36253495, 0.32820991, 0.44256010, 0.47539597, 0.58500843, 0.30979383, 0.34290614, 0.43577158, 0.50000422, 0.46030936, 0.57541364, 0.61457554, 0.71294060, 0.41727078, 0.44702687, 0.51722698, 0.57470784, 0.54004510, 0.63395887, 0.66993349, 0.75142024, 0.59639714, 0.62419461, 0.67887295, 0.73139360, 0.70128050, 0.77396945, 0.80803544, 0.86448956, // 9688- 9695
	0.51195254, 0.54319644, 0.60987279, 0.66930378, 0.63448774, 0.72224192, 0.76038993, 0.83054454, 0.67220374, 0.69335043, 0.73230774, 0.77531298, 0.75104016, 0.80794716, 0.83747780, 0.88333742, 0.80541942, 0.82056390, 0.84599126, 0.87770765, 0.86059386, 0.89803339, 0.92114988, 0.94799380, 0.88908633, 0.89631266, 0.90778894, 0.92567516, 0.91624911, 0.93622546, 0.95099463, 0.96645926, // 9720- 9727
	0.04048620, 0.05811309, 0.07493819, 0.09690155, 0.08654439, 0.10617813, 0.11843490, 0.12612042, 0.06153778, 0.09144274, 0.11721492, 0.15737234, 0.13900078, 0.17304199, 0.19945632, 0.21519130, 0.13526410, 0.18491804, 0.21683516, 0.27636580, 0.25097148, 0.29606762, 0.33551732, 0.35700234, 0.19237699, 0.26651084, 0.30684970, 0.39620754, 0.36077799, 0.42128092, 0.48714142, 0.51798280, // 9752- 9759
	0.15565535, 0.21607516, 0.25249986, 0.32754297, 0.29646240, 0.35068819, 0.40536565, 0.43315697, 0.27809207, 0.36230043, 0.39939223, 0.49333555, 0.45863580, 0.51610111, 0.58421829, 0.61311129, 0.31771461, 0.41813309, 0.45809566, 0.57135930, 0.53225242, 0.59543151, 0.68403976, 0.71561048, 0.44999655, 0.55853925, 0.59103514, 0.70066020, 0.66779626, 0.71893089, 0.79988405, 0.82300267, // 9784- 9791
	0.22606431, 0.25296220, 0.34710928, 0.40238314, 0.36687193, 0.48524989, 0.51784791, 0.62636249, 0.34381906, 0.37841160, 0.47549266, 0.53970514, 0.50001470, 0.61513394, 0.65259986, 0.74671220, 0.45623386, 0.48656443, 0.55805820, 0.61440715, 0.58035193, 0.67241901, 0.70615715, 0.78286008, 0.63127875, 0.65854643, 0.71223008, 0.76183335, 0.73340743, 0.80205163, 0.83294847, 0.88420320, // 9816- 9823
	0.54873947, 0.57989686, 0.64645973, 0.70335644, 0.67000157, 0.75410452, 0.78919404, 0.85381393, 0.70532257, 0.72568624, 0.76315399, 0.80319490, 0.78061818, 0.83356935, 0.85997645, 0.90110124, 0.82747008, 0.84169623, 0.86559794, 0.89437162, 0.87884854, 0.91281020, 0.93304083, 0.95652608, 0.90417363, 0.91076769, 0.92126244, 0.93709586, 0.92875471, 0.94643304, 0.95910101, 0.97236436, // 9848- 9855
	0.01720952, 0.02729247, 0.03689175, 0.05179335, 0.04478184, 0.05811303, 0.06760427, 0.07357961, 0.02925592, 0.04833229, 0.06475279, 0.09504710, 0.08120300, 0.10683932, 0.12957028, 0.14310055, 0.07136198, 0.10732114, 0.13041455, 0.18016537, 0.15894386, 0.19649871, 0.23365567, 0.25374776, 0.11265431, 0.17220633, 0.20452119, 0.28738458, 0.25458572, 0.31057537, 0.37910927, 0.41121167, // 9880- 9887
	0.08615607, 0.13210260, 0.15981121, 0.22593061, 0.19853029, 0.24630946, 0.30064587, 0.32827279, 0.17459349, 0.24842856, 0.28113558, 0.37439744, 0.33981904, 0.39713363, 0.47232806, 0.50424065, 0.20951443, 0.30202278, 0.33893124, 0.45817342, 0.41693198, 0.48345776, 0.58738394, 0.62441874, 0.32557168, 0.43686588, 0.47047803, 0.59685613, 0.55897166, 0.61798976, 0.72095318, 0.75036365, // 9912- 9919
	0.14234017, 0.16335281, 0.23693565, 0.28671798, 0.25477115, 0.36144258, 0.39501567, 0.50687046, 0.24522154, 0.27544229, 0.36046674, 0.42464443, 0.38492416, 0.50005082, 0.54247228, 0.64907805, 0.34341327, 0.37216491, 0.43963310, 0.49958390, 0.46357094, 0.56105803, 0.60117400, 0.69165228, 0.53023805, 0.55896960, 0.61557206, 0.67363474, 0.64037928, 0.72072422, 0.76072037, 0.82714211, // 9944- 9951
	0.44212856, 0.47355776, 0.54055414, 0.60459301, 0.56709337, 0.66173702, 0.70566926, 0.78650362, 0.60928650, 0.63198946, 0.67346382, 0.72247016, 0.69489036, 0.75938412, 0.79478975, 0.84973066, 0.76357132, 0.78043103, 0.80877509, 0.84608010, 0.82594749, 0.86997555, 0.89860183, 0.93179311, 0.86049994, 0.86885049, 0.88222769, 0.90398907, 0.89251582, 0.91682868, 0.93561042, 0.95526392, // 9976- 9983
	0.01390635, 0.02237712, 0.03046648, 0.04335771, 0.03727163, 0.04880981, 0.05731498, 0.06264957, 0.02403530, 0.04037304, 0.05447518, 0.08132654, 0.06904785, 0.09181897, 0.11270535, 0.12513928, 0.05943939, 0.09086378, 0.11109817, 0.15607121, 0.13684612, 0.17094904, 0.20576672, 0.22466467, 0.09556289, 0.14898722, 0.17799299, 0.25539734, 0.22468253, 0.27709695, 0.34377523, 0.37502836, //10008-10015
	0.07236583, 0.11303251, 0.13755005, 0.19846913, 0.17324140, 0.21724060, 0.26932021, 0.29577009, 0.14969438, 0.21716631, 0.24694948, 0.33585085, 0.30300870, 0.35746233, 0.43217174, 0.46390287, 0.18153612, 0.26772729, 0.30207917, 0.41813371, 0.37805392, 0.44290138, 0.54904931, 0.58691416, 0.28748824, 0.39477357, 0.42704153, 0.55446768, 0.51622319, 0.57565402, 0.68469395, 0.71582769, //10040-10047
	0.12280830, 0.14156208, 0.20727650, 0.25327642, 0.22368652, 0.32242270, 0.35447258, 0.46191221, 0.21773968, 0.24566644, 0.32399444, 0.38542838, 0.34742326, 0.45755686, 0.50002527, 0.60668446, 0.30838364, 0.33541588, 0.39907853, 0.45726219, 0.42213309, 0.51732639, 0.55789774, 0.65006948, 0.49255770, 0.52067710, 0.57607135, 0.63512683, 0.60128619, 0.68301491, 0.72568361, 0.79646076, //10072-10079
	0.40573521, 0.43596369, 0.50048411, 0.56459982, 0.52705107, 0.62177361, 0.66784736, 0.75263659, 0.57043121, 0.59300268, 0.63458369, 0.68506179, 0.65653028, 0.72336214, 0.76149432, 0.82083085, 0.73450390, 0.75186149, 0.78103783, 0.82101379, 0.79942315, 0.84658007, 0.87887355, 0.91630379, 0.83754642, 0.84646603, 0.86067266, 0.88472632, 0.87203676, 0.89895331, 0.92074214, 0.94357155, //10104-10111
	0.00564533, 0.01005405, 0.01425498, 0.02212789, 0.01840152, 0.02545359, 0.03138820, 0.03513897, 0.01090931, 0.02039236, 0.02857410, 0.04684552, 0.03851128, 0.05402595, 0.07028791, 0.07997901, 0.02931748, 0.04938299, 0.06227490, 0.09556904, 0.08132939, 0.10651223, 0.13560183, 0.15135377, 0.05241330, 0.09047142, 0.11115039, 0.17491364, 0.14965529, 0.19274282, 0.25494088, 0.28410949, //10136-10143
	0.03756295, 0.06508013, 0.08166473, 0.12940413, 0.10964602, 0.14417263, 0.19055339, 0.21416776, 0.08703948, 0.13875087, 0.16152247, 0.23920179, 0.21052552, 0.25801175, 0.33137039, 0.36252654, 0.11149647, 0.18154676, 0.20932566, 0.31772606, 0.28041718, 0.34084704, 0.45271661, 0.49265713, 0.19236435, 0.28901741, 0.31815538, 0.44788414, 0.40886803, 0.46948432, 0.59347513, 0.62894871, //10168-10175
	0.07359873, 0.08673691, 0.13266714, 0.16917547, 0.14569531, 0.22372973, 0.25293589, 0.34932341, 0.14872599, 0.17074499, 0.23247853, 0.28695668, 0.25331192, 0.35088136, 0.39345273, 0.50002162, 0.22046617, 0.24318805, 0.29666203, 0.35110852, 0.31820235, 0.40694748, 0.44935766, 0.54524367, 0.39771612, 0.42439638, 0.47678374, 0.53835634, 0.50308138, 0.58816975, 0.63770621, 0.71950081, //10200-10207
	0.31417175, 0.34151848, 0.39992979, 0.46410482, 0.42646414, 0.52119686, 0.57270425, 0.66719380, 0.47271552, 0.49517381, 0.53612296, 0.59094977, 0.56000096, 0.63228851, 0.67787700, 0.74856459, 0.66145157, 0.68003813, 0.71130123, 0.75801165, 0.73281320, 0.78800234, 0.82932298, 0.87738843, 0.77999228, 0.79017251, 0.80641757, 0.83631825, 0.82058763, 0.85399037, 0.88337947, 0.91420535, //10232-10239
	0.05206400, 0.07297903, 0.09291348, 0.11766677, 0.10601221, 0.12811889, 0.14132785, 0.14961306, 0.07705146, 0.11141036, 0.14099962, 0.18470923, 0.16473733, 0.20177721, 0.22925231, 0.24562217, 0.16451249, 0.21916531, 0.25433660, 0.31668093, 0.29003195, 0.33721788, 0.37693626, 0.39857011, 0.22736298, 0.30617684, 0.34899633, 0.43929129, 0.40349575, 0.46453424, 0.52830632, 0.55818384, //10264-10271
	0.18693162, 0.25237077, 0.29187876, 0.36898772, 0.33706162, 0.39282029, 0.44666540, 0.47402975, 0.32174586, 0.40833134, 0.44654429, 0.53863745, 0.50459976, 0.56102189, 0.62526166, 0.65252692, 0.36250095, 0.46390508, 0.50422767, 0.61310333, 0.57550988, 0.63626366, 0.71814622, 0.74733105, 0.49846823, 0.60346191, 0.63502657, 0.73656311, 0.70604919, 0.75348952, 0.82593555, 0.84663922, //10296-10303
	0.26199462, 0.29090084, 0.39203034, 0.44827128, 0.41207670, 0.53255483, 0.56415050, 0.66892879, 0.38377325, 0.41952295, 0.51978639, 0.58278205, 0.54381568, 0.65665317, 0.69162610, 0.77959297, 0.50000966, 0.53036524, 0.60186944, 0.65584715, 0.62323993, 0.71145913, 0.74228780, 0.81273690, 0.66788521, 0.69413209, 0.74582503, 0.79158145, 0.76533184, 0.82869030, 0.85617516, 0.90166359, //10328-10335
	0.58873012, 0.61925155, 0.68435675, 0.73775068, 0.70649242, 0.78530585, 0.81690577, 0.87499129, 0.73888425, 0.75805293, 0.79334111, 0.82968576, 0.80916523, 0.85719666, 0.88042079, 0.91656421, 0.84881503, 0.86190877, 0.88392098, 0.90949692, 0.89569211, 0.92591109, 0.94334278, 0.96359623, 0.91782718, 0.92373103, 0.93313203, 0.94690375, 0.93965000, 0.95504100, 0.96578958, 0.97705720, //10360-10367
	0.04390842, 0.06253025, 0.08031399, 0.10313526, 0.09237736, 0.11276879, 0.12531734, 0.13320242, 0.06615798, 0.09742937, 0.12437263, 0.16562542, 0.14676422, 0.18175958, 0.20849857, 0.22443491, 0.14403079, 0.19527647, 0.22817050, 0.28864397, 0.26288599, 0.30867214, 0.34819546, 0.36972457, 0.20295031, 0.27860596, 0.31975827, 0.40941144, 0.37383456, 0.43456194, 0.49979211, 0.53034682, //10392-10399
	0.16508870, 0.22708220, 0.26448391, 0.34020051, 0.30883399, 0.36356517, 0.41801379, 0.44568238, 0.29135644, 0.37642568, 0.41386959, 0.50729223, 0.47281172, 0.52998108, 0.59688425, 0.62528569, 0.33142460, 0.43219914, 0.47231425, 0.58421068, 0.54562313, 0.60806875, 0.69455606, 0.72537520, 0.46489783, 0.57246354, 0.60461781, 0.71173309, 0.67960535, 0.72963891, 0.80788627, 0.83026071, //10424-10431
	0.23693742, 0.26447617, 0.36088661, 0.41644248, 0.38071245, 0.49997549, 0.53208225, 0.63951871, 0.35601911, 0.39098468, 0.48911855, 0.55296814, 0.51348803, 0.62796220, 0.66464726, 0.75685976, 0.46964487, 0.49998752, 0.57152903, 0.62713748, 0.59361964, 0.68441284, 0.71734897, 0.79202614, 0.64253165, 0.66951381, 0.72259484, 0.77099846, 0.74323924, 0.81019667, 0.84006260, 0.88952356, //10456-10463
	0.56102911, 0.59201273, 0.65810339, 0.71395847, 0.68124286, 0.76374298, 0.79772464, 0.86027286, 0.71568864, 0.73564260, 0.77245463, 0.81134496, 0.78940472, 0.84078313, 0.86622923, 0.90573389, 0.83401055, 0.84789652, 0.87119732, 0.89898035, 0.88400386, 0.91678614, 0.93615505, 0.95862523, 0.90832969, 0.91471129, 0.92487725, 0.94005258, 0.93206332, 0.94902715, 0.96109290, 0.97373431, //10488-10495
	0.02467419, 0.03790196, 0.05053780, 0.06885753, 0.06022161, 0.07656131, 0.08760684, 0.09448535, 0.04046821, 0.06448051, 0.08513315, 0.12062961, 0.10440343, 0.13448703, 0.15955568, 0.17445943, 0.09580515, 0.13892126, 0.16669265, 0.22267193, 0.19872307, 0.24108132, 0.28044666, 0.30180674, 0.14540229, 0.21353389, 0.25054415, 0.33903316, 0.30398319, 0.36381900, 0.43255675, 0.46474416, //10520-10527
	0.11346813, 0.16738372, 0.19992619, 0.27231723, 0.24238193, 0.29458457, 0.35048510, 0.37884494, 0.21986133, 0.30104712, 0.33693121, 0.43332320, 0.39768975, 0.45672215, 0.52999422, 0.56099731, 0.25831455, 0.35746487, 0.39692716, 0.51608003, 0.47504883, 0.54145517, 0.63888511, 0.67364568, 0.38568013, 0.49908192, 0.53305057, 0.65330113, 0.61727717, 0.67320733, 0.76537968, 0.79167357, //10552-10559
	0.17793825, 0.20229097, 0.28734745, 0.34143445, 0.30670974, 0.42248005, 0.45660677, 0.57002817, 0.29056509, 0.32377909, 0.41679316, 0.48284421, 0.44201239, 0.56017291, 0.60099595, 0.70331613, 0.39810498, 0.42847201, 0.50020000, 0.55956364, 0.52361951, 0.62107580, 0.65837362, 0.74338023, 0.58284287, 0.61149443, 0.66789612, 0.72244869, 0.69115187, 0.76667877, 0.80218661, 0.86081888, //10584-10591
	0.49573103, 0.52778961, 0.59616707, 0.65787005, 0.62173473, 0.71292327, 0.75246252, 0.82558492, 0.66097836, 0.68287416, 0.72329826, 0.76806273, 0.74280931, 0.80195581, 0.83278385, 0.88052425, 0.79914660, 0.81483045, 0.84126518, 0.87421012, 0.85642216, 0.89530406, 0.91920949, 0.94692565, 0.88591991, 0.89342510, 0.90538871, 0.92390512, 0.91414036, 0.93483122, 0.95003270, 0.96592351, //10616-10623
	0.01748633, 0.02765839, 0.03737681, 0.05234862, 0.04529869, 0.05868100, 0.06819253, 0.07418765, 0.02964585, 0.04884130, 0.06538927, 0.09571076, 0.08185804, 0.10760157, 0.13029204, 0.14382677, 0.07226274, 0.10823319, 0.13137534, 0.18112473, 0.15991800, 0.19759435, 0.23458131, 0.25463928, 0.11365403, 0.17322680, 0.20558545, 0.28829304, 0.25544480, 0.31144250, 0.37979676, 0.41176616, //10648-10655
	0.08704771, 0.13301641, 0.16082225, 0.22686961, 0.19948713, 0.24725555, 0.30144849, 0.32896972, 0.17592520, 0.24944559, 0.28184572, 0.37509189, 0.34062803, 0.39775988, 0.47273608, 0.50454419, 0.21057717, 0.30307581, 0.33981710, 0.45864470, 0.41759067, 0.48392547, 0.58752023, 0.62445576, 0.32631058, 0.43762171, 0.47080925, 0.59698472, 0.55918577, 0.61798554, 0.72067961, 0.75003149, //10680-10687
	0.14324751, 0.16428991, 0.23802116, 0.28776358, 0.25577946, 0.36230108, 0.39580645, 0.50718782, 0.24604203, 0.27635726, 0.36119414, 0.42528584, 0.38565325, 0.50049059, 0.54267247, 0.64904359, 0.34419699, 0.37285534, 0.44044815, 0.50005713, 0.46406046, 0.56139828, 0.60127837, 0.69130496, 0.53053367, 0.55921522, 0.61563293, 0.67352340, 0.64037427, 0.72038776, 0.76037618, 0.82652751, //10712-10719
	0.44265709, 0.47400136, 0.54088313, 0.60472071, 0.56732671, 0.66167334, 0.70547282, 0.78605903, 0.60936900, 0.63196009, 0.67352279, 0.72217709, 0.69475024, 0.75902956, 0.79423437, 0.84911044, 0.76333872, 0.78017728, 0.80840823, 0.84567671, 0.82557769, 0.86949762, 0.89813380, 0.93132875, 0.86006710, 0.86844926, 0.88180609, 0.90351325, 0.89205914, 0.91635170, 0.93515444, 0.95484384, //10744-10751
	0.02182705, 0.03384548, 0.04532257, 0.06231383, 0.05429920, 0.06950600, 0.07989968, 0.08645465, 0.03618836, 0.05827215, 0.07731085, 0.11077429, 0.09546698, 0.12385487, 0.14795223, 0.16230996, 0.08645523, 0.12680262, 0.15264985, 0.20611588, 0.18330527, 0.22381478, 0.26223524, 0.28310431, 0.13281427, 0.19756071, 0.23273123, 0.31889487, 0.28473426, 0.34305130, 0.41162052, 0.44373261, //10776-10783
	0.10303423, 0.15381639, 0.18444394, 0.25426635, 0.22535631, 0.27584216, 0.33104777, 0.35909236, 0.20233339, 0.28066104, 0.31505482, 0.41024923, 0.37507033, 0.43342266, 0.50723278, 0.53865662, 0.23920447, 0.33588712, 0.37430691, 0.49333370, 0.45230002, 0.51871974, 0.61853319, 0.65412391, 0.36205762, 0.47473762, 0.50839945, 0.63097107, 0.59415476, 0.65135645, 0.74763310, 0.77511955, //10808-10815
	0.16417716, 0.18721717, 0.26779281, 0.32014340, 0.28650897, 0.39863674, 0.43252491, 0.54500769, 0.27294584, 0.30495037, 0.39481773, 0.45995118, 0.41965027, 0.53648738, 0.57785087, 0.68178308, 0.37669002, 0.40640637, 0.47638671, 0.53605992, 0.49998546, 0.59733881, 0.63585924, 0.72274165, 0.56209517, 0.59074144, 0.64719398, 0.70301999, 0.67100531, 0.74835370, 0.78554663, 0.84726175, //10840-10847
	0.47470114, 0.50644104, 0.57423639, 0.63681652, 0.60015133, 0.69258064, 0.73384255, 0.80991919, 0.64048899, 0.66268416, 0.70352588, 0.74993923, 0.72374646, 0.78499878, 0.81748564, 0.86796079, 0.78494935, 0.80110454, 0.82821600, 0.86287675, 0.84419328, 0.88507992, 0.91083862, 0.94073701, 0.87566456, 0.88351289, 0.89603771, 0.91583272, 0.90537217, 0.92751657, 0.94412386, 0.96152982, //10872-10879
	0.01007951, 0.01712532, 0.02386178, 0.03537178, 0.02994207, 0.04023316, 0.04824013, 0.05325235, 0.01848092, 0.03275599, 0.04507831, 0.07008218, 0.05862775, 0.07982267, 0.10017679, 0.11224604, 0.04791773, 0.07663264, 0.09514220, 0.13842287, 0.11998256, 0.15269103, 0.18733964, 0.20630985, 0.08097517, 0.13173032, 0.15927913, 0.23604491, 0.20562375, 0.25758057, 0.32549596, 0.35725922, //10904-10911
	0.05972404, 0.09761049, 0.12054978, 0.18011693, 0.15540802, 0.19841722, 0.25099203, 0.27766233, 0.13055474, 0.19624399, 0.22537711, 0.31516392, 0.28209090, 0.33700218, 0.41387782, 0.44645538, 0.16146615, 0.24680951, 0.28087468, 0.39941968, 0.35843627, 0.42465220, 0.53458472, 0.57379094, 0.26502144, 0.37386297, 0.40713824, 0.53905431, 0.49926064, 0.56103417, 0.67469200, 0.70719890, //10936-10943
	0.10763342, 0.12523874, 0.18710531, 0.23241131, 0.20325997, 0.30067628, 0.33306575, 0.44207456, 0.20025287, 0.22752698, 0.30392325, 0.36595250, 0.32766885, 0.43874171, 0.48269253, 0.59307084, 0.28876025, 0.31553632, 0.37893839, 0.43867759, 0.40264556, 0.50034620, 0.54244138, 0.63776849, 0.47678383, 0.50539966, 0.56184085, 0.62312332, 0.58803372, 0.67276607, 0.71735313, 0.79126497, //10968-10975
	0.38805565, 0.41861076, 0.48377747, 0.54993139, 0.51118130, 0.60920199, 0.65700551, 0.74529285, 0.55622476, 0.57949642, 0.62221827, 0.67485964, 0.64515195, 0.71473516, 0.75474431, 0.81665460, 0.72650923, 0.74445856, 0.77469449, 0.81630604, 0.79382292, 0.84295948, 0.87645414, 0.91531624, 0.83343717, 0.84275725, 0.85753187, 0.88252372, 0.86936938, 0.89736581, 0.91983763, 0.94343650, //11000-11007
	0.00807107, 0.01390653, 0.01947894, 0.02933905, 0.02468785, 0.03353369, 0.04058065, 0.04501554, 0.01504980, 0.02715129, 0.03757117, 0.05953121, 0.04948574, 0.06811016, 0.08657605, 0.09758860, 0.03944992, 0.06421434, 0.08019895, 0.11890125, 0.10234803, 0.13172113, 0.16383920, 0.18129171, 0.06791052, 0.11294531, 0.13732058, 0.20821120, 0.18007949, 0.22815994, 0.29341387, 0.32396182, //11032-11039
	0.04962626, 0.08285330, 0.10286413, 0.15702934, 0.13461863, 0.17383309, 0.22350168, 0.24875429, 0.11067240, 0.17014216, 0.19637294, 0.28052664, 0.24941250, 0.30115129, 0.37638887, 0.40832357, 0.13874245, 0.21721110, 0.24842152, 0.36227921, 0.32296532, 0.38650083, 0.49738255, 0.53685300, 0.23176904, 0.33549724, 0.36627089, 0.49756623, 0.45826298, 0.51937777, 0.63773957, 0.67159691, //11064-11071
	0.09214721, 0.10774587, 0.16236328, 0.20371320, 0.17707048, 0.26564103, 0.29653372, 0.40008675, 0.17682064, 0.20168729, 0.27143720, 0.33013176, 0.29382986, 0.39887167, 0.44210694, 0.55074631, 0.25762039, 0.28262701, 0.34160306, 0.39877174, 0.36429047, 0.45768373, 0.50000700, 0.59572590, 0.44080480, 0.46853165, 0.52315051, 0.58460001, 0.54940455, 0.63437177, 0.68106912, 0.75851936, //11096-11103
	0.35440876, 0.38355277, 0.44564328, 0.51088519, 0.47270782, 0.56906109, 0.61861690, 0.70976582, 0.51830855, 0.54112919, 0.58307076, 0.63671317, 0.60645102, 0.67727282, 0.71963403, 0.78570541, 0.69694629, 0.71517199, 0.74584613, 0.78978163, 0.76606964, 0.81794747, 0.85501755, 0.89805972, 0.80910933, 0.81880948, 0.83429364, 0.86149610, 0.84713381, 0.87757235, 0.90326222, 0.93013853, //11128-11135
	0.00350407, 0.00662773, 0.00960632, 0.01570850, 0.01281411, 0.01827940, 0.02322625, 0.02632916, 0.00723095, 0.01439139, 0.02056703, 0.03558170, 0.02871365, 0.04144719, 0.05573638, 0.06428119, 0.02026676, 0.03609741, 0.04637844, 0.07473074, 0.06268290, 0.08405751, 0.11042948, 0.12473605, 0.03858266, 0.07047691, 0.08781453, 0.14525400, 0.12248047, 0.16136305, 0.22070381, 0.24855096, //11160-11167
	0.02683244, 0.04925350, 0.06287270, 0.10482209, 0.08742274, 0.11776161, 0.16117958, 0.18326428, 0.06594436, 0.11069862, 0.13034614, 0.20229600, 0.17578925, 0.21996520, 0.29126135, 0.32176828, 0.08710367, 0.14970824, 0.17497829, 0.27807158, 0.24227599, 0.30003298, 0.41294903, 0.45313368, 0.15716207, 0.24766899, 0.27472491, 0.40359284, 0.36489194, 0.42503038, 0.55401646, 0.59096199, //11192-11199
	0.05705944, 0.06800110, 0.10623337, 0.13847148, 0.11767953, 0.18658968, 0.21408047, 0.30418613, 0.12371434, 0.14320833, 0.19780369, 0.24865832, 0.21723516, 0.30827885, 0.34990678, 0.45483733, 0.18735693, 0.20791580, 0.25690486, 0.30857276, 0.27720863, 0.36188654, 0.40401815, 0.50027349, 0.35946624, 0.38505055, 0.43555169, 0.49723522, 0.46181895, 0.54727378, 0.59903359, 0.68470084, //11224-11231
	0.27831250, 0.30408264, 0.35899331, 0.42229208, 0.38534970, 0.47863397, 0.53171966, 0.62914515, 0.43230135, 0.45421877, 0.49445312, 0.55011036, 0.51874302, 0.59241394, 0.64034027, 0.71533073, 0.62979664, 0.64877359, 0.68048698, 0.72974090, 0.70319718, 0.76122898, 0.80643756, 0.85904532, 0.75398171, 0.76463873, 0.78164016, 0.81385388, 0.79688904, 0.83287512, 0.86561868, 0.89999960, //11256-11263
	0.02383334, 0.03502296, 0.04570322, 0.06057625, 0.05356197, 0.06685961, 0.07577319, 0.08136660, 0.03719684, 0.05696823, 0.07400026, 0.10273467, 0.08960435, 0.11395806, 0.13458058, 0.14684884, 0.08396003, 0.11851780, 0.14073189, 0.18572156, 0.16651969, 0.20058705, 0.23333516, 0.25115405, 0.12370088, 0.17843362, 0.20820765, 0.28148342, 0.25242081, 0.30206745, 0.36277914, 0.39125049, //11288-11295
	0.09815705, 0.14138595, 0.16746213, 0.22660860, 0.20211440, 0.24485530, 0.29295273, 0.31740554, 0.18332674, 0.24878328, 0.27769987, 0.35910295, 0.32899947, 0.37885771, 0.44569376, 0.47404321, 0.21415681, 0.29578041, 0.32827118, 0.43315035, 0.39698828, 0.45551229, 0.55122229, 0.58534892, 0.31704392, 0.41390553, 0.44300177, 0.55589502, 0.52202361, 0.57471456, 0.67369609, 0.70198030, //11320-11327
	0.15013544, 0.16973488, 0.23836677, 0.28242668, 0.25411433, 0.34846981, 0.37754335, 0.47449993, 0.24326271, 0.27058399, 0.34721970, 0.40360408, 0.36875635, 0.46979413, 0.50748015, 0.60231276, 0.33212585, 0.35749271, 0.41723867, 0.46949496, 0.43795271, 0.52329094, 0.55919920, 0.64060294, 0.50001306, 0.52564087, 0.57607789, 0.62893227, 0.59864881, 0.67176782, 0.71061788, 0.77500222, //11352-11359
	0.42088205, 0.44876089, 0.50828392, 0.56569834, 0.53206498, 0.61688603, 0.65838484, 0.73470962, 0.57102827, 0.59130137, 0.62854958, 0.67366023, 0.64822715, 0.70779293, 0.74287251, 0.79747307, 0.72424019, 0.74007896, 0.76667093, 0.80378105, 0.78375334, 0.82754674, 0.85937490, 0.89631433, 0.82046740, 0.82868213, 0.84176264, 0.86475445, 0.85262831, 0.87833761, 0.90091606, 0.92456760, //11384-11391
	0.01935643, 0.02897793, 0.03816981, 0.05147660, 0.04519830, 0.05709880, 0.06533571, 0.07051446, 0.03085299, 0.04829263, 0.06331777, 0.08969719, 0.07763563, 0.09998686, 0.11957262, 0.13123086, 0.07109377, 0.10244877, 0.12258185, 0.16492182, 0.14685287, 0.17890021, 0.21066435, 0.22792763, 0.10714682, 0.15814030, 0.18587499, 0.25668638, 0.22860471, 0.27656617, 0.33706981, 0.36542437, //11416-11423
	0.08398143, 0.12370514, 0.14768863, 0.20409194, 0.18071545, 0.22148693, 0.26887729, 0.29295918, 0.16118973, 0.22351254, 0.25099285, 0.33105753, 0.30145351, 0.35047152, 0.41803352, 0.44665961, 0.19057723, 0.26931146, 0.30065308, 0.40536293, 0.36923846, 0.42768136, 0.52594429, 0.56099193, 0.28842285, 0.38443012, 0.41332764, 0.52853029, 0.49394493, 0.54775191, 0.65119622, 0.68076123, //11448-11455
	0.13234873, 0.15049032, 0.21402320, 0.25623334, 0.22910417, 0.31962724, 0.34843887, 0.44452185, 0.22114379, 0.24718858, 0.32027309, 0.37584160, 0.34147820, 0.44105794, 0.47935570, 0.57563527, 0.30588626, 0.33050623, 0.38854927, 0.44081221, 0.40929059, 0.49456883, 0.53146322, 0.61503757, 0.47437862, 0.50001908, 0.55046484, 0.60464334, 0.57358069, 0.64853838, 0.68931744, 0.75684546, //11480-11487
	0.39494162, 0.42251849, 0.48135904, 0.53967697, 0.50551205, 0.59168479, 0.63483925, 0.71429576, 0.54566568, 0.56617200, 0.60393889, 0.65057880, 0.62427643, 0.68592764, 0.72292366, 0.78049103, 0.70591058, 0.72221920, 0.74962369, 0.78859899, 0.76757140, 0.81355706, 0.84759444, 0.88708981, 0.80653977, 0.81513730, 0.82881767, 0.85326612, 0.84038508, 0.86771132, 0.89204028, 0.91754229, //11512-11519
	0.01053928, 0.01709378, 0.02335719, 0.03355869, 0.02875559, 0.03788367, 0.04480349, 0.04914137, 0.01837504, 0.03122194, 0.04228135, 0.06402563, 0.05408681, 0.07251853, 0.09004051, 0.10048119, 0.04577540, 0.07083302, 0.08694128, 0.12400561, 0.10817831, 0.13620827, 0.16604549, 0.18225387, 0.07458248, 0.11821772, 0.14192498, 0.20786353, 0.18172369, 0.22633601, 0.28645140, 0.31457519, //11544-11551
	0.05607410, 0.08891282, 0.10871622, 0.15977484, 0.13864887, 0.17550264, 0.22148787, 0.24486916, 0.11773282, 0.17382459, 0.19846268, 0.27586107, 0.24727712, 0.29464816, 0.36354516, 0.39280893, 0.14414118, 0.21724106, 0.24622283, 0.35071061, 0.31469784, 0.37288004, 0.47618799, 0.51302216, 0.23207162, 0.32653869, 0.35481332, 0.47473120, 0.43872270, 0.49461512, 0.60688405, 0.63894765, //11576-11583
	0.09735326, 0.11263519, 0.16610870, 0.20481523, 0.17992458, 0.26287952, 0.29110471, 0.38548002, 0.17760020, 0.20114308, 0.26714902, 0.32114327, 0.28779246, 0.38443805, 0.42391588, 0.52329873, 0.25420404, 0.27742346, 0.33214880, 0.38436185, 0.35281826, 0.43802226, 0.47686608, 0.56461671, 0.42392081, 0.44957328, 0.50002919, 0.55681783, 0.52424809, 0.60282464, 0.64735959, 0.72111891, //11608-11615
	0.34390097, 0.37088379, 0.42842699, 0.48845515, 0.45325891, 0.54194047, 0.58848095, 0.67413878, 0.49574299, 0.51675151, 0.55537942, 0.60518796, 0.57712214, 0.64283517, 0.68367583, 0.74720041, 0.66980212, 0.68707551, 0.71604522, 0.75872250, 0.73572575, 0.78605711, 0.82439757, 0.86890799, 0.77914290, 0.78849246, 0.80335619, 0.83068162, 0.81629449, 0.84680636, 0.87459742, 0.90370941, //11640-11647
	0.00682822, 0.01147016, 0.01590342, 0.02363273, 0.01998809, 0.02689782, 0.03248908, 0.03600015, 0.01237047, 0.02189353, 0.03009891, 0.04742229, 0.03950836, 0.05419823, 0.06914347, 0.07804686, 0.03177307, 0.05113485, 0.06359078, 0.09426606, 0.08115967, 0.10440314, 0.13069387, 0.14502962, 0.05404816, 0.08954698, 0.10885822, 0.16667381, 0.14375562, 0.18291123, 0.23928736, 0.26571592, //11672-11679
	0.03972595, 0.06579697, 0.08152130, 0.12516533, 0.10707804, 0.13864691, 0.18071930, 0.20210688, 0.08739429, 0.13462260, 0.15545167, 0.22532944, 0.19953482, 0.24229932, 0.30884950, 0.33703630, 0.10967653, 0.17324845, 0.19855239, 0.29647876, 0.26269061, 0.31732372, 0.42139553, 0.45848861, 0.18378288, 0.27031358, 0.29625088, 0.41430241, 0.37889837, 0.43392535, 0.55201153, 0.58571311, //11704-11711
	0.07340076, 0.08573597, 0.12890168, 0.16210513, 0.14076315, 0.21188584, 0.23788855, 0.32442603, 0.14248689, 0.16270660, 0.21939595, 0.26862357, 0.23819200, 0.32637620, 0.36488001, 0.46162053, 0.20841399, 0.22902013, 0.27754465, 0.32650774, 0.29696870, 0.37679485, 0.41536647, 0.50278898, 0.37109048, 0.39538680, 0.44324007, 0.50002558, 0.46744602, 0.54605805, 0.59323374, 0.67142940, //11736-11743
	0.29439874, 0.31931563, 0.37247750, 0.43121138, 0.39682369, 0.48353437, 0.53188464, 0.62085056, 0.43989169, 0.46034095, 0.49793177, 0.54867075, 0.52006092, 0.58713971, 0.63108463, 0.69953001, 0.62503742, 0.64263644, 0.67220191, 0.71786634, 0.69323101, 0.74711031, 0.79040150, 0.84064835, 0.74130189, 0.75119976, 0.76695659, 0.79715403, 0.78123711, 0.81498200, 0.84720360, 0.88096470, //11768-11775
	0.00895468, 0.01469525, 0.02017800, 0.02932432, 0.02500899, 0.03319402, 0.03954666, 0.04353056, 0.01581277, 0.02724171, 0.03709115, 0.05693406, 0.04786286, 0.06469491, 0.08112248, 0.09090527, 0.03982387, 0.06242283, 0.07696472, 0.11130777, 0.09663479, 0.12262587, 0.15095997, 0.16633367, 0.06579655, 0.10595471, 0.12774977, 0.19030385, 0.16549727, 0.20781103, 0.26632422, 0.29373803, //11800-11807
	0.04909609, 0.07903898, 0.09710854, 0.14500663, 0.12517132, 0.15979782, 0.20409305, 0.22659410, 0.10482771, 0.15703886, 0.18012395, 0.25427073, 0.22687281, 0.27230663, 0.34020849, 0.36899016, 0.12944002, 0.19844707, 0.22594005, 0.32754124, 0.29247140, 0.34917597, 0.45280439, 0.48974415, 0.21145022, 0.30250176, 0.32989856, 0.44892858, 0.41320436, 0.46874612, 0.58347068, 0.61624587, //11832-11839
	0.08712455, 0.10115620, 0.15024431, 0.18656966, 0.16319548, 0.24107565, 0.26837501, 0.35927567, 0.16261366, 0.18475116, 0.24680018, 0.29873633, 0.26661948, 0.35968610, 0.39873841, 0.49691492, 0.23468643, 0.25676647, 0.30893269, 0.35971223, 0.32903335, 0.41207706, 0.45063504, 0.53824882, 0.40138359, 0.42643539, 0.47577402, 0.53256383, 0.50001424, 0.57862110, 0.62425133, 0.69991312, //11864-11871
	0.32277789, 0.34885694, 0.40454347, 0.46402517, 0.42917424, 0.51704556, 0.56435092, 0.65143259, 0.47191397, 0.49267014, 0.53085714, 0.58105469, 0.55274145, 0.61902430, 0.66122184, 0.72685248, 0.65070288, 0.66811023, 0.69732854, 0.74128275, 0.71757855, 0.76942038, 0.80988789, 0.85685317, 0.76298974, 0.77257398, 0.78783720, 0.81636556, 0.80131979, 0.83321742, 0.86290843, 0.89398631, //11896-11903
	0.00381450, 0.00690738, 0.00986073, 0.01558814, 0.01288737, 0.01801185, 0.02250956, 0.02533969, 0.00750909, 0.01433618, 0.02022964, 0.03397236, 0.02768501, 0.03935197, 0.05220311, 0.05985866, 0.02045736, 0.03516759, 0.04468379, 0.07016006, 0.05927802, 0.07858268, 0.10209051, 0.11487527, 0.03740304, 0.06632005, 0.08208253, 0.13328682, 0.11298144, 0.14769568, 0.20106853, 0.22613986, //11928-11935
	0.02649789, 0.04705555, 0.05951452, 0.09711145, 0.08153109, 0.10873275, 0.14768867, 0.16747982, 0.06287175, 0.10293101, 0.12059954, 0.18447926, 0.16083353, 0.19988219, 0.26449346, 0.29188069, 0.08164374, 0.13759187, 0.15981099, 0.25253318, 0.22053647, 0.27230859, 0.37698148, 0.41432156, 0.14456639, 0.22477025, 0.24905242, 0.36542680, 0.33045838, 0.38482030, 0.50753817, 0.54258885, //11960-11967
	0.05398609, 0.06397335, 0.09871688, 0.12753611, 0.10906586, 0.17054971, 0.19463934, 0.27480748, 0.11402993, 0.13154598, 0.18067378, 0.22600646, 0.19799069, 0.27936135, 0.31705207, 0.41179336, 0.17133617, 0.18978372, 0.23327947, 0.27955318, 0.25168903, 0.32724557, 0.36562647, 0.45270767, 0.32827897, 0.35146581, 0.39715713, 0.45397451, 0.42140148, 0.50002526, 0.54937344, 0.63120386, //11992-11999
	0.25428916, 0.27753176, 0.32722188, 0.38485705, 0.35110066, 0.43621296, 0.48597001, 0.57773396, 0.39465615, 0.41463792, 0.45136256, 0.50286947, 0.47382453, 0.54194403, 0.58851198, 0.66094077, 0.58875350, 0.60661662, 0.63667165, 0.68474873, 0.65882425, 0.71552798, 0.76285388, 0.81773006, 0.71062585, 0.72097664, 0.73742258, 0.76998320, 0.75283283, 0.78919969, 0.82501826, 0.86248121, //12024-12031
	0.00279449, 0.00515262, 0.00740015, 0.01193077, 0.00979329, 0.01384307, 0.01754369, 0.01986762, 0.00561042, 0.01095053, 0.01555640, 0.02677169, 0.02164283, 0.03114611, 0.04210435, 0.04863030, 0.01546713, 0.02721762, 0.03477888, 0.05598564, 0.04692331, 0.06299784, 0.08334431, 0.09441744, 0.02898732, 0.05281318, 0.06578825, 0.11001140, 0.09246945, 0.12242420, 0.17096147, 0.19372489, //12056-12063
	0.02030097, 0.03699072, 0.04705935, 0.07904766, 0.06579713, 0.08888783, 0.12371325, 0.14138005, 0.04924481, 0.08282015, 0.09765554, 0.15381933, 0.13301883, 0.16741439, 0.22710714, 0.25240668, 0.06508658, 0.11303701, 0.13212249, 0.21607765, 0.18711889, 0.23394515, 0.33462722, 0.37052574, 0.11778244, 0.18891737, 0.21026955, 0.31927921, 0.28655935, 0.33746169, 0.45944795, 0.49430113, //12088-12095
	0.04303961, 0.05121418, 0.07981493, 0.10425369, 0.08853571, 0.14098697, 0.16227648, 0.23333422, 0.09461236, 0.10960500, 0.15168606, 0.19198187, 0.16707974, 0.23927316, 0.27430537, 0.36243543, 0.14386219, 0.15993722, 0.19786910, 0.23970875, 0.21444325, 0.28266532, 0.31888731, 0.40109750, 0.28940920, 0.31072266, 0.35266800, 0.40678600, 0.37577479, 0.45063183, 0.50002505, 0.58184560, //12120-12127
	0.22079028, 0.24178013, 0.28661692, 0.34051200, 0.30894439, 0.38855381, 0.43749427, 0.52758226, 0.35097035, 0.36957486, 0.40378821, 0.45367524, 0.42554200, 0.49145779, 0.53860300, 0.61193686, 0.54808331, 0.56557020, 0.59494401, 0.64380891, 0.61745254, 0.67509153, 0.72566676, 0.78435183, 0.67187467, 0.68227193, 0.69886267, 0.73288489, 0.71497009, 0.75301063, 0.79236558, 0.83360772, //12152-12159
	0.00109823, 0.00224081, 0.00333097, 0.00586569, 0.00466736, 0.00694772, 0.00930741, 0.01079101, 0.00246366, 0.00534040, 0.00781435, 0.01482514, 0.01162440, 0.01755584, 0.02535875, 0.02999980, 0.00723546, 0.01404224, 0.01841876, 0.03246279, 0.02646757, 0.03708577, 0.05231482, 0.06059176, 0.01508130, 0.03041842, 0.03874756, 0.07144300, 0.05847603, 0.08059861, 0.12106008, 0.14001594, //12184-12191
	0.01003160, 0.02029823, 0.02647603, 0.04910626, 0.03972334, 0.05602322, 0.08398327, 0.09815270, 0.02681551, 0.04963044, 0.05968298, 0.10302554, 0.08702141, 0.11349815, 0.16508109, 0.18691645, 0.03757308, 0.07235596, 0.08619440, 0.15563811, 0.13166756, 0.17048594, 0.26445712, 0.29799169, 0.07356628, 0.12948256, 0.14617683, 0.24285685, 0.21385934, 0.25899452, 0.37976746, 0.41427192, //12216-12223
	0.02488265, 0.03009609, 0.04840658, 0.06570739, 0.05457723, 0.09159137, 0.10843568, 0.16441951, 0.06243906, 0.07328284, 0.10363246, 0.13551548, 0.11580665, 0.17297360, 0.20355722, 0.28045839, 0.09830867, 0.11050754, 0.13918676, 0.17348541, 0.15274687, 0.20887585, 0.24143132, 0.31532984, 0.22500375, 0.24316376, 0.27891591, 0.32859673, 0.30011247, 0.36881758, 0.41820971, 0.50001981, //12248-12255
	0.16528147, 0.18251652, 0.21922835, 0.26707729, 0.23903016, 0.30968422, 0.35717691, 0.44459282, 0.27860319, 0.29492493, 0.32500053, 0.37205890, 0.34553893, 0.40775578, 0.45586871, 0.53082328, 0.48071157, 0.49753311, 0.52577140, 0.57593844, 0.54888067, 0.60805182, 0.66404288, 0.72906232, 0.60762528, 0.61817058, 0.63492123, 0.67144086, 0.65217459, 0.69298606, 0.73828083, 0.78572575, //12280-12287
	0.03714499, 0.05291766, 0.06797486, 0.08749045, 0.07829022, 0.09574898, 0.10667961, 0.11354027, 0.05597707, 0.08262639, 0.10559673, 0.14138381, 0.12502391, 0.15536225, 0.17921301, 0.19341650, 0.12190620, 0.16596274, 0.19425783, 0.24746847, 0.22474325, 0.26502226, 0.30104344, 0.32062231, 0.17256790, 0.23866267, 0.27458719, 0.35588217, 0.32363333, 0.37868020, 0.44082610, 0.46994090, //12288-12319
	0.14001096, 0.19372471, 0.22612218, 0.29372998, 0.26573140, 0.31458719, 0.36542868, 0.39125696, 0.24851838, 0.32398317, 0.35726979, 0.44373801, 0.41179553, 0.46473461, 0.53034408, 0.55818382, 0.28407251, 0.37501565, 0.41120488, 0.51799667, 0.48117475, 0.54073092, 0.62991250, 0.66171571, 0.40259786, 0.50329960, 0.53347437, 0.64107413, 0.60879185, 0.65899992, 0.74546557, 0.77019664, //12320-12351
	0.20286252, 0.22685371, 0.31081170, 0.36059342, 0.32860491, 0.43525539, 0.46550936, 0.56615534, 0.30948428, 0.34079174, 0.42855975, 0.48805839, 0.45128602, 0.55789389, 0.59430167, 0.68583463, 0.41127985, 0.43897634, 0.50429558, 0.55733911, 0.52532094, 0.61192293, 0.64555615, 0.72177828, 0.57914161, 0.60507283, 0.65608639, 0.70561566, 0.67723246, 0.74573755, 0.77923694, 0.83475501, //12352-12383
	0.50000979, 0.52913123, 0.59128020, 0.64681423, 0.61428728, 0.69629985, 0.73312853, 0.80082369, 0.65017753, 0.66989983, 0.70621041, 0.74722105, 0.72412583, 0.77823277, 0.80770028, 0.85355492, 0.78296258, 0.79751291, 0.82192780, 0.85361944, 0.83653071, 0.87391246, 0.89896008, 0.92803842, 0.86636426, 0.87349771, 0.88483892, 0.90349172, 0.89365029, 0.91451482, 0.93150059, 0.94930939, //12384-12415
	0.03057879, 0.04432095, 0.05743420, 0.07509286, 0.06676992, 0.08254708, 0.09276212, 0.09917338, 0.04699019, 0.07075574, 0.09123487, 0.12441476, 0.10924497, 0.13736979, 0.16023157, 0.17384361, 0.10443913, 0.14488249, 0.17090182, 0.22148882, 0.19989478, 0.23821178, 0.27355580, 0.29275491, 0.15095619, 0.21321736, 0.24708055, 0.32653459, 0.29500164, 0.34880628, 0.41146463, 0.44082756, //12416-12447
	0.12105489, 0.17096846, 0.20109718, 0.26632574, 0.23929386, 0.28645732, 0.33706512, 0.36278615, 0.22069851, 0.29343168, 0.32547999, 0.41163033, 0.37977483, 0.43252058, 0.49979649, 0.52832972, 0.25496484, 0.34378068, 0.37913979, 0.48714307, 0.44987736, 0.51015136, 0.60302845, 0.63614252, 0.36911770, 0.47037030, 0.50073309, 0.61229151, 0.57883400, 0.63087356, 0.72283737, 0.74910845, //12448-12479
	0.18017310, 0.20260666, 0.28116689, 0.32948340, 0.29841730, 0.40207687, 0.43243211, 0.53367891, 0.28280526, 0.31293505, 0.39745237, 0.45682239, 0.42011839, 0.52651440, 0.56408143, 0.65848680, 0.38077873, 0.40800586, 0.47226216, 0.52601386, 0.49357669, 0.58142414, 0.61647326, 0.69597752, 0.55126445, 0.57749919, 0.62915860, 0.68069494, 0.65115673, 0.72247546, 0.75823863, 0.81751804, //12480-12511
	0.47087872, 0.50000449, 0.56218271, 0.61937129, 0.58586545, 0.67034387, 0.70931172, 0.78105006, 0.62340165, 0.64367471, 0.68096344, 0.72410484, 0.69978732, 0.75673728, 0.78843717, 0.83777295, 0.76447495, 0.77967747, 0.80520542, 0.83906390, 0.82080489, 0.86076046, 0.88808357, 0.91980572, 0.85307840, 0.86066509, 0.87273730, 0.89294539, 0.88229732, 0.90488804, 0.92361374, 0.94324356, //12512-12543
	0.01657423, 0.02598370, 0.03496842, 0.04862228, 0.04218597, 0.05439348, 0.06306321, 0.06849348, //12544-12551
	0.02780169, 0.04540107, 0.06057031, 0.08819524, 0.07556806, 0.09898611, 0.11971451, 0.13205551, //12552-12559
	0.06717394, 0.09989087, 0.12099011, 0.16604264, 0.14682746, 0.18100652, 0.21489610, 0.23327260, //12560-12567
	0.10485590, 0.15887186, 0.18827971, 0.26386712, 0.23388967, 0.28506843, 0.34879575, 0.37868830, //12568-12575
	0.08057954, 0.12241680, 0.14770820, 0.20782126, 0.18289085, 0.22639500, 0.27653591, 0.30203119, //12576-12583
	0.16133373, 0.22811064, 0.25762563, 0.34307755, 0.31151848, 0.36389913, 0.43459050, 0.46461272, //12584-12591
	0.19282122, 0.27709986, 0.31069607, 0.42132763, 0.38307388, 0.44488519, 0.54563904, 0.58157331, //12592-12599
	0.29797450, 0.40020737, 0.43097678, 0.55090106, 0.51495846, 0.57086299, 0.67451841, 0.70413552, //12600-12607
	0.13176414, 0.15089048, 0.21790209, 0.26309339, 0.23404369, 0.33092722, 0.36173378, 0.46441461, //12608-12615
	0.22581567, 0.25350542, 0.33104693, 0.39014112, 0.35361507, 0.45950581, 0.49950297, 0.60008828, //12616-12623
	0.31558864, 0.34189802, 0.40384656, 0.45921632, 0.42582979, 0.51604255, 0.55434681, 0.64085161, //12624-12631
	0.49173953, 0.51866738, 0.57163429, 0.62748531, 0.59551314, 0.67271408, 0.71342509, 0.78077565, //12632-12639
	0.40872310, 0.43784287, 0.50000706, 0.56075641, 0.52519511, 0.61493019, 0.65854194, 0.73879369, //12640-12647
	0.56625554, 0.58766556, 0.62700239, 0.67473406, 0.64787183, 0.71082337, 0.74730341, 0.80408075, //12648-12655
	0.72500322, 0.74161755, 0.76951518, 0.80800697, 0.78722562, 0.83269440, 0.86487345, 0.90222229, //12656-12663
	0.82471763, 0.83329305, 0.84689564, 0.87046783, 0.85804668, 0.88439154, 0.90679564, 0.93026675, //12664-12671
	0.01105198, 0.01789083, 0.02442252, 0.03501941, 0.03002192, 0.03949564, 0.04666012, 0.05114017, //12672-12679
	0.01922387, 0.03258131, 0.04409324, 0.06658896, 0.05629837, 0.07536946, 0.09339542, 0.10412643, //12680-12687
	0.04781132, 0.07377263, 0.09046931, 0.12869931, 0.11237867, 0.14134399, 0.17190072, 0.18851168, //12688-12695
	0.07767105, 0.12273351, 0.14718719, 0.21489642, 0.18805328, 0.23389375, 0.29499496, 0.32363350, //12696-12703
	0.05847676, 0.09247201, 0.11297973, 0.16550851, 0.14372813, 0.18171624, 0.22858818, 0.25242355, //12704-12711
	0.12246765, 0.18015228, 0.20565055, 0.28478195, 0.25553431, 0.30395616, 0.37383763, 0.40349593, //12712-12719
	0.14964573, 0.22470703, 0.25457695, 0.36077954, 0.32412287, 0.38341935, 0.48727396, 0.52429952, //12720-12727
	0.24032905, 0.33661407, 0.36561139, 0.48678059, 0.45036309, 0.50698627, 0.61889743, 0.65087828, //12728-12735
	0.10110636, 0.11690243, 0.17214925, 0.21194409, 0.18638117, 0.27157537, 0.30045784, 0.39683395, //12736-12743
	0.18355938, 0.20776678, 0.27565423, 0.33072741, 0.29667052, 0.39540586, 0.43542463, 0.53596789, //12744-12751
	0.26223614, 0.28605860, 0.34215936, 0.39529759, 0.36320312, 0.45001425, 0.48917266, 0.57761361, //12752-12759
	0.43431345, 0.46034375, 0.51159414, 0.56879779, 0.53602686, 0.61522934, 0.65950606, 0.73296180, //12760-12767
	0.35319677, 0.38065415, 0.43926994, 0.50000792, 0.46444151, 0.55415250, 0.60070581, 0.68634690, //12768-12775
	0.50710264, 0.52837268, 0.56752492, 0.61748214, 0.58934193, 0.65531412, 0.69580424, 0.75882004, //12776-12783
	0.67988305, 0.69718206, 0.72623554, 0.76854978, 0.74571519, 0.79564798, 0.83312814, 0.87661049, //12784-12791
	0.78840185, 0.79768978, 0.81251373, 0.83933841, 0.82519115, 0.85519057, 0.88207916, 0.91024623, //12792-12799
	0.01428335, 0.02263347, 0.03059208, 0.04298605, 0.03714370, 0.04824137, 0.05627055, 0.06131239, //12800-12807
	0.02424942, 0.04009362, 0.05375253, 0.07924996, 0.06758804, 0.08919186, 0.10880738, 0.12048364, //12808-12815
	0.05914635, 0.08906623, 0.10831666, 0.15063183, 0.13257810, 0.16456189, 0.19705919, 0.21471530, //12816-12823
	0.09356352, 0.14389812, 0.17126833, 0.24358747, 0.21491115, 0.26386632, 0.32650912, 0.35588946, //12824-12831
	0.07142846, 0.11000857, 0.13329211, 0.19030147, 0.16668383, 0.20788393, 0.25668514, 0.28148555, //12832-12839
	0.14525960, 0.20825383, 0.23606746, 0.31893368, 0.28829573, 0.33904485, 0.40940507, 0.43929860, //12840-12847
	0.17487981, 0.25538898, 0.28743838, 0.39621970, 0.35867477, 0.41941018, 0.52145564, 0.55785718, //12848-12855
	0.27389547, 0.37378830, 0.40377244, 0.52435512, 0.48816845, 0.54440342, 0.65148057, 0.68207633, //12856-12863
	0.11905207, 0.13679899, 0.19893303, 0.24193607, 0.21426772, 0.30636493, 0.33640373, 0.43626485, //12864-12871
	0.20831368, 0.23454076, 0.30808817, 0.36554612, 0.33000519, 0.43292579, 0.47295335, 0.57349125, //12872-12879
	0.29350874, 0.31875122, 0.37820884, 0.43267749, 0.39990005, 0.48871033, 0.52734332, 0.61470750, //12880-12887
	0.46795690, 0.49451196, 0.54674052, 0.60320485, 0.57086079, 0.64889654, 0.69108010, 0.76095054, //12888-12895
	0.38570741, 0.41415789, 0.47484586, 0.53560294, 0.50000172, 0.58974403, 0.63457869, 0.71705307, //12896-12903
	0.54176622, 0.56311095, 0.60243113, 0.65104427, 0.62360841, 0.68781198, 0.72599177, 0.78527754, //12904-12911
	0.70631686, 0.72320246, 0.75156496, 0.79165577, 0.77002986, 0.81734680, 0.85171564, 0.89160790, //12912-12919
	0.80967603, 0.81854603, 0.83266665, 0.85756463, 0.84443650, 0.87228723, 0.89654811, 0.92198074, //12920-12927
	0.00613216, 0.01068867, 0.01501701, 0.02289353, 0.01917509, 0.02622384, 0.03202041, 0.03566318, //12928-12935
	0.01155999, 0.02114515, 0.02939168, 0.04733409, 0.03913596, 0.05432211, 0.06992964, 0.07924210, //12936-12943
	0.03058823, 0.05048875, 0.06324829, 0.09538050, 0.08168702, 0.10604221, 0.13359327, 0.14860185, //12944-12951
	0.05352494, 0.09047527, 0.11055493, 0.17125523, 0.14720690, 0.18829906, 0.24704517, 0.27459248, //12952-12959
	0.03881168, 0.06575481, 0.08201707, 0.12777605, 0.10882501, 0.14187698, 0.18585943, 0.20820930, //12960-12967
	0.08770604, 0.13738659, 0.15933769, 0.23275713, 0.20556643, 0.25072091, 0.31969882, 0.34902526, //12968-12975
	0.11119860, 0.17797431, 0.20460399, 0.30683840, 0.27150389, 0.32852777, 0.43524447, 0.47326514, //12976-12983
	0.18889701, 0.27992374, 0.30736777, 0.42967926, 0.39290649, 0.44996121, 0.56930623, 0.60343045, //12984-12991
	0.07381347, 0.08658721, 0.13139992, 0.16628064, 0.14390273, 0.21857215, 0.24590369, 0.33685412, //12992-12999
	0.14588535, 0.16700365, 0.22623120, 0.27776926, 0.24594306, 0.33829845, 0.37826412, 0.47873102, //13000-13007
	0.21466741, 0.23624859, 0.28718262, 0.33850958, 0.30754022, 0.39115883, 0.43094177, 0.52134929, //13008-13015
	0.38312256, 0.40835406, 0.45803565, 0.51647646, 0.48302518, 0.56376547, 0.61145246, 0.69030084, //13016-13023
	0.30372571, 0.32969015, 0.38511212, 0.44585066, 0.41025153, 0.49999296, 0.54920403, 0.63961519, //13024-13031
	0.45437361, 0.47553458, 0.51436669, 0.56639089, 0.53707100, 0.60582399, 0.64995291, 0.71857564, //13032-13039
	0.63966474, 0.65758478, 0.68766735, 0.73336983, 0.70871814, 0.76265893, 0.80482581, 0.85375796, //13040-13047
	0.75601100, 0.76595959, 0.78182068, 0.81157518, 0.79591879, 0.82919332, 0.86004617, 0.89241554, //13048-13055
	0.00461388, 0.00817250, 0.01156883, 0.01794005, 0.01493631, 0.02063622, 0.02549401, 0.02854620, //13056-13063
	0.00886678, 0.01653835, 0.02315571, 0.03807465, 0.03125273, 0.04389185, 0.05745240, 0.06552307, //13064-13071
	0.02375086, 0.03999946, 0.05048358, 0.07770983, 0.06605722, 0.08671075, 0.11113502, 0.12436752, //13072-13079
	0.04243710, 0.07357277, 0.09045687, 0.14392128, 0.12271014, 0.15889093, 0.21321838, 0.23865624, //13080-13087
	0.03041547, 0.05279921, 0.06634799, 0.10594270, 0.08954893, 0.11819531, 0.15813679, 0.17844161, //13088-13095
	0.07044133, 0.11294894, 0.13164487, 0.19757986, 0.17320195, 0.21349523, 0.27859574, 0.30618152, //13096-13103
	0.09044149, 0.14894251, 0.17227315, 0.26651557, 0.23399948, 0.28661857, 0.39078629, 0.42795818, //13104-13111
	0.15710508, 0.23947519, 0.26424589, 0.38076557, 0.34580775, 0.40030438, 0.52099360, 0.55548044, //13112-13119
	0.05997723, 0.07071759, 0.10828868, 0.13855319, 0.11908353, 0.18398296, 0.20858227, 0.29091234, //13120-13127
	0.12296646, 0.14138634, 0.19301794, 0.23964158, 0.21082100, 0.29439101, 0.33222534, 0.42727497, //13128-13135
	0.18309724, 0.20226661, 0.24748236, 0.29453681, 0.26614306, 0.34309171, 0.38137367, 0.46836006, //13136-13143
	0.34165078, 0.36519091, 0.41150297, 0.46811978, 0.43568312, 0.51405591, 0.56254174, 0.64287782, //13144-13151
	0.26690229, 0.29070241, 0.34150400, 0.39931767, 0.36544106, 0.45085625, 0.50002486, 0.59048293, //13152-13159
	0.40868676, 0.42878787, 0.46560955, 0.51683445, 0.48798071, 0.55551096, 0.60110160, 0.67190530, //13160-13167
	0.59943210, 0.61719485, 0.64701886, 0.69418785, 0.66873155, 0.72439150, 0.77029389, 0.82359875, //13168-13175
	0.71922082, 0.72939198, 0.74557393, 0.77729080, 0.76058659, 0.79601232, 0.83065043, 0.86691939, //13176-13183
	0.00181564, 0.00355301, 0.00520571, 0.00881092, 0.00711105, 0.01034566, 0.01347607, 0.01545646, //13184-13191
	0.00389122, 0.00805904, 0.01164388, 0.02105150, 0.01675410, 0.02473441, 0.03448114, 0.04028208, //13192-13199
	0.01114463, 0.02066368, 0.02681253, 0.04515892, 0.03732828, 0.05121610, 0.06974507, 0.07982904, //13200-13207
	0.02211605, 0.04243685, 0.05347308, 0.09356073, 0.07767246, 0.10485265, 0.15094703, 0.17257708, //13208-13215
	0.01506119, 0.02898137, 0.03739812, 0.06581118, 0.05403590, 0.07458998, 0.10714171, 0.12368912, //13216-13223
	0.03857884, 0.06799740, 0.08102928, 0.13283466, 0.11363459, 0.14539133, 0.20301456, 0.22737477, //13224-13231
	0.05240255, 0.09555837, 0.11279853, 0.19238444, 0.16489278, 0.20932450, 0.30898831, 0.34452956, //13232-13239
	0.09866688, 0.16482676, 0.18476143, 0.29113374, 0.25920150, 0.30887982, 0.43202051, 0.46731250, //13240-13247
	0.03453661, 0.04149238, 0.06581884, 0.08747501, 0.07351212, 0.12021101, 0.14009086, 0.20650808, //13248-13255
	0.08082777, 0.09425325, 0.13188176, 0.16944379, 0.14625500, 0.21350435, 0.24750511, 0.33280231, //13256-13263
	0.12498550, 0.13973649, 0.17434799, 0.21397731, 0.19010110, 0.25469522, 0.29015129, 0.37072915, //13264-13271
	0.26533948, 0.28572733, 0.32589604, 0.37921173, 0.34861691, 0.42240511, 0.47239281, 0.55534726, //13272-13279
	0.19917454, 0.21893978, 0.26125018, 0.31365651, 0.28296103, 0.36039202, 0.40951476, 0.49999878, //13280-13287
	0.32469125, 0.34273263, 0.37590808, 0.42555677, 0.39761563, 0.46320314, 0.51121610, 0.58624062, //13288-13295
	0.52547296, 0.54291059, 0.57221469, 0.62206512, 0.59516420, 0.65398036, 0.70676251, 0.76803918, //13296-13303
	0.65154352, 0.66210961, 0.67890832, 0.71417839, 0.69559099, 0.73495425, 0.77649494, 0.82001533, //13304-13311
	0.01188617, 0.01896438, 0.02570803, 0.03641717, 0.03137516, 0.04094761, 0.04803631, 0.05249528, //13312-13319
	0.02034104, 0.03393811, 0.04564269, 0.06805196, 0.05780965, 0.07678664, 0.09451189, 0.10506090, //13320-13327
	0.04987006, 0.07592713, 0.09265279, 0.13031786, 0.11422623, 0.14271551, 0.17256920, 0.18878659, //13328-13335
	0.07985013, 0.12439848, 0.14864841, 0.21472157, 0.18851320, 0.23324041, 0.29275401, 0.32059483, //13336-13343
	0.06058804, 0.09443229, 0.11484883, 0.16635880, 0.14501202, 0.18220781, 0.22792922, 0.25112554, //13344-13351
	0.12474560, 0.18126654, 0.20630100, 0.28311309, 0.25474803, 0.30174002, 0.36971037, 0.39853518, //13352-13359
	0.15140496, 0.22466003, 0.25385009, 0.35702667, 0.32142046, 0.37898930, 0.48060497, 0.51683791, //13360-13367
	0.24008126, 0.33364980, 0.36167431, 0.47949861, 0.44407974, 0.49911402, 0.60928051, 0.64076533, //13368-13375
	0.10280179, 0.11847858, 0.17333781, 0.21222885, 0.18720503, 0.27068185, 0.29863041, 0.39191489, //13376-13383
	0.18378135, 0.20756761, 0.27424993, 0.32780504, 0.29468330, 0.39074881, 0.42960919, 0.52726652, //13384-13391
	0.26111616, 0.28437102, 0.33911506, 0.39068164, 0.35947198, 0.44368046, 0.48175552, 0.56787943, //13392-13399
	0.42900448, 0.45436182, 0.50430810, 0.56013065, 0.52810987, 0.60535720, 0.64903772, 0.72143134, //13400-13407
	0.34985055, 0.37660955, 0.43374432, 0.49295740, 0.45828804, 0.54558350, 0.59132024, 0.67536035, //13408-13415
	0.50002056, 0.52072872, 0.55886087, 0.60766562, 0.58016674, 0.64458882, 0.68470095, 0.74709451, //13416-13423
	0.67152097, 0.68854063, 0.71708424, 0.75904425, 0.73639374, 0.78593836, 0.82385098, 0.86784912, //13424-13431
	0.77925563, 0.78844664, 0.80309579, 0.82998369, 0.81580563, 0.84588359, 0.87347253, 0.90234961, //13432-13439
	0.00976951, 0.01588523, 0.02172515, 0.03130713, 0.02679032, 0.03535900, 0.04191719, 0.04603763, //13440-13447
	0.01707218, 0.02911200, 0.03948223, 0.06005147, 0.05064546, 0.06808223, 0.08485831, 0.09484819, //13448-13455
	0.04263691, 0.06620762, 0.08137729, 0.11660852, 0.10157635, 0.12826811, 0.15699149, 0.17256846, //13456-13463
	0.06975213, 0.11112730, 0.13360977, 0.19707969, 0.17189784, 0.21485120, 0.27352895, 0.30104178, //13464-13471
	0.05232249, 0.08334769, 0.10207964, 0.15096525, 0.13071326, 0.16603716, 0.21066031, 0.23333674, //13472-13479
	0.11038265, 0.16381580, 0.18742014, 0.26222471, 0.23457514, 0.28046944, 0.34818753, 0.37694381, //13480-13487
	0.13562178, 0.20573369, 0.23365902, 0.33555420, 0.30038580, 0.35724167, 0.46011806, 0.49681272, //13488-13495
	0.21940102, 0.31105775, 0.33847462, 0.45701339, 0.42147877, 0.47674538, 0.59001735, 0.62238061, //13496-13503
	0.09151462, 0.10600410, 0.15667409, 0.19368399, 0.16989608, 0.24928524, 0.27668765, 0.36830739, //13504-13511
	0.16847585, 0.19105812, 0.25432945, 0.30667369, 0.27432857, 0.36809801, 0.40703072, 0.50492182, //13512-13519
	0.24196723, 0.26433160, 0.31708278, 0.36809657, 0.33733116, 0.42051775, 0.45894313, 0.54598373, //13520-13527
	0.40873421, 0.43384911, 0.48326437, 0.53969105, 0.50734910, 0.58542796, 0.63043068, 0.70503848, //13528-13535
	0.33009780, 0.35634751, 0.41234904, 0.47165538, 0.43691779, 0.52447871, 0.57124353, 0.65726620, //13536-13543
	0.47928503, 0.50000036, 0.53812810, 0.58785367, 0.55981492, 0.62541939, 0.66698565, 0.73159043, //13544-13551
	0.65589562, 0.67316045, 0.70213692, 0.74550002, 0.72211755, 0.77327400, 0.81301819, 0.85914833, //13552-13559
	0.76680489, 0.77626843, 0.79134022, 0.81937241, 0.80461283, 0.83595074, 0.86509310, 0.89560575, //13560-13567
	0.00587243, 0.01022350, 0.01437313, 0.02191303, 0.01835830, 0.02509677, 0.03067261, 0.03414770, //13568-13575
	0.01107312, 0.02024086, 0.02815504, 0.04532551, 0.03745536, 0.05204461, 0.06710272, 0.07604934, //13576-13583
	0.02929044, 0.04834152, 0.06057578, 0.09139756, 0.07828499, 0.10156902, 0.12827174, 0.14277124, //13584-13591
	0.05113767, 0.08672571, 0.10606862, 0.16451961, 0.14134202, 0.18098963, 0.23820789, 0.26502582, //13592-13599
	0.03708230, 0.06296553, 0.07858022, 0.12264239, 0.10439530, 0.13626710, 0.17887109, 0.20057537, //13600-13607
	0.08410210, 0.13175359, 0.15259670, 0.22389346, 0.19750213, 0.24109777, 0.30863227, 0.33724938, //13608-13615
	0.10647973, 0.17092048, 0.19651962, 0.29602977, 0.26172593, 0.31722924, 0.42244970, 0.45995212, //13616-13623
	0.18142000, 0.26950869, 0.29589126, 0.41565483, 0.37978981, 0.43570852, 0.55460171, 0.58855518, //13624-13631
	0.07073741, 0.08303001, 0.12614261, 0.15961674, 0.13802073, 0.20997677, 0.23639215, 0.32445283, //13632-13639
	0.14029218, 0.16067346, 0.21785943, 0.26771664, 0.23685714, 0.32641527, 0.36549626, 0.46385206, //13640-13647
	0.20668208, 0.22754534, 0.27677022, 0.32645786, 0.29653240, 0.37763175, 0.41695114, 0.50557353, //13648-13655
	0.37148660, 0.39608917, 0.44463380, 0.50217023, 0.46916804, 0.54883613, 0.59621001, 0.67502115, //13656-13663
	0.29377651, 0.31903799, 0.37299632, 0.43250102, 0.39766610, 0.48564332, 0.53443322, 0.62413336, //13664-13671
	0.44116975, 0.46191596, 0.49996756, 0.55134284, 0.52236269, 0.59027345, 0.63445766, 0.70318214, //13672-13679
	0.62711327, 0.64487729, 0.67472503, 0.72057632, 0.69584855, 0.75001870, 0.79310927, 0.84315285, //13680-13687
	0.74387608, 0.75383377, 0.76968556, 0.79986634, 0.78398911, 0.81772720, 0.84967745, 0.88317173, //13688-13695
	0.00383857, 0.00693094, 0.00986561, 0.01557845, 0.01288884, 0.01798502, 0.02247175, 0.02528366, //13696-13703
	0.00752614, 0.01433512, 0.02019777, 0.03387660, 0.02761576, 0.03922177, 0.05199409, 0.05959998, //13704-13711
	0.02044775, 0.03512332, 0.04452355, 0.06990970, 0.05906101, 0.07825864, 0.10158975, 0.11427618, //13712-13719
	0.03731993, 0.06607111, 0.08168548, 0.13256945, 0.11239537, 0.14679252, 0.19988304, 0.22475706, //13720-13727
	0.02647217, 0.04693081, 0.05927340, 0.09663250, 0.08116600, 0.10816206, 0.14685876, 0.16649772, //13728-13735
	0.06262555, 0.10237705, 0.11994297, 0.18335974, 0.15994768, 0.19875584, 0.26286133, 0.29004855, //13736-13743
	0.08137091, 0.13681609, 0.15892245, 0.25097980, 0.21923434, 0.27056026, 0.37479391, 0.41193788, //13744-13751
	0.14388287, 0.22342799, 0.24733237, 0.36310999, 0.32830761, 0.38235989, 0.50469689, 0.53965393, //13752-13759
	0.05381513, 0.06371317, 0.09830947, 0.12685408, 0.10848508, 0.16973795, 0.19352876, 0.27312190, //13760-13767
	0.11342243, 0.13082634, 0.17961191, 0.22465821, 0.19681662, 0.27749607, 0.31497137, 0.40906900, //13768-13775
	0.17036319, 0.18868569, 0.23187491, 0.27782962, 0.25011712, 0.32513840, 0.36331663, 0.44997606, //13776-13783
	0.32638842, 0.34942695, 0.39483426, 0.45133564, 0.41897171, 0.49713766, 0.54634722, 0.62791507, //13784-13791
	0.25280441, 0.27593097, 0.32531201, 0.38255618, 0.34900813, 0.43359359, 0.48316992, 0.57449457, //13792-13799
	0.39234994, 0.41215820, 0.44865340, 0.50002026, 0.47104775, 0.53899082, 0.58536536, 0.65760084, //13800-13807
	0.58622613, 0.60405071, 0.63400104, 0.68201159, 0.65612102, 0.71279925, 0.76019212, 0.81518890, //13808-13815
	0.70800253, 0.71831438, 0.73475500, 0.76729994, 0.75015584, 0.78654643, 0.82252339, 0.86022386, //13816-13823
	0.00498860, 0.00878663, 0.01241139, 0.01914704, 0.01597628, 0.02199290, 0.02709364, 0.03028628, //13824-13831
	0.00952485, 0.01766645, 0.02466390, 0.04033404, 0.03317393, 0.04644307, 0.06050531, 0.06886893, //13832-13839
	0.02540554, 0.04258087, 0.05360706, 0.08204256, 0.06991984, 0.09141002, 0.11664955, 0.13033701, //13840-13847
	0.04513583, 0.07769471, 0.09538238, 0.15060022, 0.12871876, 0.16610643, 0.22148729, 0.24744169, //13848-13855
	0.03246214, 0.05597326, 0.07018535, 0.11129839, 0.09425159, 0.12399276, 0.16492429, 0.18572848, //13856-13863
	0.07479803, 0.11894900, 0.13837282, 0.20619069, 0.18112358, 0.22260486, 0.28867129, 0.31663377, //13864-13871
	0.09557443, 0.15608416, 0.18013698, 0.27638359, 0.24318832, 0.29696786, 0.40165787, 0.43902496, //13872-13879
	0.16495286, 0.24932578, 0.27483433, 0.39279912, 0.35737622, 0.41233924, 0.53284862, 0.56720565, //13880-13887
	0.06335838, 0.07459665, 0.11395151, 0.14529002, 0.12515261, 0.19235722, 0.21772031, 0.30219645, //13888-13895
	0.12858470, 0.14765088, 0.20114081, 0.24894239, 0.21940656, 0.30512332, 0.34346050, 0.43983320, //13896-13903
	0.19084161, 0.21059570, 0.25718961, 0.30534983, 0.27628869, 0.35471704, 0.39360707, 0.48124244, //13904-13911
	0.35179224, 0.37574580, 0.42291538, 0.48001643, 0.44726704, 0.52624270, 0.57446990, 0.65443905, //13912-13919
	0.27593140, 0.30024582, 0.35217721, 0.41071554, 0.37640576, 0.46285149, 0.51204917, 0.60248889, //13920-13927
	0.41991088, 0.44019841, 0.47761840, 0.52895921, 0.49996585, 0.56784905, 0.61304405, 0.68330362, //13928-13935
	0.60929972, 0.62707267, 0.65696118, 0.70377840, 0.67852415, 0.73375642, 0.77872357, 0.83099072, //13936-13943
	0.72822290, 0.73833992, 0.75446781, 0.78566370, 0.76921959, 0.80409480, 0.83783548, 0.87315708, //13944-13951
	0.00229534, 0.00443227, 0.00645562, 0.01078575, 0.00875312, 0.01260053, 0.01627314, 0.01856314, //13952-13959
	0.00484210, 0.00986382, 0.01419429, 0.02519681, 0.02016303, 0.02949019, 0.04055790, 0.04713407, //13960-13967
	0.01375619, 0.02509122, 0.03242807, 0.05361083, 0.04453691, 0.06057513, 0.08137088, 0.09265523, //13968-13975
	0.02681888, 0.05043187, 0.06328705, 0.10833840, 0.09046918, 0.12095753, 0.17087065, 0.19428534, //13976-13983
	0.01843176, 0.03479210, 0.04463690, 0.07695782, 0.06357436, 0.08689491, 0.12258631, 0.14071926, //13984-13991
	0.04644736, 0.08018907, 0.09511921, 0.15262324, 0.13142082, 0.16673301, 0.22825921, 0.25429443, //13992-13999
	0.06223667, 0.11105085, 0.13042324, 0.21685581, 0.18708871, 0.23516933, 0.33870081, 0.37555646, //14000-14007
	0.11539124, 0.18863951, 0.21041067, 0.32318433, 0.28937356, 0.34203328, 0.46692330, 0.50256792, //14008-14015
	0.04098190, 0.04906346, 0.07734918, 0.10197118, 0.08614955, 0.13911150, 0.16119447, 0.23448482, //14016-14023
	0.09308477, 0.10823829, 0.15068551, 0.19209039, 0.16646914, 0.24053177, 0.27671745, 0.36776243, //14024-14031
	0.14280430, 0.15919383, 0.19799178, 0.24090130, 0.21495910, 0.28537423, 0.32268827, 0.40753912, //14032-14039
	0.29217049, 0.31409337, 0.35723853, 0.41285424, 0.38093805, 0.45807304, 0.50858216, 0.59224421, //14040-14047
	0.22173556, 0.24325791, 0.28914782, 0.34481861, 0.31221746, 0.39428164, 0.44444766, 0.53692098, //14048-14055
	0.35540181, 0.37455851, 0.40983820, 0.46109314, 0.43208504, 0.49997565, 0.54818705, 0.62306459, //14056-14063
	0.55528133, 0.57314154, 0.60318868, 0.65279740, 0.62603077, 0.68463808, 0.73525370, 0.79405022, //14064-14071
	0.68081593, 0.69140772, 0.70829149, 0.74262687, 0.72453012, 0.76297356, 0.80198538, 0.84287139, //14072-14079
	0.00167417, 0.00328445, 0.00483094, 0.00821898, 0.00661903, 0.00965013, 0.01263909, 0.01451052, //14080-14087
	0.00359931, 0.00749989, 0.01086672, 0.01978972, 0.01571264, 0.02326785, 0.03263138, 0.03820596, //14088-14095
	0.01035486, 0.01934120, 0.02508082, 0.04257088, 0.03512027, 0.04833604, 0.06621752, 0.07591774, //14096-14103
	0.02069114, 0.03998847, 0.05048506, 0.08908656, 0.07377089, 0.09990368, 0.14488069, 0.16595679, //14104-14111
	0.01403040, 0.02721991, 0.03517730, 0.06240495, 0.05113036, 0.07082179, 0.10245855, 0.11850374, //14112-14119
	0.03613932, 0.06424568, 0.07662723, 0.12676228, 0.10827435, 0.13895842, 0.19527176, 0.21919495, //14120-14127
	0.04939987, 0.09084769, 0.10732015, 0.18491448, 0.15817736, 0.20145395, 0.29994819, 0.33508225, //14128-14135
	0.09355999, 0.15764295, 0.17695167, 0.28130334, 0.24995424, 0.29866597, 0.42144642, 0.45655527, //14136-14143
	0.03257966, 0.03917182, 0.06224998, 0.08308093, 0.06965453, 0.11434074, 0.13359983, 0.19782572, //14144-14151
	0.07710753, 0.08998948, 0.12614291, 0.16252596, 0.14004084, 0.20526631, 0.23858582, 0.32225020, //14152-14159
	0.11962204, 0.13377852, 0.16730096, 0.20575826, 0.18245657, 0.24544656, 0.28035612, 0.35957261, //14160-14167
	0.25714207, 0.27707423, 0.31636331, 0.36891659, 0.33875611, 0.41151882, 0.46141962, 0.54419390, //14168-14175
	0.19229656, 0.21156063, 0.25273841, 0.30419470, 0.27405383, 0.35007941, 0.39894630, 0.48877093, //14176-14183
	0.31537168, 0.33298866, 0.36555846, 0.41467631, 0.38701140, 0.45181532, 0.50000332, 0.57487699, //14184-14191
	0.51637016, 0.53368838, 0.56281688, 0.61270737, 0.58576903, 0.64462321, 0.69810743, 0.76014038, //14192-14199
	0.64264046, 0.65318318, 0.66997996, 0.70549698, 0.68679754, 0.72648069, 0.76874695, 0.81302053, //14200-14207
	0.00069958, 0.00151100, 0.00228062, 0.00422630, 0.00330971, 0.00504953, 0.00699109, 0.00820632, //14208-14215
	0.00166928, 0.00383138, 0.00569674, 0.01136837, 0.00877627, 0.01358020, 0.02029229, 0.02428707, //14216-14223
	0.00505942, 0.01033734, 0.01373685, 0.02540949, 0.02044631, 0.02929270, 0.04264300, 0.04988388, //14224-14231
	0.01114334, 0.02374227, 0.03055461, 0.05910537, 0.04780719, 0.06714719, 0.10443831, 0.12191769, //14232-14239
	0.00722714, 0.01546037, 0.02041908, 0.03980717, 0.03178027, 0.04578252, 0.07109891, 0.08396628, //14240-14247
	0.02026555, 0.03949007, 0.04788466, 0.08647436, 0.07216982, 0.09584396, 0.14397433, 0.16443867, //14248-14255
	0.02928253, 0.05946826, 0.07138690, 0.13523938, 0.11319592, 0.14882141, 0.23968252, 0.27213260, //14256-14263
	0.05941960, 0.10962981, 0.12477740, 0.21609302, 0.18862507, 0.23129943, 0.35071977, 0.38484769, //14264-14271
	0.01949938, 0.02377401, 0.03877617, 0.05361069, 0.04405473, 0.07573045, 0.09081388, 0.14113723, //14272-14279
	0.05223782, 0.06162988, 0.08793653, 0.11669504, 0.09892934, 0.15036067, 0.17910617, 0.25147569, //14280-14287
	0.08348301, 0.09424284, 0.11950225, 0.15100732, 0.13203542, 0.18329206, 0.21436897, 0.28469120, //14288-14295
	0.20255510, 0.21953182, 0.25281887, 0.30048753, 0.27316432, 0.33904775, 0.38812013, 0.46931705, //14296-14303
	0.14645188, 0.16225367, 0.19595502, 0.24118221, 0.21468509, 0.28143468, 0.32806431, 0.41388242, //14304-14311
	0.25297351, 0.26834274, 0.29681835, 0.34244286, 0.31673378, 0.37697664, 0.42517081, 0.49998517, //14312-14319
	0.45587052, 0.47230367, 0.49999315, 0.55029250, 0.52316794, 0.58251638, 0.64029470, 0.70734546, //14320-14327
	0.58323330, 0.59373552, 0.61045047, 0.64768697, 0.62805315, 0.66972810, 0.71707300, 0.76664685, //14328-14335
	0.00550297, 0.00899773, 0.01233437, 0.01795258, 0.01530340, 0.02032895, 0.02431240, 0.02681673, //14336-14343
	0.00967719, 0.01667317, 0.02270380, 0.03511250, 0.02943613, 0.03995077, 0.05063967, 0.05700475, //14344-14351
	0.02429161, 0.03820627, 0.04714239, 0.06888379, 0.05960940, 0.07605561, 0.09485262, 0.10507624, //14352-14359
	0.04029366, 0.06552752, 0.07923550, 0.12050371, 0.10413584, 0.13207637, 0.17384601, 0.19342951, //14360-14367
	0.03000349, 0.04863528, 0.05988908, 0.09091236, 0.07805451, 0.10047863, 0.13123570, 0.14686111, //14368-14375
	0.06427035, 0.09759471, 0.11229179, 0.16233680, 0.14382889, 0.17447672, 0.22443824, 0.24563266, //14376-14383
	0.07998289, 0.12514661, 0.14311961, 0.21520153, 0.19033613, 0.23055150, 0.31422097, 0.34405332, //14384-14391
	0.13228793, 0.19438177, 0.21304527, 0.30285098, 0.27590964, 0.31779052, 0.41841035, 0.44716459, //14392-14399
	0.05402948, 0.06283822, 0.09366823, 0.11720295, 0.10207601, 0.15246019, 0.17107480, 0.23324583, //14400-14407
	0.10375142, 0.11824992, 0.15892462, 0.19459815, 0.17254584, 0.23645162, 0.26548975, 0.33856150, //14408-14415
	0.15120010, 0.16600190, 0.20087393, 0.23668099, 0.21506975, 0.27347176, 0.30311517, 0.37017273, //14416-14423
	0.27577339, 0.29412304, 0.33019334, 0.37499257, 0.34931195, 0.41128796, 0.45194321, 0.51931919, //14424-14431
	0.21704757, 0.23553583, 0.27500574, 0.32013233, 0.29370992, 0.36036110, 0.40056530, 0.47457632, //14432-14439
	0.32846602, 0.34412585, 0.37291689, 0.41378858, 0.39074022, 0.44474017, 0.48365488, 0.54418188, //14440-14447
	0.50001130, 0.51492576, 0.53997011, 0.58183531, 0.55925845, 0.60864751, 0.65461150, 0.70795910, //14448-14455
	0.60773970, 0.61663310, 0.63079354, 0.66058231, 0.64487748, 0.67819891, 0.71557716, 0.75475992, //14456-14463
	0.00441271, 0.00736733, 0.01018700, 0.01511312, 0.01279266, 0.01720037, 0.02081371, 0.02308204, //14464-14471
	0.00794244, 0.01400836, 0.01923715, 0.03040822, 0.02529897, 0.03477018, 0.04471636, 0.05063755, //14472-14479
	0.02028847, 0.03263268, 0.04057365, 0.06050183, 0.05199921, 0.06710355, 0.08486243, 0.09451957, //14480-14487
	0.03448468, 0.05745202, 0.06993689, 0.10881940, 0.09340572, 0.11972064, 0.16023443, 0.17922342, //14488-14495
	0.02536100, 0.04210790, 0.05221923, 0.08112960, 0.06914653, 0.09005581, 0.11957897, 0.13458389, //14496-14503
	0.05577482, 0.08658480, 0.10018644, 0.14796494, 0.13031746, 0.15957722, 0.20850375, 0.22926668, //14504-14511
	0.07028994, 0.11271672, 0.12959528, 0.19947390, 0.17535700, 0.21436516, 0.29759086, 0.32726590, //14512-14519
	0.11868057, 0.17838874, 0.19626420, 0.28484688, 0.25826353, 0.29960441, 0.40097104, 0.42993797, //14520-14527
	0.04725354, 0.05524324, 0.08321757, 0.10516073, 0.09104105, 0.13807322, 0.15597658, 0.21563485, //14528-14535
	0.09373802, 0.10728013, 0.14525989, 0.17946490, 0.15830799, 0.21959317, 0.24814997, 0.31998018, //14536-14543
	0.13810264, 0.15212100, 0.18515311, 0.21982978, 0.19891040, 0.25554563, 0.28482139, 0.35124371, //14544-14551
	0.25993787, 0.27779286, 0.31294138, 0.35738401, 0.33191833, 0.39340216, 0.43444697, 0.50248479, //14552-14559
	0.20250421, 0.22033640, 0.25841038, 0.30283334, 0.27680996, 0.34243777, 0.38283366, 0.45711496, //14560-14567
	0.31148758, 0.32686138, 0.35513830, 0.39596822, 0.37295109, 0.42688584, 0.46632528, 0.52769084, //14568-14575
	0.48510279, 0.50001526, 0.52505611, 0.56747991, 0.54460228, 0.59464473, 0.64180269, 0.69652116, //14576-14583
	0.59413998, 0.60312436, 0.61742769, 0.64785024, 0.63182190, 0.66582489, 0.70437511, 0.74475745, //14584-14591
	0.00258510, 0.00462901, 0.00658063, 0.01035243, 0.00857496, 0.01194532, 0.01493337, 0.01681076, //14592-14599
	0.00502607, 0.00952987, 0.01341332, 0.02251512, 0.01835026, 0.02606314, 0.03476815, 0.03995055, //14600-14607
	0.01358656, 0.02327168, 0.02952350, 0.04645643, 0.03922307, 0.05203888, 0.06808107, 0.07680214, //14608-14615
	0.02472266, 0.04388909, 0.05432490, 0.08919209, 0.07536558, 0.09899038, 0.13736816, 0.15536049, //14616-14623
	0.01757016, 0.03114285, 0.03935235, 0.06469361, 0.05419528, 0.07250390, 0.09998863, 0.11396311, //14624-14631
	0.04143961, 0.06810432, 0.07986790, 0.12385079, 0.10761712, 0.13455815, 0.18174940, 0.20180667, //14632-14639
	0.05398931, 0.09181110, 0.10686957, 0.17305464, 0.15021185, 0.18713931, 0.26966550, 0.29909539, //14640-14647
	0.09586834, 0.15140396, 0.16804053, 0.25461528, 0.22864192, 0.26898273, 0.37167313, 0.40101799, //14648-14655
	0.03588064, 0.04248444, 0.06558945, 0.08492033, 0.07248941, 0.11390328, 0.13057634, 0.18620171, //14656-14663
	0.07692933, 0.08885550, 0.12229118, 0.15403535, 0.13442481, 0.19123998, 0.21899383, 0.28876990, //14664-14671
	0.11608887, 0.12879673, 0.15879517, 0.19157762, 0.17178143, 0.22535087, 0.25420528, 0.31957152, //14672-14679
	0.23333936, 0.25039821, 0.28393342, 0.32781775, 0.30266697, 0.36332297, 0.40507522, 0.47425093, //14680-14687
	0.17806693, 0.19482404, 0.23054565, 0.27379951, 0.24845787, 0.31232276, 0.35300752, 0.42778932, //14688-14695
	0.28293294, 0.29787228, 0.32530842, 0.36601813, 0.34306499, 0.39682877, 0.43722886, 0.50003476, //14696-14703
	0.46005596, 0.47497095, 0.49999685, 0.54337868, 0.51999739, 0.57114078, 0.62030632, 0.67731069, //14704-14711
	0.57127869, 0.58042865, 0.59497023, 0.62647499, 0.60987255, 0.64508943, 0.68555164, 0.72798766, //14712-14719
	0.00155562, 0.00289849, 0.00418042, 0.00683278, 0.00558250, 0.00795853, 0.01020470, 0.01161794, //14720-14727
	0.00316029, 0.00626703, 0.00894520, 0.01571121, 0.01261694, 0.01835201, 0.02530014, 0.02943550, //14728-14735
	0.00877450, 0.01571340, 0.02018166, 0.03317324, 0.02762747, 0.03747202, 0.05065043, 0.05780955, //14736-14743
	0.01675093, 0.03125223, 0.03913096, 0.06759319, 0.05630610, 0.07556446, 0.10924627, 0.12502970, //14744-14751
	0.01162426, 0.02164570, 0.02768084, 0.04786943, 0.03950531, 0.05409092, 0.07763739, 0.08959787, //14752-14759
	0.02870689, 0.04948998, 0.05867383, 0.09548211, 0.08185868, 0.10442146, 0.14676927, 0.16474358, //14760-14767
	0.03849223, 0.06905047, 0.08122574, 0.13901880, 0.11906746, 0.15134027, 0.22903541, 0.25674672, //14768-14775
	0.07113777, 0.11822157, 0.13236139, 0.21118675, 0.18754458, 0.22432881, 0.32450384, 0.35312541, //14776-14783
	0.02556194, 0.03053829, 0.04795644, 0.06334433, 0.05345280, 0.08642160, 0.10055581, 0.14759315, //14784-14791
	0.05879258, 0.06840056, 0.09535257, 0.12229752, 0.10563573, 0.15393717, 0.17899664, 0.24200526, //14792-14799
	0.09050546, 0.10102030, 0.12582662, 0.15435665, 0.13712980, 0.18372584, 0.21023110, 0.27034854, //14800-14807
	0.19623367, 0.21142590, 0.24129493, 0.28214840, 0.25874642, 0.31523309, 0.35623229, 0.42410026, //14808-14815
	0.14639028, 0.16095016, 0.19201296, 0.23148293, 0.20836357, 0.26664419, 0.30584907, 0.37795546, //14816-14823
	0.24097189, 0.25451754, 0.27944067, 0.31801906, 0.29626567, 0.34719043, 0.38733006, 0.44975817, //14824-14831
	0.41819407, 0.43255508, 0.45665839, 0.50001369, 0.47661933, 0.52778200, 0.57915240, 0.63875364, //14832-14839
	0.52949765, 0.53856157, 0.55300391, 0.58530263, 0.56828308, 0.60438502, 0.64754971, 0.69278051, //14840-14847
	0.00210957, 0.00383333, 0.00547766, 0.00873101, 0.00719593, 0.01010893, 0.01275533, 0.01441822, //14848-14855
	0.00416680, 0.00802583, 0.01135417, 0.01938227, 0.01571033, 0.02251411, 0.03040805, 0.03510755, //14856-14863
	0.01136988, 0.01979604, 0.02520666, 0.04034137, 0.03388406, 0.04534041, 0.06005792, 0.06805015, //14864-14871
	0.02105660, 0.03807706, 0.04732509, 0.07925132, 0.06658474, 0.08820168, 0.12441815, 0.14139221, //14872-14879
	0.01483161, 0.02676866, 0.03397671, 0.05694590, 0.04743184, 0.06403580, 0.08969867, 0.10274603, //14880-14887
	0.03556463, 0.05952131, 0.07011726, 0.11077181, 0.09573992, 0.12065849, 0.16563186, 0.18472305, //14888-14895
	0.04687642, 0.08133474, 0.09505115, 0.15737411, 0.13588210, 0.17065095, 0.25095578, 0.27957973, //14896-14903
	0.08448684, 0.13615259, 0.15166289, 0.23460767, 0.20971172, 0.24844898, 0.34995851, 0.37896783, //14904-14911
	0.03112900, 0.03698744, 0.05746592, 0.07499957, 0.06373408, 0.10124396, 0.11675638, 0.16844367, //14912-14919
	0.06857267, 0.07943171, 0.10989924, 0.13941952, 0.12116387, 0.17406360, 0.20058022, 0.26721156, //14920-14927
	0.10431783, 0.11601143, 0.14358626, 0.17444710, 0.15582437, 0.20616918, 0.23395398, 0.29687367, //14928-14935
	0.21625798, 0.23244864, 0.26430348, 0.30678013, 0.28244078, 0.34120159, 0.38257659, 0.45117231, //14936-14943
	0.16348141, 0.17921002, 0.21278052, 0.25430413, 0.22998938, 0.29131979, 0.33128126, 0.40486807, //14944-14951
	0.26360878, 0.27789104, 0.30417423, 0.34391631, 0.32149852, 0.37403315, 0.41423925, 0.47687738, //14952-14959
	0.44077060, 0.45543408, 0.48003848, 0.52340642, 0.50001093, 0.55118284, 0.60133472, 0.65956128, //14960-14967
	0.55204120, 0.56114458, 0.57564823, 0.60750768, 0.59071581, 0.62634975, 0.66805833, 0.71175932, //14968-14975
	0.00089843, 0.00179047, 0.00264471, 0.00458121, 0.00366936, 0.00539921, 0.00717477, 0.00829166, //14976-14983
	0.00196551, 0.00417864, 0.00608359, 0.01135364, 0.00894208, 0.01341344, 0.01923382, 0.02270070, //14984-14991
	0.00569891, 0.01087144, 0.01417979, 0.02468959, 0.02020501, 0.02814788, 0.03948105, 0.04564862, //14992-14999
	0.01164445, 0.02314352, 0.02940720, 0.05374814, 0.04409823, 0.06058127, 0.09123366, 0.10559283, //15000-15007
	0.00782059, 0.01555770, 0.02022859, 0.03709187, 0.03010390, 0.04229149, 0.06332168, 0.07400687, //15008-15015
	0.02056302, 0.03756626, 0.04502963, 0.07731709, 0.06539250, 0.08513926, 0.12436334, 0.14102082, //15016-15023
	0.02858867, 0.05446766, 0.06474968, 0.11722135, 0.09912010, 0.12840993, 0.20301736, 0.22961438, //15024-15031
	0.05528170, 0.09694634, 0.10948028, 0.18336316, 0.16118871, 0.19569457, 0.29429431, 0.32245518, //15032-15039
	0.01896754, 0.02289815, 0.03663571, 0.04951181, 0.04124094, 0.06884753, 0.08131210, 0.12287065, //15040-15047
	0.04717195, 0.05529286, 0.07809887, 0.10197271, 0.08719650, 0.13000672, 0.15340074, 0.21206058, //15048-15055
	0.07410312, 0.08321859, 0.10474044, 0.13046400, 0.11494178, 0.15699815, 0.18211397, 0.23879562, //15056-15063
	0.17247607, 0.18645233, 0.21394806, 0.25291066, 0.23057630, 0.28446474, 0.32493251, 0.39197629, //15064-15071
	0.12610872, 0.13925045, 0.16733937, 0.20435721, 0.18266412, 0.23733358, 0.27562495, 0.34605556, //15072-15079
	0.21408124, 0.22676123, 0.25003183, 0.28725734, 0.26628766, 0.31538548, 0.35535815, 0.41752489, //15080-15087
	0.39139225, 0.40538563, 0.42886673, 0.47224396, 0.44886307, 0.50001734, 0.55280329, 0.61407307, //15088-15095
	0.50272765, 0.51174785, 0.52611376, 0.55894143, 0.54164231, 0.57834549, 0.62322529, 0.67021327, //15096-15103
	0.00058250, 0.00118654, 0.00176350, 0.00312897, 0.00248562, 0.00370666, 0.00501267, 0.00583527, //15104-15111
	0.00130461, 0.00284795, 0.00417767, 0.00802723, 0.00626695, 0.00953141, 0.01400735, 0.01667447, //15112-15119
	0.00383029, 0.00750320, 0.00986294, 0.01766489, 0.01433738, 0.02024816, 0.02911480, 0.03393667, //15120-15127
	0.00805608, 0.01653936, 0.02115108, 0.04010169, 0.03258696, 0.04540674, 0.07075640, 0.08263701, //15128-15135
	0.00534012, 0.01095249, 0.01434190, 0.02724369, 0.02189988, 0.03121728, 0.04829420, 0.05696926, //15136-15143
	0.01438951, 0.02714982, 0.03278013, 0.05828802, 0.04885109, 0.06448683, 0.09742925, 0.11140661, //15144-15151
	0.02039963, 0.04038123, 0.04833386, 0.09145738, 0.07657336, 0.10064044, 0.16667948, 0.19022692, //15152-15159
	0.04043007, 0.07364767, 0.08361310, 0.14644501, 0.12760106, 0.15691669, 0.24704901, 0.27279974, //15160-15167
	0.01361820, 0.01651572, 0.02666054, 0.03648716, 0.03017172, 0.05121695, 0.06118790, 0.09443100, //15168-15175
	0.03563972, 0.04194251, 0.05962351, 0.07885291, 0.06695996, 0.10142893, 0.12114129, 0.17070222, //15176-15183
	0.05665975, 0.06384603, 0.08079877, 0.10188561, 0.08916269, 0.12358017, 0.14497686, 0.19358799, //15184-15191
	0.14063824, 0.15242292, 0.17562549, 0.20962342, 0.19013592, 0.23717641, 0.27435946, 0.33598223, //15192-15199
	0.10105243, 0.11192832, 0.13515485, 0.16690470, 0.14830423, 0.19520482, 0.22972908, 0.29327452, //15200-15207
	0.17617127, 0.18700544, 0.20692071, 0.23983569, 0.22127342, 0.26478820, 0.30191263, 0.35971995, //15208-15215
	0.34542088, 0.35823328, 0.37975071, 0.42088817, 0.39868918, 0.44723417, 0.50001836, 0.56129803, //15216-15223
	0.45171123, 0.46019411, 0.47369565, 0.50571491, 0.48884029, 0.52464657, 0.57067024, 0.61887727, //15224-15231
	0.00021555, 0.00048594, 0.00074385, 0.00144498, 0.00111404, 0.00173921, 0.00250497, 0.00298351, //15232-15239
	0.00053811, 0.00130387, 0.00196514, 0.00416852, 0.00316037, 0.00502782, 0.00794142, 0.00967800, //15240-15247
	0.00166868, 0.00359858, 0.00484505, 0.00953071, 0.00752745, 0.01107280, 0.01707367, 0.02034200, //15248-15255
	0.00388668, 0.00886592, 0.01157081, 0.02425185, 0.01922220, 0.02781244, 0.04699687, 0.05598679, //15256-15263
	0.00246150, 0.00560883, 0.00751514, 0.01581323, 0.01237394, 0.01838068, 0.03085488, 0.03720086, //15264-15271
	0.00722574, 0.01504000, 0.01849653, 0.03619794, 0.02964766, 0.04047326, 0.06616878, 0.07704820, //15272-15279
	0.01091129, 0.02403274, 0.02923994, 0.06155668, 0.05041070, 0.06844886, 0.12450939, 0.14451018, //15280-15287
	0.02319410, 0.04665184, 0.05363738, 0.10358827, 0.08859740, 0.11190825, 0.19222566, 0.21516934, //15288-15295
	0.00741602, 0.00912084, 0.01508542, 0.02136714, 0.01733579, 0.03079486, 0.03780908, 0.06128841, //15296-15303
	0.02225043, 0.02643978, 0.03819793, 0.05201277, 0.04347937, 0.06822539, 0.08370474, 0.12260199, //15304-15311
	0.03639961, 0.04137227, 0.05308069, 0.06866502, 0.05925864, 0.08474221, 0.10196876, 0.14102269, //15312-15319
	0.10369992, 0.11292991, 0.13107721, 0.15937975, 0.14317109, 0.18229030, 0.21567744, 0.27096526, //15320-15327
	0.07196959, 0.08020646, 0.09780147, 0.12340618, 0.10840672, 0.14625587, 0.17644244, 0.23198789, //15328-15335
	0.13217568, 0.14085808, 0.15684706, 0.18481111, 0.16906322, 0.20599478, 0.23987835, 0.29260303, //15336-15343
	0.29206624, 0.30351188, 0.32271214, 0.36128727, 0.34047910, 0.38594515, 0.43875462, 0.50003493, //15344-15351
	0.39250733, 0.40036334, 0.41288853, 0.44395384, 0.42757583, 0.46230218, 0.50967123, 0.55930435, //15352-15359
	0.00149105, 0.00274226, 0.00393145, 0.00635454, 0.00521234, 0.00737998, 0.00941433, 0.01068964, //15360-15367
	0.00298250, 0.00583416, 0.00829119, 0.01441993, 0.01161811, 0.01681129, 0.02308315, 0.02681821, //15368-15375
	0.00820517, 0.01451142, 0.01856022, 0.03029224, 0.02528038, 0.03414430, 0.04603487, 0.05250392, //15376-15383
	0.01545920, 0.02855404, 0.03567017, 0.06131397, 0.05114414, 0.06851095, 0.09917718, 0.11354270, //15384-15391
	0.01079181, 0.01987228, 0.02534532, 0.04353236, 0.03600023, 0.04914770, 0.07051292, 0.08137213, //15392-15399
	0.02634540, 0.04502280, 0.05326397, 0.08645929, 0.07419816, 0.09449244, 0.13321597, 0.14962700, //15400-15407
	0.03514501, 0.06265448, 0.07362407, 0.12614291, 0.10801790, 0.13733358, 0.20972341, 0.23554095, //15408-15415
	0.06447624, 0.10695581, 0.11970061, 0.19189861, 0.17022908, 0.20392864, 0.29852805, 0.32555836, //15416-15423
	0.02339827, 0.02789871, 0.04365111, 0.05751151, 0.04860414, 0.07830149, 0.09103314, 0.13347427, //15424-15431
	0.05346670, 0.06215863, 0.08652138, 0.11092614, 0.09584402, 0.13955268, 0.16247582, 0.22007892, //15432-15439
	0.08217457, 0.09167410, 0.11406342, 0.13993621, 0.12433903, 0.16656682, 0.19092804, 0.24606343, //15440-15447
	0.17954469, 0.19347516, 0.22085120, 0.25872045, 0.23702208, 0.28939541, 0.32816769, 0.39240096, //15448-15455
	0.13364572, 0.14692664, 0.17529930, 0.21161561, 0.19035336, 0.24398859, 0.28079718, 0.34849464, //15456-15463
	0.22074570, 0.23320585, 0.25612908, 0.29202915, 0.27180054, 0.31923121, 0.35737949, 0.41678636, //15464-15471
	0.39228718, 0.40589308, 0.42874759, 0.47054322, 0.44799069, 0.49731113, 0.54832867, 0.60753679, //15472-15479
	0.50002630, 0.50871991, 0.52256948, 0.55418569, 0.53752727, 0.57284637, 0.61641815, 0.66204026, //15480-15487
	0.00122346, 0.00229532, 0.00332111, 0.00547212, 0.00445749, 0.00638051, 0.00824418, 0.00941264, //15488-15495
	0.00250488, 0.00501466, 0.00717656, 0.01275449, 0.01020507, 0.01493473, 0.02081311, 0.02431276, //15496-15503
	0.00698874, 0.01263226, 0.01626980, 0.02709683, 0.02247382, 0.03067544, 0.04192442, 0.04803814, //15504-15511
	0.01348028, 0.02549558, 0.03202770, 0.05626271, 0.04665251, 0.06307165, 0.09276815, 0.10667965, //15512-15519
	0.00930947, 0.01754614, 0.02252119, 0.03954572, 0.03249247, 0.04479879, 0.06533908, 0.07577731, //15520-15527
	0.02323014, 0.04059125, 0.04825283, 0.07990396, 0.06822109, 0.08760569, 0.12532369, 0.14133953, //15528-15535
	0.03139867, 0.05731490, 0.06761925, 0.11844438, 0.10090973, 0.12927407, 0.20080306, 0.22631225, //15536-15543
	0.05865559, 0.09935829, 0.11153609, 0.18229342, 0.16106867, 0.19408077, 0.28840155, 0.31534622, //15544-15551
	0.02090033, 0.02502129, 0.03945250, 0.05241817, 0.04409055, 0.07188353, 0.08409158, 0.12474655, //15552-15559
	0.04922546, 0.05739508, 0.08030332, 0.10370177, 0.08924062, 0.13116370, 0.15355556, 0.20983626, //15560-15567
	0.07627016, 0.08530031, 0.10657630, 0.13157464, 0.11649530, 0.15730674, 0.18118747, 0.23536009, //15568-15575
	0.17133820, 0.18487943, 0.21153815, 0.24882993, 0.22744119, 0.27906178, 0.31775829, 0.38188415, //15576-15583
	0.12652351, 0.13934482, 0.16671986, 0.20233124, 0.18146847, 0.23405273, 0.27061633, 0.33794351, //15584-15591
	0.21156645, 0.22375158, 0.24618563, 0.28171295, 0.26166909, 0.30861579, 0.34685155, 0.40628973, //15592-15599
	0.38339230, 0.39690598, 0.41960243, 0.46147590, 0.43888124, 0.48830102, 0.53984876, 0.59967051, //15600-15607
	0.49130885, 0.50001227, 0.51387442, 0.54569621, 0.52892747, 0.56451286, 0.60867360, 0.65493373, //15608-15615
	0.00079483, 0.00158675, 0.00234086, 0.00406420, 0.00325122, 0.00479214, 0.00637996, 0.00737664, //15616-15623
	0.00174229, 0.00370537, 0.00539751, 0.01010836, 0.00795594, 0.01194696, 0.01719782, 0.02033006, //15624-15631
	0.00504894, 0.00964192, 0.01261279, 0.02200329, 0.01798451, 0.02509710, 0.03536863, 0.04095295, //15632-15639
	0.01033753, 0.02063740, 0.02622168, 0.04823186, 0.03950529, 0.05439651, 0.08256159, 0.09575299, //15640-15647
	0.00694456, 0.01383704, 0.01800967, 0.03319805, 0.02690163, 0.03787170, 0.05710589, 0.06687133, //15648-15655
	0.01827753, 0.03351343, 0.04027313, 0.06950229, 0.05865818, 0.07658274, 0.11278678, 0.12814029, //15656-15663
	0.02547030, 0.04881890, 0.05814819, 0.10620210, 0.08960535, 0.11645893, 0.18660643, 0.21161770, //15664-15671
	0.04940395, 0.08722276, 0.09861758, 0.16704045, 0.14649958, 0.17844979, 0.27230921, 0.29909183, //15672-15679
	0.01692053, 0.02043539, 0.03277871, 0.04432747, 0.03690138, 0.06171849, 0.07306895, 0.11090154, //15680-15687
	0.04247482, 0.04982331, 0.07039874, 0.09221107, 0.07874409, 0.11782841, 0.13938425, 0.19354836, //15688-15695
	0.06686971, 0.07513903, 0.09460044, 0.11825072, 0.10400593, 0.14251199, 0.16570461, 0.21837526, //15696-15703
	0.15825973, 0.17119140, 0.19663512, 0.23307215, 0.21218172, 0.26260965, 0.30117782, 0.36509292, //15704-15711
	0.11516643, 0.12726377, 0.15309396, 0.18753175, 0.16735067, 0.21821573, 0.25445073, 0.32109377, //15712-15719
	0.19694424, 0.20868917, 0.23032391, 0.26526343, 0.24556958, 0.29178142, 0.33004299, 0.38957562, //15720-15727
	0.36922192, 0.38260931, 0.40504082, 0.44702951, 0.42438838, 0.47391817, 0.52634114, 0.58714805, //15728-15735
	0.47747206, 0.48616987, 0.50004503, 0.53222574, 0.51524873, 0.55124618, 0.59634712, 0.64360078, //15736-15743
	0.00048579, 0.00101152, 0.00151108, 0.00273748, 0.00216057, 0.00325315, 0.00445728, 0.00521295, //15744-15751
	0.00111393, 0.00248444, 0.00366872, 0.00719481, 0.00558255, 0.00857503, 0.01279185, 0.01530342, //15752-15759
	0.00331234, 0.00662042, 0.00874894, 0.01597203, 0.01288596, 0.01836653, 0.02678728, 0.03137511, //15760-15767
	0.00711237, 0.01493625, 0.01917932, 0.03714736, 0.03002456, 0.04218679, 0.06677223, 0.07828793, //15768-15775
	0.00467127, 0.00979440, 0.01289089, 0.02501338, 0.01999088, 0.02874919, 0.04520425, 0.05356078, //15776-15783
	0.01283092, 0.02469145, 0.02993791, 0.05431684, 0.04530013, 0.06023333, 0.09236360, 0.10601864, //15784-15791
	0.01842068, 0.03727848, 0.04480542, 0.08656761, 0.07214935, 0.09546898, 0.16073849, 0.18402237, //15792-15799
	0.03710308, 0.06892629, 0.07847258, 0.14001500, 0.12154694, 0.15029089, 0.24009038, 0.26572353, //15800-15807
	0.01232838, 0.01500007, 0.02435797, 0.03358026, 0.02764395, 0.04742578, 0.05696514, 0.08885740, //15808-15815
	0.03318640, 0.03914751, 0.05586460, 0.07432644, 0.06291185, 0.09601710, 0.11528429, 0.16373197, //15816-15823
	0.05310523, 0.05994856, 0.07610977, 0.09648819, 0.08418722, 0.11745585, 0.13846944, 0.18620779, //15824-15831
	0.13526480, 0.14674505, 0.16935853, 0.20287582, 0.18366391, 0.23004877, 0.26712885, 0.32859165, //15832-15839
	0.09653408, 0.10707083, 0.12955048, 0.16068441, 0.14245725, 0.18845624, 0.22275737, 0.28587924, //15840-15847
	0.17003990, 0.18062232, 0.20016386, 0.23271535, 0.21433365, 0.25736631, 0.29454741, 0.35235855, //15848-15855
	0.33945186, 0.35217866, 0.37357093, 0.41472393, 0.39250529, 0.44111938, 0.49432839, 0.55609699, //15856-15863
	0.44585783, 0.45433481, 0.46784875, 0.50003100, 0.48305330, 0.51902947, 0.56556996, 0.61431601, //15864-15871
	0.00064937, 0.00131333, 0.00195044, 0.00343786, 0.00273655, 0.00406604, 0.00547280, 0.00635426, //15872-15879
	0.00144331, 0.00312914, 0.00458085, 0.00872979, 0.00683228, 0.01035257, 0.01511399, 0.01795236, //15880-15887
	0.00422750, 0.00821588, 0.01077807, 0.01915171, 0.01557823, 0.02191991, 0.03130837, 0.03641736, //15888-15895
	0.00881782, 0.01793979, 0.02289686, 0.04298060, 0.03502043, 0.04862351, 0.07508887, 0.08749943, //15896-15903
	0.00586578, 0.01192639, 0.01558529, 0.02932916, 0.02363177, 0.03356000, 0.05147726, 0.06057437, //15904-15911
	0.01570212, 0.02935109, 0.03539096, 0.06229973, 0.05233796, 0.06884486, 0.10312907, 0.11768452, //15912-15919
	0.02212649, 0.04336578, 0.05180780, 0.09691282, 0.08134833, 0.10651311, 0.17437641, 0.19856628, //15920-15927
	0.04358510, 0.07857956, 0.08910551, 0.15426177, 0.13469715, 0.16512667, 0.25704817, 0.28331551, //15928-15935
	0.01475200, 0.01786499, 0.02877136, 0.03923818, 0.03251949, 0.05495239, 0.06545878, 0.10045222, //15936-15943
	0.03808113, 0.04476768, 0.06353013, 0.08375622, 0.07124206, 0.10747749, 0.12796669, 0.17944992, //15944-15951
	0.06035893, 0.06795564, 0.08586654, 0.10793692, 0.09462343, 0.13067627, 0.15286750, 0.20312849, //15952-15959
	0.14737323, 0.15963470, 0.18372779, 0.21879168, 0.19869836, 0.24718955, 0.28507274, 0.34782962, //15960-15967
	0.10635119, 0.11771496, 0.14195283, 0.17482529, 0.15558171, 0.20412899, 0.23943981, 0.30442255, //15968-15975
	0.18420928, 0.19540725, 0.21606934, 0.24987501, 0.23078905, 0.27541509, 0.31324505, 0.37199358, //15976-15983
	0.35514808, 0.36821214, 0.39016071, 0.43175865, 0.40932144, 0.45839925, 0.51118815, 0.57245477, //15984-15991
	0.46252345, 0.47110430, 0.48480178, 0.51699297, 0.50001689, 0.53600693, 0.58179461, 0.62974100, //15992-15999
	0.00030270, 0.00067204, 0.00102374, 0.00195015, 0.00151303, 0.00234239, 0.00331907, 0.00393152, //16000-16007
	0.00074453, 0.00176497, 0.00264440, 0.00547609, 0.00418040, 0.00657755, 0.01018831, 0.01233416, //16008-16015
	0.00228054, 0.00482897, 0.00646463, 0.01240918, 0.00988094, 0.01436673, 0.02172123, 0.02571624, //16016-16023
	0.00521269, 0.01156936, 0.01501599, 0.03058955, 0.02441965, 0.03497539, 0.05743818, 0.06797023, //16024-16031
	0.00333295, 0.00740032, 0.00986573, 0.02018406, 0.01589618, 0.02335521, 0.03817487, 0.04570427, //16032-16039
	0.00960843, 0.01948808, 0.02386359, 0.04532358, 0.03737162, 0.05051827, 0.08032257, 0.09295290, //16040-16047
	0.01426815, 0.03045411, 0.03693702, 0.07495911, 0.06183271, 0.08305816, 0.14545526, 0.16770640, //16048-16055
	0.02976654, 0.05809053, 0.06650624, 0.12405118, 0.10679197, 0.13362534, 0.22104099, 0.24602369, //16056-16063
	0.00961279, 0.01179053, 0.01940774, 0.02719808, 0.02218888, 0.03898881, 0.04747806, 0.07577256, //16064-16071
	0.02770266, 0.03283688, 0.04725151, 0.06377513, 0.05356778, 0.08313176, 0.10103305, 0.14606425, //16072-16079
	0.04496374, 0.05098560, 0.06515756, 0.08365974, 0.07250086, 0.10267532, 0.12241935, 0.16719042, //16080-16087
	0.12168615, 0.13229767, 0.15322866, 0.18503109, 0.16677826, 0.21083298, 0.24699943, 0.30703160, //16088-16095
	0.08552719, 0.09511681, 0.11563455, 0.14481945, 0.12772808, 0.17088319, 0.20399724, 0.26504429, //16096-16103
	0.15412184, 0.16405130, 0.18227827, 0.21347742, 0.19587038, 0.23713395, 0.27358996, 0.33039433, //16104-16111
	0.32185423, 0.33420017, 0.35494371, 0.39564121, 0.37369282, 0.42167136, 0.47540022, 0.53772671, //16112-16119
	0.42718231, 0.43552089, 0.44881663, 0.48100010, 0.46403288, 0.49999813, 0.54738378, 0.59697490, //16120-16127
	0.00019357, 0.00043833, 0.00067195, 0.00131502, 0.00101152, 0.00158754, 0.00229573, 0.00274040, //16128-16135
	0.00048563, 0.00118679, 0.00179223, 0.00383127, 0.00289874, 0.00462595, 0.00736816, 0.00899728, //16136-16143
	0.00150881, 0.00328835, 0.00442977, 0.00878264, 0.00692596, 0.01022736, 0.01588698, 0.01895688, //16144-16151
	0.00355275, 0.00817399, 0.01068374, 0.02263288, 0.01789156, 0.02599099, 0.04432807, 0.05291671, //16152-16159
	0.00223950, 0.00515284, 0.00691122, 0.01469787, 0.01147208, 0.01710133, 0.02898352, 0.03502661, //16160-16167
	0.00661439, 0.01390265, 0.01711902, 0.03384989, 0.02767126, 0.03792225, 0.06255095, 0.07298602, //16168-16175
	0.01005092, 0.02238609, 0.02729781, 0.05812524, 0.04748644, 0.06469190, 0.11915317, 0.13857142, //16176-16183
	0.02148823, 0.04369335, 0.05033598, 0.09835480, 0.08395673, 0.10638032, 0.18485495, 0.20726984, //16184-16191
	0.00685482, 0.00843910, 0.01397969, 0.01987680, 0.01608580, 0.02868679, 0.03536250, 0.05752936, //16192-16199
	0.02085255, 0.02480709, 0.03589270, 0.04901011, 0.04089908, 0.06441121, 0.07927065, 0.11663929, //16200-16207
	0.03420964, 0.03891565, 0.04997103, 0.06484259, 0.05587856, 0.08015529, 0.09674542, 0.13437591, //16208-16215
	0.09910208, 0.10796450, 0.12541245, 0.15281635, 0.13711767, 0.17502473, 0.20765089, 0.26174013, //16216-16223
	0.06850943, 0.07639382, 0.09322540, 0.11792325, 0.10345891, 0.13996839, 0.16937572, 0.22353904, //16224-16231
	0.12655425, 0.13492835, 0.15036267, 0.17748514, 0.16217682, 0.19804721, 0.23129778, 0.28299066, //16232-16239
	0.28445488, 0.29565781, 0.31448106, 0.35249322, 0.33197875, 0.37682984, 0.42937308, 0.49038178, //16240-16247
	0.38362260, 0.39136787, 0.40368513, 0.43447672, 0.41825626, 0.45267332, 0.50003661, 0.54963562, //16248-16255
	0.00007766, 0.00019350, 0.00030434, 0.00064946, 0.00048620, 0.00079358, 0.00122355, 0.00149240, //16256-16263
	0.00021562, 0.00058208, 0.00089726, 0.00210930, 0.00155641, 0.00258292, 0.00441331, 0.00550436, //16264-16271
	0.00070380, 0.00167207, 0.00229437, 0.00498693, 0.00384057, 0.00586591, 0.00977435, 0.01189345, //16272-16279
	0.00181971, 0.00461433, 0.00613004, 0.01429055, 0.01105312, 0.01657389, 0.03058304, 0.03714654, //16280-16287
	0.00109867, 0.00279428, 0.00382049, 0.00895456, 0.00682831, 0.01054172, 0.01935864, 0.02383741, //16288-16295
	0.00349876, 0.00806835, 0.01008075, 0.02182055, 0.01749066, 0.02468272, 0.04391715, 0.05207114, //16296-16303
	0.00565131, 0.01391101, 0.01718943, 0.04049849, 0.03246029, 0.04545949, 0.09159145, 0.10804569, //16304-16311
	0.01285735, 0.02860206, 0.03334933, 0.07143677, 0.06000038, 0.07778636, 0.14693055, 0.16668406, //16312-16319
	0.00396457, 0.00492919, 0.00832035, 0.01216567, 0.00969043, 0.01792375, 0.02268372, 0.03846689, //16320-16327
	0.01368039, 0.01638936, 0.02399369, 0.03353428, 0.02763990, 0.04476199, 0.05643250, 0.08579554, //16328-16335
	0.02295712, 0.02627046, 0.03406013, 0.04515450, 0.03848243, 0.05659449, 0.06982756, 0.10003365, //16336-16343
	0.07544299, 0.08246717, 0.09630659, 0.11905799, 0.10601317, 0.13752746, 0.16641405, 0.21431445, //16344-16351
	0.05069852, 0.05677310, 0.06974080, 0.08975334, 0.07802557, 0.10762843, 0.13310640, 0.18005529, //16352-16359
	0.09765874, 0.10440666, 0.11680282, 0.13979971, 0.12685280, 0.15719355, 0.18696297, 0.23337312, //16360-16367
	0.24528088, 0.25527148, 0.27205481, 0.30726141, 0.28827943, 0.32981264, 0.38117138, 0.44076495, //16368-16375
	0.33798840, 0.34510993, 0.35643682, 0.38575104, 0.37029208, 0.40305350, 0.45040664, 0.50002851, //16376-16383
};
typedef int _compile_time_check_lperm16_ratios[sizeof(_lperm16_ratio) / sizeof(_lperm16_ratio[0]) == 16384 ? 1 : -1];
static int _lperm16_8(const Uint8 data[16]) {
	//single words
	int bit0 = (data[0] < data[1]), bit1 = (data[2] < data[3]), bit3 = (data[4] < data[5]), bit4 = (data[6] < data[7]);
	int bit7 = (data[8] < data[9]), bit8 = (data[10] < data[11]), bit10 = (data[12] < data[13]), bit11 = (data[14] < data[15]);
	//double words
	int bit2 = (data[1] < data[3]), bit5 = (data[5] < data[7]), bit9 = (data[9] < data[11]), bit12 = (data[13] < data[15]);
	//quadruple words & oct words
	int bit6 = (data[3] < data[7]), bit13 = (data[11] < data[15]), bit14 = (data[7] < data[15]);
	return bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3) | (bit4 << 4) | (bit5 << 5) | (bit6 << 6) | (bit7 << 7) | (bit8 << 8) | (bit9 << 9) | (bit10 << 10) | (bit11 << 11) | (bit12 << 12) | (bit13 << 13) | (bit14 << 14);
}
static int _lperm16_16(const Uint16 data[16]) {
	//single words
	int bit0 = (data[0] < data[1]), bit1 = (data[2] < data[3]), bit3 = (data[4] < data[5]), bit4 = (data[6] < data[7]);
	int bit7 = (data[8] < data[9]), bit8 = (data[10] < data[11]), bit10 = (data[12] < data[13]), bit11 = (data[14] < data[15]);
	//double words
	int bit2 = (data[1] < data[3]), bit5 = (data[5] < data[7]), bit9 = (data[9] < data[11]), bit12 = (data[13] < data[15]);
	//quadruple words & oct words
	int bit6 = (data[3] < data[7]), bit13 = (data[11] < data[15]), bit14 = (data[7] < data[15]);
	return bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3) | (bit4 << 4) | (bit5 << 5) | (bit6 << 6) | (bit7 << 7) | (bit8 << 8) | (bit9 << 9) | (bit10 << 10) | (bit11 << 11) | (bit12 << 12) | (bit13 << 13) | (bit14 << 14);
}
static int _lperm16_32(const Uint32 data[16]) {
	//single words
	int bit0 = (data[0] < data[1]), bit1 = (data[2] < data[3]), bit3 = (data[4] < data[5]), bit4 = (data[6] < data[7]);
	int bit7 = (data[8] < data[9]), bit8 = (data[10] < data[11]), bit10 = (data[12] < data[13]), bit11 = (data[14] < data[15]);
	//double words
	int bit2 = (data[1] < data[3]), bit5 = (data[5] < data[7]), bit9 = (data[9] < data[11]), bit12 = (data[13] < data[15]);
	//quadruple words & oct words
	int bit6 = (data[3] < data[7]), bit13 = (data[11] < data[15]), bit14 = (data[7] < data[15]);
	return bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3) | (bit4 << 4) | (bit5 << 5) | (bit6 << 6) | (bit7 << 7) | (bit8 << 8) | (bit9 << 9) | (bit10 << 10) | (bit11 << 11) | (bit12 << 12) | (bit13 << 13) | (bit14 << 14);
}
static int _lperm16_64(const Uint64 data[16]) {
	//single words
	int bit0 = (data[0] < data[1]), bit1 = (data[2] < data[3]), bit3 = (data[4] < data[5]), bit4 = (data[6] < data[7]);
	int bit7 = (data[8] < data[9]), bit8 = (data[10] < data[11]), bit10 = (data[12] < data[13]), bit11 = (data[14] < data[15]);
	//double words
	int bit2 = (data[1] < data[3]), bit5 = (data[5] < data[7]), bit9 = (data[9] < data[11]), bit12 = (data[13] < data[15]);
	//quadruple words & oct words
	int bit6 = (data[3] < data[7]), bit13 = (data[11] < data[15]), bit14 = (data[7] < data[15]);
	return bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3) | (bit4 << 4) | (bit5 << 5) | (bit6 << 6) | (bit7 << 7) | (bit8 << 8) | (bit9 << 9) | (bit10 << 10) | (bit11 << 11) | (bit12 << 12) | (bit13 << 13) | (bit14 << 14);
}
static int _lperm8_8(const Uint8 data[8]) {
	//single words
	int bit0 = (data[0] < data[1]), bit1 = (data[2] < data[3]), bit3 = (data[4] < data[5]), bit4 = (data[6] < data[7]);
	//double words
	int bit2 = (data[1] < data[3]), bit5 = (data[5] < data[7]);
	//quadruple words & oct words
	int bit6 = (data[3] < data[7]);
	return bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3) | (bit4 << 4) | (bit5 << 5) | (bit6 << 6);
}

void PractRand::Tests::LPerm16::init(PractRand::RNGs::vRNG *known_good) {
	lperm_counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next_pass = blocks_per_pass - 1;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::LPerm16::get_name() const {
	std::ostringstream buf;
	buf << "LPerm16(" << word_bits;
	if (passes_at_once > 1 || blocks_per_pass > 1) buf << ",";
	if (passes_at_once > 1) buf << passes_at_once;
	if (blocks_per_pass > 1) buf << "/" << blocks_per_pass;
	buf << ")";
	return buf.str();
}
void PractRand::Tests::LPerm16::get_results(std::vector<TestResult> &results) {
	int perms_per_block = !passes_at_once ? (8 * TestBlock::SIZE / word_bits / 16) : passes_at_once;
	if (blocks_tested * perms_per_block / (blocks_per_pass ? blocks_per_pass : 1) < LPERM_BUCKETS * 32) return;
	if (blocks_tested >= 1ull << 39) return;// believed to produce false positives when test length approachs 2**50 bytes
	const int fact4 = 24;
	const int fact8 = 40320;
	//double lperm4_chances[8] = { 2 / 24., 4 / 24., 2 / 24., 4 / 24., 4 / 24., 2 / 24., 4 / 24., 2 / 24. };
	double lperm4_chances[8] = { 3 / 24., 5 / 24., 1 / 24., 3 / 24., 3 / 24., 1 / 24., 5 / 24., 3 / 24. };
	int lperm4_greaterthan[8] = { 0, 0, 1, 1, 1, 2, 2, 3 };
	int lperm4_lessthan[8] = { 3, 2, 2, 1, 1, 1, 0, 0 };
	int lperm4_incomparable[8] = { 0, 1, 0, 1, 1, 0, 1, 0 };
	double lperm8_chances[128];
	int lperm8_greaterthan[128];
	int lperm8_lessthan[128];
	int lperm8_incomparable[128];
	for (int i = 0; i < 128; i++) lperm8_chances[i] = 0;
	for (int i = 0; i < fact8; i++) {
		Uint8 rawperm[8];
		Uint8 used = 0;
		if (i == 1736) {
			std::printf("");
		}
		int u = i;
		int d;
		for (int x = 0; x < 8; x++) {
			d = u % (8 - x);
			u /= 8 - x;
			int td = 0;
			int d2 = d;
			for (int y = 0; true; y++) {
				if ((used>>y) & 1) continue;
				if (!d2--) {
					td = y;
					break;
				}
			}
			used |= 1 << td;
			rawperm[x] = td;
		}
		if (used != 255) std::printf("duplicated values\n");
		int limited = _lperm8_8(rawperm);
		lperm8_chances[limited] += 1.0 / fact8;
	}
	for (int i = 0; i < 128; i++) {
		int low = i & 7;
		int high = (i >> 3) & 7;
		int which = i >> 6;
		double low_chance = lperm4_chances[low];
		double high_chance = lperm4_chances[high];
		double chance = low_chance * high_chance;
		double actual_chance = lperm8_chances[i] + lperm8_chances[i ^ 64];
		if (actual_chance - chance > 0.0000001) {
			issue_error("chances don't add up");
		}
		double observed_ratio = lperm8_chances[i] / chance;
		int low_odds = lperm4_greaterthan[low] - lperm4_lessthan[low];
		int high_odds = lperm4_greaterthan[high] - lperm4_lessthan[high];
		lperm8_greaterthan[i] = which ? lperm4_greaterthan[high] : lperm4_greaterthan[high] + lperm4_greaterthan[low] + 1;
		lperm8_lessthan[i] = which ? lperm4_lessthan[high] + lperm4_lessthan[low] + 1 : lperm4_lessthan[high];
		//if (!which) std::printf("%3d (%d:%d): %.5f %d\n", i, low, high, observed_ratio, high_odds - low_odds);
	}
	/*for (int i = 0; i < LPERM_BUCKETS; i++) {
		int low = i & 127;
		int high = (i >> 7) & 127;
		int meta = (((i >> 2) & 1) << 0) | (((i >> 5) & 1) << 1) | (((i >> 6) & 1) << 2) | (((i >> 9) & 1) << 3) | (((i >> 12) & 1) << 4) | (((i >> 13) & 1) << 5) | (((i >> 14) & 1) << 6);
		int which = i >> 14;
		lperm16_chances[i] = lperm8_chances[low] * lperm8_chances[high];
		double ratio = (lperm8_chances[meta] / (lperm8_chances[meta] + lperm8_chances[meta ^ 64]));
		lperm16_chances[i] *= ratio;// 
	}
	double sum = 0;
	for (int i = 0; i < 128; i++) sum += lperm8_chances[i];
	//std::printf("%f, ", sum);
	sum = 0;
	for (int i = 0; i < LPERM_BUCKETS; i++) sum += lperm16_chances[i];*/
	//std::printf("%f\n", sum);
	//double lowest = 1.0;
	//for (int i = 0; i < 128; i++) if (lowest > lperm8_chances[i]) lowest = lperm8_chances[i];
	//std::printf("%f\n", lowest);
	/*
		lperm8 chances are exact, but lperm16 chances are a crude approximation
	*/
	/*const Uint64 *counts_ = lperm_counts.get_array();
	double odds_sum[60] = { 0 };
	double odds_sum2[60] = { 0 };
	Uint64 odds_count[60] = { 0 };
	for (int i = 0; i < LPERM_BUCKETS / 2; i++) {
		int low = i & 127;
		int high = (i >> 7) & 127;
		int meta = (((i >> 2) & 1) << 0) | (((i >> 5) & 1) << 1) | (((i >> 6) & 1) << 2) | (((i >> 9) & 1) << 3) | (((i >> 12) & 1) << 4) | (((i >> 13) & 1) << 5) | (((i >> 14) & 1) << 6);
		int which = i >> 14;
		double base_chance = lperm8_chances[low] * lperm8_chances[high];
		double observed_ratio = counts_[i] / double(counts_[i] + counts_[i ^ 16384]);
		if (!double(counts_[i] + counts_[i ^ 16384])) continue;
		int low_odds = lperm8_greaterthan[low] - lperm8_lessthan[low];
		int high_odds = lperm8_greaterthan[high] - lperm8_lessthan[high];
		int odds = high_odds - low_odds;
		//if (odds < 0) { continue; odds = -odds; observed_ratio = 1 - observed_ratio; }
		//std::printf("%f\n", observed_ratio);
		if (std::fabs(observed_ratio - 0.5) > 0.2) {
			std::printf("");
		}
		odds += 30;
		if (odds < 0 || odds >= 60) continue;
		odds_count[odds]++;
		odds_sum[odds] += observed_ratio;
		odds_sum2[odds] += observed_ratio * observed_ratio;
		//if (!which) std::printf("%5d (%2x:%2x): %.6f %d\n", i, low, high, observed_ratio, high_odds - low_odds)
	}
	for (int i = 0; i < 60; i++) {
		if (odds_count[i]) {
			double mean = odds_sum[i] / odds_count[i];
			double meansqr = odds_sum2[i] / odds_count[i];
			std::printf("odds %+d: ratio %.8f +/- %.8f\n", i - 30, mean, std::sqrt(meansqr - mean * mean));
		}
	}
	for (int i = 0; i < LPERM_BUCKETS / 2; i++) {
		double observed_ratio = counts_[i] / double(counts_[i] + counts_[i ^ 16384]);
		std::printf("%.8f, ", observed_ratio);
		if ((i % 8) == 7) {
			std::printf("//%5d-%5d\n", i - 7, i);
		}
	}
	std::printf("\n");*/

	std::vector<double> lperm16_chances;
	lperm16_chances.resize(LPERM_BUCKETS);
	for (int i = 0; i < LPERM_BUCKETS / 2; i++) {
		int low = i & 127;
		int high = (i >> 7) & 127;
		double base_chance = lperm8_chances[low] * lperm8_chances[high];
		lperm16_chances[i] = base_chance * _lperm16_ratio[i];
		lperm16_chances[i + LPERM_BUCKETS/2] = base_chance * (1 - _lperm16_ratio[i]);
	}

	const Uint64 *counts_ = lperm_counts.get_array();
	double chisqr = g_test(LPERM_BUCKETS, &lperm16_chances[0], counts_);
	//double chisqr = g_test_flat(LPERM_BUCKETS, counts_);
	double n = math_chisquared_to_normal(chisqr, LPERM_BUCKETS - 1);
	//double n = g_test_flat_merge_normal(LPERM_BUCKETS, counts_);
	results.push_back(TestResult(get_name(), n, n, TestResult::TYPE_RAW_NORMAL, 0.05));
}
void PractRand::Tests::LPerm16::test_blocks(TestBlock *data, int numblocks) {
	blocks_tested += numblocks;
	while (numblocks > blocks_till_next_pass) {
		data += blocks_till_next_pass;
		numblocks -= blocks_till_next_pass;

		if (word_bits == 8) {
			//int perms_per_block = !passes_at_once ? (8 * TestBlock::SIZE / word_bits / 16) : passes_at_once;
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 7) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _lperm16_8(&data[0].as8[pos * 16]);
				lperm_counts.increment(code);
			}
		}
		if (word_bits == 16) {
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 8) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _lperm16_16(&data[0].as16[pos * 16]);
				lperm_counts.increment(code);
			}
		}
		if (word_bits == 32) {
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 9) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _lperm16_32(&data[0].as32[pos * 16]);
				lperm_counts.increment(code);
			}
		}
		if (word_bits == 64) {
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 10) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _lperm16_64(&data[0].as64[pos * 16]);
				lperm_counts.increment(code);
			}
		}

		data += 1;
		numblocks -= 1;
		blocks_till_next_pass = blocks_per_pass - 1;
	}
	blocks_till_next_pass -= numblocks;
}


template<typename T, typename CompareFunctor> int _generic_lperm16(const T data[16], const CompareFunctor &Compare) {
	//single words
	int bit0 = Compare(data[0], data[1]), bit1 = Compare(data[2], data[3]), bit3 = Compare(data[4], data[5]), bit4 = Compare(data[6], data[7]);
	int bit7 = Compare(data[8], data[9]), bit8 = Compare(data[10], data[11]), bit10 = Compare(data[12], data[13]), bit11 = Compare(data[14], data[15]);
	//double words
	int bit2 = Compare(data[1], data[3]), bit5 = Compare(data[5], data[7]), bit9 = Compare(data[9], data[11]), bit12 = Compare(data[13], data[15]);
	//quadruple words & oct words
	int bit6 = Compare(data[3], data[7]), bit13 = Compare(data[11], data[15]), bit14 = Compare(data[7], data[15]);
	return bit0 | (bit1 << 1) | (bit2 << 2) | (bit3 << 3) | (bit4 << 4) | (bit5 << 5) | (bit6 << 6) | (bit7 << 7) | (bit8 << 8) | (bit9 << 9) | (bit10 << 10) | (bit11 << 11) | (bit12 << 12) | (bit13 << 13) | (bit14 << 14);
}
template<typename T> class _laperm16_compare_functor {
public:
	bool operator() (T a, T b) const {
		int aw = count_ones_(a);
		int bw = count_ones_(b);
		if (aw == bw) return a < b;
		else return aw < bw;
	}
};

void PractRand::Tests::LPerm16A::init(PractRand::RNGs::vRNG *known_good) {
	lperm_counts.reset_counts();
	blocks_tested = 0;
	blocks_till_next_pass = blocks_per_pass - 1;
	TestBaseclass::init(known_good);
}
std::string PractRand::Tests::LPerm16A::get_name() const {
	std::ostringstream buf;
	buf << "LPerm16A(" << word_bits;
	if (passes_at_once > 1 || blocks_per_pass > 1) buf << ",";
	if (passes_at_once > 1) buf << passes_at_once;
	if (blocks_per_pass > 1) buf << "/" << blocks_per_pass;
	buf << ")";
	return buf.str();
}
void PractRand::Tests::LPerm16A::get_results(std::vector<TestResult> &results) {
	int perms_per_block = !passes_at_once ? (8 * TestBlock::SIZE / word_bits / 16) : passes_at_once;
	if (blocks_tested * perms_per_block / (blocks_per_pass ? blocks_per_pass : 1) < LPERM_BUCKETS * 32) return;
	const int fact4 = 24;
	const int fact8 = 40320;
	//double lperm4_chances[8] = { 2 / 24., 4 / 24., 2 / 24., 4 / 24., 4 / 24., 2 / 24., 4 / 24., 2 / 24. };
	double lperm4_chances[8] = { 3 / 24., 5 / 24., 1 / 24., 3 / 24., 3 / 24., 1 / 24., 5 / 24., 3 / 24. };
	int lperm4_greaterthan[8] = { 0, 0, 1, 1, 1, 2, 2, 3 };
	int lperm4_lessthan[8] = { 3, 2, 2, 1, 1, 1, 0, 0 };
	int lperm4_incomparable[8] = { 0, 1, 0, 1, 1, 0, 1, 0 };
	double lperm8_chances[128];
	int lperm8_greaterthan[128];
	int lperm8_lessthan[128];
	int lperm8_incomparable[128];
	for (int i = 0; i < 128; i++) lperm8_chances[i] = 0;
	for (int i = 0; i < fact8; i++) {
		Uint8 rawperm[8];
		Uint8 used = 0;
		if (i == 1736) {
			std::printf("");
		}
		int u = i;
		int d;
		for (int x = 0; x < 8; x++) {
			d = u % (8 - x);
			u /= 8 - x;
			int td = 0;
			int d2 = d;
			for (int y = 0; true; y++) {
				if ((used>>y) & 1) continue;
				if (!d2--) {
					td = y;
					break;
				}
			}
			used |= 1 << td;
			rawperm[x] = td;
		}
		if (used != 255) std::printf("duplicated values\n");
		int limited = _lperm8_8(rawperm);
		lperm8_chances[limited] += 1.0 / fact8;
	}
	for (int i = 0; i < 128; i++) {
		int low = i & 7;
		int high = (i >> 3) & 7;
		int which = i >> 6;
		double low_chance = lperm4_chances[low];
		double high_chance = lperm4_chances[high];
		double chance = low_chance * high_chance;
		double actual_chance = lperm8_chances[i] + lperm8_chances[i ^ 64];
		if (actual_chance - chance > 0.0000001) {
			issue_error("chances don't add up");
		}
		double observed_ratio = lperm8_chances[i] / chance;
		int low_odds = lperm4_greaterthan[low] - lperm4_lessthan[low];
		int high_odds = lperm4_greaterthan[high] - lperm4_lessthan[high];
		lperm8_greaterthan[i] = which ? lperm4_greaterthan[high] : lperm4_greaterthan[high] + lperm4_greaterthan[low] + 1;
		lperm8_lessthan[i] = which ? lperm4_lessthan[high] + lperm4_lessthan[low] + 1 : lperm4_lessthan[high];
	}

	std::vector<double> lperm16_chances;
	lperm16_chances.resize(LPERM_BUCKETS);
	for (int i = 0; i < LPERM_BUCKETS / 2; i++) {
		int low = i & 127;
		int high = (i >> 7) & 127;
		double base_chance = lperm8_chances[low] * lperm8_chances[high];
		lperm16_chances[i] = base_chance * _lperm16_ratio[i];
		lperm16_chances[i + LPERM_BUCKETS/2] = base_chance * (1 - _lperm16_ratio[i]);
	}

	const Uint64 *counts_ = lperm_counts.get_array();
	double chisqr = g_test(LPERM_BUCKETS, &lperm16_chances[0], counts_);
	//double chisqr = g_test_flat(LPERM_BUCKETS, counts_);
	double n = math_chisquared_to_normal(chisqr, LPERM_BUCKETS - 1);
	//double n = g_test_flat_merge_normal(LPERM_BUCKETS, counts_);
	results.push_back(TestResult(get_name(), n, n, TestResult::TYPE_RAW_NORMAL, 0.05));
}
void PractRand::Tests::LPerm16A::test_blocks(TestBlock *data, int numblocks) {
	blocks_tested += numblocks;
	while (numblocks > blocks_till_next_pass) {
		data += blocks_till_next_pass;
		numblocks -= blocks_till_next_pass;

		if (word_bits == 8) {
			//int perms_per_block = !passes_at_once ? (8 * TestBlock::SIZE / word_bits / 16) : passes_at_once;
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 7) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _generic_lperm16(&data[0].as8[pos * 16], _laperm16_compare_functor<Uint8>());
				lperm_counts.increment(code);
			}
		}
		if (word_bits == 16) {
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 8) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _generic_lperm16(&data[0].as16[pos * 16], _laperm16_compare_functor<Uint16>());
				lperm_counts.increment(code);
			}
		}
		if (word_bits == 32) {
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 9) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _generic_lperm16(&data[0].as32[pos * 16], _laperm16_compare_functor<Uint32>());
				lperm_counts.increment(code);
			}
		}
		if (word_bits == 64) {
			int perms_per_block = !passes_at_once ? ((8 * TestBlock::SIZE) >> 10) : passes_at_once;
			for (int pos = 0; pos < perms_per_block; pos++) {
				int code = _generic_lperm16(&data[0].as64[pos * 16], _laperm16_compare_functor<Uint64>());
				lperm_counts.increment(code);
			}
		}

		data += 1;
		numblocks -= 1;
		blocks_till_next_pass = blocks_per_pass - 1;
	}
	blocks_till_next_pass -= numblocks;
}










PractRand::Tests::Transforms::multiplex::multiplex(const char *name_, const ListOfTests &testlist)
:
	subtests(testlist)
{
//	for (unsigned int i = 0; i < testlist.tests.size(); i++) subtests.push_back(testlist.tests[i]);
	if (name_) name = name_;
	//else if (subtests.tests.size() == 1) name = subtests.tests[0]->get_name();
	else {
		//std::ostringstream str;
		//str << "{" << int(subtests.size()) << "}";
		//name = str.str();
		name = "";
	}
//	else name = make_string("");
}
bool PractRand::Tests::Transforms::multiplex::recommend_subtest_tree_descent() const {
	return true;
}
void PractRand::Tests::Transforms::multiplex::deinit() {
	for (std::vector<Tests::TestBaseclass*>::iterator it = subtests.tests.begin(); it != subtests.tests.end(); it++)
		(*it)->deinit();
}
void PractRand::Tests::Transforms::multiplex::init( RNGs::vRNG *known_good ) {
	blocks_already = 0;
	for (std::vector<Tests::TestBaseclass*>::iterator it = subtests.tests.begin(); it != subtests.tests.end(); it++)
		(*it)->init(known_good);
}
PractRand::Tests::Transforms::multiplex::~multiplex ( ) {
	for (std::vector<Tests::TestBaseclass*>::iterator it = subtests.tests.begin(); it != subtests.tests.end(); it++)
		delete (*it);
	subtests.tests.clear();
}
std::string PractRand::Tests::Transforms::multiplex::get_name() const {
	return name.c_str();
}
void PractRand::Tests::Transforms::multiplex::test_blocks(TestBlock *data, int numblocks) {
	for (std::vector<Tests::TestBaseclass*>::iterator it = subtests.tests.begin(); it != subtests.tests.end(); it++)
		(*it)->test_blocks(data, numblocks);
	blocks_already += numblocks;
}
static std::pair<int,std::pair<int,int> > extract_low_transform_params(const std::string name) {
	std::pair<int,std::pair<int,int> > fail(0, std::pair<int,int>(0,0));
	int first, last;
	char termination;
	const char *c = name.c_str();
	int r = std::sscanf(c, "[Low%d/%d%c", &first, &last, &termination);
	if (r != 3 || termination != ']') return fail;
	return std::pair<int,std::pair<int,int> >(int(strchr(c, ']') - c + 1), std::pair<int,int>(first,last));
}
static std::string combine_transform_names(const std::string &prefix, const std::string &name) {
	std::string fail = prefix + name;
	std::pair<int,std::pair<int,int> > a = extract_low_transform_params(prefix);
	std::pair<int,std::pair<int,int> > b = extract_low_transform_params(name);
	if (a.second.first != b.second.second) return fail;
	if (!a.second.second || !b.second.first) return fail;
	if (a.first != prefix.length()) return fail;
	if (b.first > name.length()) return fail;
	std::ostringstream buf;
	buf << "[Low" << b.second.first << "/" << a.second.second << "]" << name.substr(b.first);
	return buf.str();
}
void PractRand::Tests::Transforms::multiplex::get_results(std::vector<TestResult> &results) {
	size_t old_size = results.size();
	for (std::vector<Tests::TestBaseclass*>::iterator it = subtests.tests.begin(); it != subtests.tests.end(); it++) {
		(*it)->get_results(results);
	}
	for (size_t i = old_size; i < results.size(); i++) {
		results[i].name = combine_transform_names(get_name(), results[i].name);
		//results[i].name = get_name() + results[i].name;
	}
}
int PractRand::Tests::Transforms::multiplex::get_blocks_to_repeat() const {
	int rv = 0;
	for (std::vector<Tests::TestBaseclass*>::const_iterator it = subtests.tests.begin(); it != subtests.tests.end(); it++) {
		int x = (*it)->get_blocks_to_repeat();
		if (rv < x) rv = x;
	}
	return rv;
}
int PractRand::Tests::Transforms::multiplex::get_num_children() const { return subtests.tests.size(); }
Tests::TestBaseclass *PractRand::Tests::Transforms::multiplex::get_child  (int index) const {return subtests.tests[index];}
//std::string PractRand::Tests::Transforms::multiplex::get_child_name  (int index) const {return subtests[index]->get_name();}
//double      PractRand::Tests::Transforms::multiplex::get_child_result(int index) {return subtests[index]->get_result();}

Uint64 PractRand::Tests::Transforms::multiplex::get_blocks_passed_through(int index) const {
	return blocks_already;
}
Uint64 PractRand::Tests::Transforms::switching::get_blocks_passed_through(int index) const {
	return blocks_already_per[index];
}


PractRand::Tests::Transforms::switching::switching(
	const char *name_, 
	const ListOfTests &testlist, 
	std::vector<Uint64> lengths_)
:
	multiplex(name_, testlist),
	lengths(lengths_),
	phase(0),
	which(0)
{
	if (lengths.size() != testlist.tests.size()) issue_error();
	blocks_already_per.resize(lengths.size());
	for (unsigned long i = 0; i < blocks_already_per.size(); i++)
		blocks_already_per[i] = 0;
	total_length = 0;
	for (unsigned long i = 0; i < blocks_already_per.size(); i++) total_length += lengths[i];
}
PractRand::Tests::Transforms::switching::switching(
	const char *name_, 
	const ListOfTests &testlist, 
	Uint64 length)
:
	multiplex(name_, testlist)
{
	lengths.resize(testlist.tests.size());
	blocks_already_per.resize(lengths.size());
	for (unsigned long i = 0; i < blocks_already_per.size(); i++)
		lengths[i] = length;
	total_length = length * blocks_already_per.size();
}
void PractRand::Tests::Transforms::switching::init( RNGs::vRNG *known_good ) {
	for (unsigned long i = 0; i < blocks_already_per.size(); i++)
		blocks_already_per[i] = 0;
	phase = 0;
	which = 0;
	multiplex::init(known_good);
}
void PractRand::Tests::Transforms::switching::test_blocks( TestBlock *data, int numblocks_ ) {
	Uint64 numblocks = numblocks_;
	if (phase + numblocks < lengths[which]) {
		phase += numblocks;
		subtests.tests[which]->test_blocks(data, numblocks_);
		return;
	}
	Uint64 part = lengths[which] - phase;
	subtests.tests[which]->test_blocks(data, Uint32(part));
	phase = 0;
	if (++which == subtests.tests.size()) which = 0;
	data += part;
	numblocks_ -= part;
	test_blocks(data, numblocks_);
}
//double PractRand::Tests::Transforms::switching::get_result() {
//}
	
void PractRand::Tests::Transforms::Transform_Baseclass::init( RNGs::vRNG *known_good ) {
	Transforms::multiplex::init(known_good);
	leftovers = 0;
//	buffered.reserve( flush_size * TESTBLOCK_SIZE );
	buffered.clear();
}
void PractRand::Tests::Transforms::Transform_Baseclass::flush(bool aggressive) {
	int blocks_to_repeat = get_blocks_to_repeat();
	int minblocks = aggressive ? 1 : flush_size;
	int numblocks = buffered.size() - (leftovers?1:0);
	int old = (blocks_already > blocks_to_repeat) ? blocks_to_repeat : (int)blocks_already;
	numblocks -= old;
	if (numblocks < minblocks) return;
	Transforms::multiplex::test_blocks((TestBlock*)&buffered[old], numblocks);
//	int newold = blocks_already;
//	if (blocks_already > TestBaseclass::REPEATED_BLOCKS) newold = TestBaseclass::REPEATED_BLOCKS;
	int newold = (blocks_already > blocks_to_repeat) ? blocks_to_repeat : (int)blocks_already;
	int blocks_to_move = newold + (leftovers?1:0);
	int how_far = buffered.size() - blocks_to_move;
	if (blocks_to_move && how_far) std::memmove(&buffered[0], &buffered[how_far], blocks_to_move * TestBlock::SIZE);
	buffered.resize(blocks_to_move);
}
