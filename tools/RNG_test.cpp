
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
//#include <string>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <sstream>
#include <memory>

#define MULTITHREADING_SUPPORTED
#define CPP_2011_SUPPORTED

//master header, includes everything in PractRand for both 
//  practical usage and research... 
//  EXCEPT it does not include specific algorithms
#include "PractRand_full.h"

//specific algorithms: all recommended RNGs
#include "PractRand/RNGs/all.h"

//specific algorithms: non-recommended RNGs
#include "PractRand/RNGs/other/mt19937.h"
#include "PractRand/RNGs/other/transform.h"
#include "PractRand/RNGs/other/lcgish.h"
#include "PractRand/RNGs/other/oddball.h"
#include "PractRand/RNGs/other/simple.h"
#include "PractRand/RNGs/other/cbuf.h"
#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/other/special.h"

#ifdef WIN32 // needed to allow binary stdin on windows
#include <fcntl.h>
#include <io.h>
#endif

using namespace PractRand;

PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);

//helpers for the test programs, to deal with RNG names, test usage, etc
#include "RNG_from_name.h"

#include "TestManager.h"
#ifdef MULTITHREADING_SUPPORTED
#include "multithreading.h"
#include "MultithreadedTestManager.h"
#endif
#include "Candidate_RNGs.h"

#if defined CPP_2011_SUPPORTED
#include <chrono>
typedef std::chrono::system_clock::rep TimeUnit;
TimeUnit get_time() { return std::chrono::system_clock::now().time_since_epoch().count(); }
double get_time_period() { return std::chrono::system_clock::period().num / (double)std::chrono::system_clock::period().den; }
#else
typedef std::clock_t TimeUnit;
std::clock_t get_time() {return std::clock();}
double get_time_period() { return 1.0 / CLOCKS_PER_SEC; }
#endif

#include "dummy_rng.h"

void print_128bit_hex(Uint64 low, Uint64 high) {
	//two complications:
	//1. printing things as longs in 32 bit chunks - I don't want to assume llx support, and while I'm willing to assume long is at least 32 bits I don't want to assume anything more specific
	//2. I don't want to print more digits than necessary
	if (high >> 32) std::printf("0x%lx%08lx%08lx%08lx", long(high >> 32), long(Uint32(high)), long(low >> 32), long(Uint32(low)));
	else if (high) std::printf("0x%lx%08lx%08lx", long(high), long(low >> 32), long(Uint32(low)));
	else if (low >> 32) std::printf("0x%lx%08lx", long(low >> 32), long(Uint32(low)));
	else std::printf("0x%lx", long(low));
}

double print_result(const PractRand::TestResult &result, bool print_header = false, int table_format = 0) {

	double raw_suspicion = result.get_suspicion();
	double suspicion_level = std::fabs(raw_suspicion) + 1 - 1;// +1 for suspicion conversion, -1 to account for there being 2 failure regions (near-zero and near-1)
	double weight = result.get_weight();//already normalized by caller
	double weight_adjustment = -std::log2(weight);
	//double as = (std::fabs(result.get_suspicion()) + 1 - 1) * dec;// +1 for suspicion conversion, -1 to account for there being 2 failure regions (near-zero and near-1)
	//double wmod = std::log(result.get_weight()) / std::log(0.5) * dec;
	//double rs = as - wmod;
	//double ap = std::fabs(0.5 - result.get_pvalue());
	double decimal_conversion_factor = std::log(2.) / std::log(10.0);
	double score = (suspicion_level - weight_adjustment) * decimal_conversion_factor;
	if (score < 0) score = 0;

	if (print_header) {
		if (table_format == 2) {//alternative format, format 2
			//looking at alternative table format ideas:
			std::printf("  Test Name                         Raw   p-value     W%%  Scor Evaluation       \n");
			//           00        10        20        30        40        50   %     60        70        80
			/*
				"  "               2 characters for initial indentation
				"DC6-9x1bytes-1"   34 characters for test name
				"+0.3  "           6 characters for "raw" result
				"1-2.3e-456 "      11 characters for p-value or alternative
				"135 "             5 characters for weight (normalized)
				"9.7 "             5 characters for "score" (normalized decimal suspicion)
				"mildly suspcious" 17 characters for evalation
				2+34+6+11+5+5+17 = 80 characters total
				so if there's a line wrap at 80 characters, then it doesn't quite fit
				I can replace Weight% with WMod to save 1 more character?
			*/
		}
		else if (table_format == 3) {//scoring format, format 3
			std::printf("  Test Name                         Raw   Sus  WMod Scor Evaluation       \n");
			//           00        10        20        30        40        50        60        70        80
			/*
				"  "               2 characters for initial indentation
				"DC6-9x1bytes-1"   34 characters for test name
				"+0.3  "           6 characters for "raw" result
				"+0.12 "           5 characters for (decimal) suspicion (same information as p-value)
				"1.3  "            5 characters for weight modifier (decimal)
				"9.7  "            5 characters for score (normalized decimal suspicion)
				"mildly suspcious" 17 characters for evalation
				2+34+6+5+5+5+17 = 74 characters total
				so if there's a line wrap at 80 characters, then it easily fits
			*/
		}
		else if (table_format == 0) {//original format, format 0
			std::printf("  Test Name                         Raw       Processed     Evaluation       \n");
			//           00        10        20        30        40        50        60        70        80
			/*
				2 characters for initial indentation
				34 characters for test name
				10 characters for "raw" result
				14 characters for p-value or alternative
				17 characters for evalation
				2+34+10+14+17 = 77 characters total
				so if there's a line wrap at 80 characters, there's still room for 2 more
			*/
		}
		else if (table_format == 1) {//wide format, format 1
			std::printf("  Test Name                         Raw       Processed     Weight%%  Score  Evaluation       \n");
			//           00        10        20        30        40        50        60    %    70        80        90        100
			/*
				"  "               2 characters for initial indentation
				"DC6-9x1bytes-1"   34 characters for test name
				"R = +0.3  "       10 characters for "raw" result
				"p =1-2.3e-456 "   14 characters for p-value or alternative
				"W=0.135  "        9 characters for weight (normalized)
				"S=9.7  "          7 characters for "score" (normalized decimal suspicion)
				"mildly suspcious" 17 characters for evalation
				2+34+10+14+9+7+17 = 93 characters total
				so if there's a line wrap at 100 characters, there's still room for 4 more
			*/
		}
		else {
			issue_error("unrecognized table format");
		}
	}
	std::printf("  ");//2 characters, used by all formats
	//NAME
	if (true) {// 34 characters, used by all formats
		std::printf("%s", result.name.c_str());
		int len = result.name.length();
		for (int i = len; i < 34; i++) std::printf(" ");
	}

	//RAW TEST RESULT
	if (table_format == 0 || table_format == 1) {// 10 characters?
		double raw = result.get_raw();
		if (std::fabs(raw) < 999.95) std::printf("R=%+6.1f  ", raw);
		else if (std::fabs(raw) <= 99999) std::printf("R=%+6.0f  ", raw);
		else std::printf("R%s99999  ", raw < 0 ? "<-" : ">+");
	}
	else if (table_format == 2 || table_format == 3) {// 6 characters
		double raw = result.get_raw();
		if (std::fabs(raw) < 99.95) std::printf("%+5.1f ", raw);
		else if (std::fabs(raw) < 9999) std::printf("%+5.0f ", raw);
		else std::printf("%s9999 ", raw < 0 ? "-" : "+");
	}

	//RESULT (as a p-value if reasonable, otherwise as "pass" or "fail")
	if (table_format == 0 || table_format == 1) {// 14 characters
		if (result.type == result.TYPE_PASSFAIL)
			std::printf("  %s      ", result.get_pvalue() ? "\"pass\"" : "\"fail\"");
		else if (result.type == result.TYPE_RAW || result.type == result.TYPE_UNKNOWN)
			std::printf("              ");
		else if (result.type == result.TYPE_BAD_P || result.type == result.TYPE_GOOD_P || result.type == result.TYPE_RAW_NORMAL) {
			double p = result.get_pvalue();
			double a = std::fabs(p-0.5);
			std::printf((result.type != result.TYPE_GOOD_P) ? "p~= " : "p = ");
			if (a > 0.49) {
				double s = result.get_suspicion();
				double ns = std::fabs(s) + 1;
				double dec = ns / (std::log(10.0) / std::log(2.0));
				double dig = std::ceil(dec);
				double sig = std::floor(std::pow(0.1, dec - dig));
				if (dig > 999) std::printf(" %d        ", (s > 0) ? 1 : 0);
				else {
					if (s > 0) std::printf("1-%1.0fe-%.0f  ", sig, dig);
					else       std::printf("  %1.0fe-%.0f  ", sig, dig);
					if (dig < 100) std::printf(" ");
					if (dig < 10) std::printf(" ");
				}
			}
			else if (result.type == result.TYPE_GOOD_P) std::printf("%5.3f     ", p);
			else if (a >= 0.4)                          std::printf("%4.2f      ", p);
			else                                        std::printf("%3.1f       ", p);
		}
		else if (result.type == result.TYPE_BAD_S || result.type == result.TYPE_GOOD_S) {
			double s = result.get_suspicion();
			double p = result.get_pvalue();
			std::printf((result.type == result.TYPE_BAD_S || result.type == result.TYPE_RAW_NORMAL) ? "p~=" : "p =");
			if (p >= 0.01 && p <= 0.99) std::printf(" %.3f     ", p);
			else {
				double ns = std::fabs(s) + 1;
				double dec = ns / (std::log(10.0) / std::log(2.0));
				double dig = std::ceil(dec);
				double sig = std::pow(0.1, dec - dig);
				sig = std::floor(sig * 10) * 0.1;
				if (dig > 9999) std::printf(" %d         ", (s > 0) ? 1 : 0);
				else if (dig > 999) {
					sig = std::floor(sig);
					if (s > 0) std::printf("1-%1.0fe-%.0f  ", sig, dig);
					else       std::printf("  %1.0fe-%.0f  ", sig, dig);
				}
				else {
					if (s > 0) std::printf("1-%3.1fe-%.0f ", sig, dig);
					else       std::printf("  %3.1fe-%.0f ", sig, dig);
					if (dig < 100) std::printf(" ");
					if (dig < 10) std::printf(" ");
				}
			}
		}
	}
	else if (table_format == 2) {// 11 characters
		if (result.type == result.TYPE_PASSFAIL)
			std::printf("  %s   ", result.get_pvalue() ? "\"pass\"" : "\"fail\"");
		else if (result.type == result.TYPE_RAW || result.type == result.TYPE_UNKNOWN)
			std::printf("           ");
		else {
			double s = result.get_suspicion();
			double p = result.get_pvalue();
			double a = std::fabs(p - 0.5);
			// TODO
			std::printf("           ");
			if (result.type == result.TYPE_RAW_NORMAL || result.type == result.TYPE_BAD_P) {
			}
			else if (result.type == result.TYPE_GOOD_P || result.type == result.TYPE_BAD_S || result.type == result.TYPE_GOOD_S) {
			}
		}
	}

	//Sus (suspicioun, converted to decimal)
	if (table_format == 3) {// 5 characters
		double dsus = -raw_suspicion * decimal_conversion_factor;
		if (false) ;
		else if (std::fabs(dsus) < 9.95) std::printf("%+4.1f ", dsus);
		//else if (std::fabs(dsus) < 99.95) std::printf("%+5.1f ", dsus);
		else if (std::fabs(dsus) < 999.5) std::printf("%+4.0f ", dsus);
		else if (dsus < 0) std::printf("-999 ");
		else if (dsus > 0) std::printf("+999 ");
		else std::printf("nan? ");
	}

	//WEIGHT (normalized)
	if (table_format == 1) {// 9 characters
		std::printf("W=%6.3f ", weight * 100);
	}
	else if (table_format == 2) {// 5 characters
		//	"**   "
		if (weight >= 0.995) std::printf("**   ");
		//	"12   "
		else if (weight >= 0.1) std::printf("%2.0f   ", weight * 100);
		//	" 2.3 "
		else if (weight >= 0.0095) std::printf(" %3.1f ", weight * 100);
		//	" .34 "
		//	" .04 "
		else if (weight >= 0.000095) std::printf(" .%02d ", int(weight * 10000 + 0.5));
		//	"5e-3 "
		else if (weight >= 0.0000000001) {
			// might mangle an edge case once in a rare while?
			double dl = -std::log10(weight);
			double idl;
			double dlf = std::modf(dl, &idl);
			if (dlf < 0.022276394711152234) {// -log10(0.95)
				idl -= 1;
			}
			double dig = weight * std::pow(10.0, idl + 1);
			int idig = dig + 0.5;
			if (idig > 9) idig = 9;
			else if (idig < 1) idig = 1;
			int e = idl - 1;
			std::printf("%de-%d ", idig, e);
		}
		//  " 0   "
		else std::printf(" 0   ");//shouldn't be necessary, right?
	}

	//WMod (Weight-based score modifier)
	if (table_format == 3) {// 5 characters
		if (weight > 0.0000000001) {
			double wmod = -std::log10(weight);
			std::printf("%4.2f ", wmod);
		}
		else std::printf("???  "); // ...weight should NOT be that low, right?
	}

	//SCORE (normalized decimal suspicion)
	if (table_format == 1) {// 7 characters
		//	"S=0    "
		if (score == 0) {
			std::printf("S=0    ");
		}
		//	"S=0.0  "
		else if (score < 0.95) {
			std::printf("S=%.1f  ", score);
		}
		//	"S=1.0  "
		else if (score < 9.995) {
			std::printf("S=%.2f ", score);
		}
		//	"S=21.0 "
		else if (score < 99.95) {
			std::printf("S=%.1f  ", score);
		}
		//	"S=321  "
		else if (score < 999.5) {
			std::printf("S=%.0f  ", score);
		}
		//	"S=4321 "
		else if (score < 9999) {
			std::printf("S=%.0f ", score);
		}
		//	"S>9999 "
		else {
			std::printf("S>9999 ");
		}
	}
	else if (table_format == 2 || table_format == 3) {//5 characters
		if (score == 0) {
			std::printf("0    ");
		}
		else if (score < 9.95) {
			std::printf("%.1f  ", score);
		}
		else if (score < 99.5) {
			std::printf("%.0f.  ", score);
		}
		else if (score < 999.5) {
			std::printf("%.0f  ", score);
		}
		else if (score < 9999.5) {
			std::printf("%.0f ", score);
		}
		else {
			std::printf("9999+");//doesn't even need a seperater from the next column, consider that's guaranteed to be indented under the circumstances
		}
	}

	//EVALUATION - MESSAGE DESCRIBING RESULT IN ENGLISH
	if (true) {// 17 characters
		/*
			Threshold Values:
			The idea is to assign a suspicioun level based not just upon the 
			p-value but also the number of p-values and their relative importance.  
			If there are a million p-values then we probably don't care about 
			anything less extreme than a one in ten million event.  
			But if there's one important p-value and a million unimportant ones then 
			the important one doesn't have to be that extreme to rouse our suspicion.  

			Output Format:
			unambiguous failures are indented 2 spaces to make them easier to spot
			probable failures with (barely) enough room for ambiguity are indented 1 space
			the most extreme failures get a sequence of exclamation marks to distinguish them
		*/
		if (false) ;
		else if (score >999) std::printf("  FAIL !!!!!!!!  ");
		else if (score >325) std::printf("  FAIL !!!!!!!   ");
		else if (score >165) std::printf("  FAIL !!!!!!    ");
		else if (score > 85) std::printf("  FAIL !!!!!     ");
		else if (score > 45) std::printf("  FAIL !!!!      ");
		else if (score > 25) std::printf("  FAIL !!!       ");
		else if (score > 17) std::printf("  FAIL !!        ");
		else if (score > 12) std::printf("  FAIL !         ");
		else if (score >8.5) std::printf("  FAIL           ");
		else if (score >6.0) std::printf(" VERY SUSPICIOUS ");
		else if (score >4.0) std::printf("very suspicious  ");
		else if (score >3.0) std::printf("suspicious       ");
		else if (score >2.0) std::printf("mildly suspicious");
		else if (score >1.0) std::printf("unusual          ");
		else if (score >0.0) std::printf("normalish        ");
		else                 std::printf("normal           ");
		/*
			what if I shortened this to make room for more columns?
			in particular the score (normalized decimal suspicion) column is a partial replacement for this
			this would save 6 characters, from 17 down to 11
		if (false) ;
		else if (score <=  0) std::printf("normal     ");
		else if (score <   2) std::printf("normalish  ");//   0-2
		else if (score <   4) std::printf("unusual    ");//   2-4
		else if (score <   6) std::printf("suspicious ");//   4-6
		else if (score <   9) std::printf(" SUSPICIOUS");//   6-9
		else if (score <  25) std::printf("  FAIL     ");//   9-25
		else if (score < 125) std::printf("  FAIL !   ");//  25-125
		else if (score < 625) std::printf("  FAIL !!  ");// 125-625
		else if (score <3125) std::printf("  FAIL !!! ");// 625-3125
		else                  std::printf("  FAIL!!!!!");//   3125+
		*/
	}
	std::printf("\n");
	return score;
}

#include "SeedingTester.h"

const char *seed_str = NULL;
Uint64 seed_low, seed_high;

void show_checkpoint(TestManager &tman, int mode, double time, bool smart_thresholds, double threshold, bool end_on_failure, bool do_metatest, int table_format) {
	double metatest_result = 0;
	std::printf("rng=%s", tman.get_rng()->get_name().c_str());

	std::printf(", seed=");
	if (tman.get_rng()->get_flags() & PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED) {
		if (seed_str) std::printf("%s", seed_str);
		else std::printf("unknown");
	}
	else {
		print_128bit_hex(seed_low, seed_high);
	}
	std::printf("\n");

	std::printf("length= ");
	Uint64 length = tman.get_blocks_so_far() * Tests::TestBlock::SIZE;
	double log2b = std::log(double(length)) / std::log(2.0);
	const char *unitstr[6] = {"kilobyte", "megabyte", "gigabyte", "terabyte", "petabyte", "exabyte"};
	int units = int(std::floor(log2b / 10)) - 1;
	if (units < 0 || units > 5) {std::printf("internal error: length out of bounds?\n");std::exit(1);}
	if (length & (length-1))
		std::printf("%.3f %ss", length * std::pow(0.5,units*10.0+10), unitstr[units] );
	else std::printf("%.0f %s%s", length * std::pow(0.5,units*10.0+10), unitstr[units], length != (Uint64(1024)<<(units*10)) ? "s" : "" );
	if (length & (length-1)) std::printf(" (2^%.3f", log2b - (mode?3:0)); else std::printf(" (2^%.0f", log2b - (mode?3:0));
	const char *mode_unit_names[3] = {"bytes", "seeds", "entropy strings"};
	std::printf(" %s), time= ", mode_unit_names[mode]);
	if (time < 99.95) std::printf("%.1f seconds\n", time);
	else std::printf("%.0f seconds\n", time);

	std::vector<PractRand::TestResult> results;
	tman.get_results(results);
	for (int i = 0; i < results.size(); i++) {//sanity check results
		double weight = results[i].get_weight();
		if (std::isnan(weight) || std::isinf(weight) || weight <= 0) {
			std::printf("\nerror: garbage weight in test results\n");
			//std::fprintf(stderr, "error: garbage weight in test results\n");
			std::printf("\ttest=%s\n", results[i].name.c_str());
			//std::fprintf(stderr, "\ttest=%s\n", results[i].name.c_str());
			std::printf("\tweight=%f\n", weight);
			//std::fprintf(stderr, "\tweight=%f\n", weight);
			std::exit(1);
		}
		double p = results[i].get_pvalue();
		if (std::isnan(p) || std::isinf(p) || p < 0 || p > 1) {
			std::printf("\nerror: garbage p-value in test results\n");
			//std::fprintf(stderr, "error: garbage p-value in test results\n");
			std::printf("\ttest=%s\n", results[i].name.c_str());
			//std::fprintf(stderr, "\ttest=%s\n", results[i].name.c_str());
			std::printf("\tpvalue=%f\n", p);
			//std::fprintf(stderr, "\tpvalue=%f\n", p);
			std::exit(1);
		}
	}
	double total_weight = 0;
	for (int i = 0; i < results.size(); i++) {
		double weight = results[i].get_weight();
		total_weight += weight;
	}
	std::vector<int> marked;
	for (int i = 0; i < results.size(); i++) {
		results[i].set_weight(results[i].get_weight() / total_weight);
		if (!smart_thresholds) {
			if (std::fabs(0.5 - results[i].get_pvalue()) < 0.5 - threshold) continue;
		}
		else {
			double T = threshold * results[i].get_weight() * 0.5;
			if (std::fabs(0.5 - results[i].get_pvalue()) < (0.5 - T)) continue;
		}
		marked.push_back(i);
	}
	double biggest_decimal_suspicion = 0;
	for (int i = 0; i < marked.size(); i++) {
		double decimal_suspicion = print_result(results[marked[i]], i == 0, table_format);
		if (decimal_suspicion > biggest_decimal_suspicion) biggest_decimal_suspicion = decimal_suspicion;
		if (decimal_suspicion > 0) metatest_result += decimal_suspicion * decimal_suspicion;
	}
	if (marked.size() == results.size())
		;
	else if (marked.size() == 0)
		std::printf("  no anomalies in %d test result(s)\n", int(results.size()));
	else
		std::printf("  ...and %d test result(s) without anomalies\n", int(results.size() - marked.size()));
	if (results.size() != 0 && do_metatest) {
		metatest_result = std::sqrt(metatest_result);
		std::string metatest_eval = "";
		if (metatest_result < 0.25) metatest_eval = "no problems";
		else if (metatest_result < 2) metatest_eval = "minor blip";
		else if (metatest_result < 4) metatest_eval = "major blip";
		else if (metatest_result < 7) metatest_eval = "suspicious";
		else if (metatest_result < 10) metatest_eval = "very suspicious";
		else metatest_eval = "extreme ; probably a failure";
		std::printf("  metatest score: %.1f, evaluated as: %s\n", metatest_result, metatest_eval.c_str());
	}
	std::printf("\n");
	std::fflush(stdout);
	if (end_on_failure && biggest_decimal_suspicion > 8.5) {
		std::exit(0);
	}
}
double interpret_length(const std::string &lengthstr, bool normal_mode) {
	//(0-9)*[.(0-9)*][((K|M|G|T|P)[B])|(s|m|h|d)]
	int mode_factor = normal_mode ? 1 : 8;
	int pos = 0;
	double value = 0;
	for (; pos < lengthstr.size(); pos++) {
		char c = lengthstr[pos];
		if (c < '0') break;
		if (c > '9') break;
		value = value * 10 + (c - '0');
	}
	if (!pos) return 0;
	if (pos == lengthstr.size()) return std::pow(2.0,value) * mode_factor;
	if (lengthstr[pos] == '.') {
		pos++;
		double sig = 0.1;
		for (; pos < lengthstr.size(); pos++,sig*=0.1) {
			char c = lengthstr[pos];
			if (c < '0') break;
			if (c > '9') break;
			value += (c - '0') * sig;
		}
		if (pos == lengthstr.size()) return std::pow(2.0,value) * mode_factor;
	}
	double scale = 0;
	bool expect_B = true;
	char c = lengthstr[pos];
	switch (c) {
	case 'K':
		scale = 1024.0;
		break;
	case 'M':
		scale = 1024.0 * 1024.0;
		break;
	case 'G':
		scale = 1024.0 * 1024.0 * 1024.0;
		break;
	case 'T':
		scale = 1024.0 * 1024.0 * 1024.0 * 1024.0;
		break;
	case 'P':
		scale = 1024.0 * 1024.0 * 1024.0 * 1024.0 * 1024.0;
		break;
	case 's':
		scale = -1;//one second
		expect_B = false;
		break;
	case 'm':
		scale = -60;//one minute
		expect_B = false;
		break;
	case 'h':
		scale = -3600;//one hour
		expect_B = false;
		break;
	case 'd':
		scale = -86400;//one day
		expect_B = false;
		break;
	}
	pos++;
	if (pos == lengthstr.size()) {
		if (scale < 0) {
			if (value < 0.05) value = 0.05;
			return value * scale;
		}
		return value * scale * mode_factor;
	}
	if (!expect_B) return 0;
	if (lengthstr[pos++] != 'B') return 0;
	if (pos != lengthstr.size()) return 0;
	return value * scale;
}
bool interpret_seed(const std::string &seedstr, Uint64 &seed_low, Uint64 &seed_high) {
	//would prefer strtol, but that is insufficiently portable when it has to handle 64 bit values
	//and now it's 128 bit values, so even less so
	Uint64 value_low = 0, value_high = 0;
	int position = 0;
	int used_characters = 0;
	if (seedstr.length() >= 3 && seedstr[0] == '0' && seedstr[1] == 'x') position = 2;
	while (position < seedstr.length()) {
		int c = seedstr[position++];
		if (value_high >> 60) return false;//too long
		value_high *= 16; value_high |= value_low >> 60;
		value_low *= 16;
		if (c >= '0' && c <= '9') value_low += (c - '0');
		else if (c >= 'a' && c <= 'f') value_low += (c - 'a') + 10;
		else if (c >= 'A' && c <= 'F') value_low += (c - 'A') + 10;
		else return false;//invalid character
		used_characters++;
	}
	seed_low = value_low;
	seed_high = value_high;
	return used_characters >= 1;
}

#include "PractRand/Tests/Birthday.h"
#include "PractRand/Tests/FPMulti.h"
#include "PractRand/Tests/FPF.h"
#include "PractRand/Tests/DistFreq4.h"
#include "PractRand/Tests/Gap16.h"
#include "PractRand/Tests/BRank.h"
#include "PractRand/Tests/NearSeq.h"
PractRand::Tests::ListOfTests testset_BirthdaySystematic() {
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayAlt(10), new PractRand::Tests::Birthday32());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayAlt(22));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdaySystematic128(28));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Birthday32());
	return PractRand::Tests::ListOfTests(new PractRand::Tests::Birthday64());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayLamda1(20));
}
PractRand::Tests::ListOfTests testset_experimental() {
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::FPMulti(3,0));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayAlt(10), new PractRand::Tests::Birthday32());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayAlt(22));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdaySystematic128(25));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Birthday32());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Birthday64());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayLamda1(20));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Rep16());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::FPMulti());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::NearSeq2());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::NearSeq3());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Gap16());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::LimitedBigGapPrototype(), new PractRand::Tests::Gap16());
	return PractRand::Tests::ListOfTests(new PractRand::Tests::LimitedBigGapPrototype());
	//return PractRand::Tests::ListOfTests(new Tests::mergewalk());
}
struct UnfoldedTestSet {
	int number;
	PractRand::Tests::ListOfTests(*callback)();
	const char *name;
};
UnfoldedTestSet test_sets[] = {
	{ 0, PractRand::Tests::Batteries::get_core_tests, "core" },//default value must come first
	{ 1, PractRand::Tests::Batteries::get_expanded_core_tests, "expanded" },
	{ 10, testset_BirthdaySystematic, "special (Birthday)" },
	{ 20, testset_experimental, "experimental" },
	{ -1, NULL, NULL }
};
int lookup_te_value(int te) {
	for (int i = 0; true; i++) {
		if (test_sets[i].number == te) return i;
		if (test_sets[i].number == -1) return -1;
	}
}
void print_rng_list(const char **list, int list_length, const char *name) {
	std::printf("%s = {\n", name);
	const int NAMES_PER_LINE = 4;
	const int NAME_SPACING = 18;
	const int LEADING_SPACES = 6;
	const bool USE_QUOTATION_MARKS = false;
	// characters printed per line is (at most): LEADING_SPACES + NAMES_PER_LINE * NAME_SPACING
	// unless one of the names doesn't fit in to the space provided...
	int so_far_in_line = 0;
	for (int index = 0; index < list_length; index++) {
		int len = std::strlen(list[index]);
		if (!so_far_in_line) std::printf("%*s", LEADING_SPACES, "");
		std::printf(USE_QUOTATION_MARKS ? "\"%s\"" : "%s", list[index]);
		if (index != list_length - 1) {
			std::printf(",");
			if (++so_far_in_line < NAMES_PER_LINE) {
				std::printf("%*s", NAME_SPACING - len - (USE_QUOTATION_MARKS ? 3 : 1), "");
			}
			else {
				std::printf("\n");
				so_far_in_line = 0;
			}
		}
		else {
			std::printf("\n");
			so_far_in_line = 0;
		}
	}
	std::printf("}\n");
}
void print_recommended_rng_list() {
	print_rng_list(PractRand::RNG_Sets::recommended_rngs, PractRand::RNG_Sets::num_recommended_rngs, "recommended_rngs");
	/*std::printf("recommended_rngs = {\n");
	int n = 0;
	const int NAMES_PER_LINE = 4;
	const int NAME_SPACING = 18;
	const int LEADING_SPACES = 6;
	int pos = 0;
	while (PractRand::RNG_Sets::recommended_rngs[n]) {
		const char *name = PractRand::RNG_Sets::recommended_rngs[n];
		int len = std::strlen(name);
		std::printf("%s%*s", name, NAME_WIDTH - len, "");
		if (len < NAME_WIDTH) len = NAME_WIDTH;
		pos += len;
		if (pos + NAME_WIDTH > LINE_WIDTH) {
			std::printf("\n");
			pos = 0;
		}
		n++;
	}
	if (pos) std::printf("\n");
	std::printf("}\n");*/
}
void print_reference_rng_list() {
	print_rng_list(PractRand::RNG_Sets::nonrecommended_simple, PractRand::RNG_Sets::num_nonrecommended_simple, "nonrecommended_rngs_simple");
	print_rng_list(PractRand::RNG_Sets::nonrecommended_oddball, PractRand::RNG_Sets::num_nonrecommended_oddball, "nonrecommended_rngs_oddball");
	print_rng_list(PractRand::RNG_Sets::nonrecommended_lcgish, PractRand::RNG_Sets::num_nonrecommended_lcgish, "nonrecommended_rngs_lcgish");
	print_rng_list(PractRand::RNG_Sets::nonrecommended_cbuf, PractRand::RNG_Sets::num_nonrecommended_cbuf, "nonrecommended_rngs_cbuf");
	print_rng_list(PractRand::RNG_Sets::nonrecommended_indirect, PractRand::RNG_Sets::num_nonrecommended_indirect, "nonrecommended_rngs_indirect");
	std::printf("... and nonrecommended_rngs includes all of the above, in order\n");
}

const char *get_exec_name(const char *argv0) {
	const char *p = std::strpbrk(argv0, "/\\");
	if (p) return get_exec_name(p + 1);
	else return argv0;
}

int main(int argc, char **argv) {
	PractRand::initialize_PractRand();
	const char *name = get_exec_name(argv[0]);
	std::printf("%s using PractRand version %s\n", name, PractRand::version_str);
#ifdef WIN32 // needed to allow binary stdin on windows
	_setmode( _fileno(stdin), _O_BINARY);
#endif
	using std::strcmp;
	if (argc <= 1) {
		std::printf("usage: %s RNG_name [options]   (runs tests on RNG_name)\n", name);
		std::printf("or: %s   (use no command line options to see this message)\n", name);
		std::printf("or: %s --help   (displays more instructions)\n", name);
		std::printf("or: %s --version   (displays version information)\n", name);
		std::printf("or: %s --self_test   (to verify PractRand was built correctly)\n", name);
		std::printf("or: %s --list_recommended_rngs   (lists recommended RNGs)\n", name);
		std::printf("or: %s --list_reference_rngs   (lists of all RNGs reference sets; these are bad RNGs used for comparison purposes)\n", name);
		std::printf("RNG_name can be the name of any PractRand recommended RNG (example: sfc16) or\n");
		std::printf("non-recommended RNG (example: mm32) or transformed RNG (exmple: SShrink(sfc16).\n");
		//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
		std::printf("Alternatively, use stdin as an RNG name to read raw binay data piped in from an\n");
		std::printf("external RNG.\n");
		std::printf("Options available include -a, -e, -p, -tf, -te, -ttnormal, -ttseed, -ttep,\n");
		std::printf("-tlmin, -tlmax, -tlshow, -tlfail, -tlmax_only, -multithreaded, -singlethreaded, and -seed.\n");
		std::printf("For more information run: %s --help\n\n", name);
		std::exit(0);
	}
	if (argv[1][0] == '-') {
		const char *cmd = argv[1] + 1;
		if (*cmd == '-') cmd++;
		if (false);
		else if (!strcmp(cmd, "version") || !strcmp(cmd, "v")) {
			std::printf("%s version %s\n", name, PractRand::version_str);
			// arbitrarily declaring the version number of RNG_test to match the version number of PractRand
			std::printf("A command line tool for testing RNGs with the PractRand library.\n");
			std::exit(0);
		}
		else if (!strcmp(cmd, "help") || !strcmp(cmd, "h")) {
			std::printf("syntax: %s RNG_name [options]\n", name);
			std::printf("or: %s   (use no command line options to see brief overview)\n", name);
			std::printf("or: %s --help   (to see this message)\n", name);
			std::printf("or: %s --version   (to see version number)\n", name);
			std::printf("or: %s --self_test   (to verify PractRand was built correctly)\n", name);
			std::printf("or: %s --list_recommended_rngs   (to see a list of all recommended rngs)\n", name);
			std::printf("or: %s --list_reference_rngs   (to see a list of all rngs currently in the reference set)\n", name);
			std::printf("A command line tool for testing RNGs with the PractRand library.\n");
			std::printf("RNG names:\n");
			std::printf("  To use an external RNG, use stdin as an RNG name and pipe in the random\n");
			std::printf("  numbers.  stdin8, stdin16, stdin32, and stdin64 also work, each interpretting\n");
			std::printf("  the input in slightly different ways.  Use stdin if you're uncertain how many\n");
			//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::printf("  bits the RNG produces at a time, or if it's not one of those options.\n");
			std::printf("  The lowest quality recommended RNGs are sfc16 and mt19937.\n");
			std::printf("  The entropy pooling RNGs available are arbee and sha2_basd_pool.\n");
			std::printf("  Small recommended RNGs include sfc16, sfc32, sfc64, jsf32, jsf64, .\n");
			std::printf("threshold options:\n");
			std::printf(" At most one threshold option should be specified.\n");
			//std::printf(" The default threshold setting is '-e 0.1', an alternative is '-p 0.001'\n");
			std::printf(" The default threshold setting is '-e 0.1', alternatives are '-p 0.001' or '-a'\n");
			std::printf("  -a             no threshold - display all test results.\n");
			std::printf("  -e EXPECTED    sets intelligent p-value thesholds to display an expected\n");
			std::printf("                 number of test results equal to EXPECTED.  If EXPECTED is zero\n");
			std::printf("                 or less then intelligent p-value thresholds will be disabled\n");
			std::printf("                 EXPECTED is a float with default value 0.1\n");
			std::printf("  -p THRESHOLD   sets simple p-value thresholds to display any test results\n");
			std::printf("                 within THRESHOLD of an extrema.  If THRESHOLD is zero or less \n");
			std::printf("                 then simple p-value thresholds will be disabled\n");
			std::printf("                 THRESHOLD is a float with recommended value 0.001\n");
			std::printf("test set options:\n");
			std::printf(" The default test set options are '-tf 1' and '-te 0'\n");
			std::printf("  -tf FOLDING    FOLDING may be 0, 1, or 2.  0 means that the base tests are \n");
			std::printf("                 run on only the raw test data.  1 means that the base tests \n");
			std::printf("                 are run on the raw test data and also on a simple transform \n");
			std::printf("                 that emphasizes the lowest bits.  2 means that the base tests\n");
			std::printf("                 are run on a wider variety of transforms of the test data.\n");
			std::printf("  -te EXPANDED   EXPANDED may be 0 or 1.  0 means that the base tests used are\n");
			std::printf("                 the normal ones for PractRand, optimized for sensitivity per \n");
			std::printf("                 time.  1 means that the expanded test set is used, optimized \n");
			std::printf("                 for sensitivity per bit.\n");
			std::printf("                 ... and now additional value(s) are supported.  Setting this \n");
			std::printf("                 to 10 will use a systematically expanding Birthday Spacings \n");
			std::printf("                 Test in place of a normal test set.  This test is separate \n");
			std::printf("                 because it uses too much memory to run concurrently with other\n");
			std::printf("                 tests\n");
			std::printf("test target options:\n");
			std::printf(" At most one test target option should be specified.\n");
			std::printf(" The default test target option is '-ttnormal'\n");
			std::printf("  -ttnormal      Test target: normal - the testing is done on the RNGs output.\n");
			std::printf("  -ttseed#       Test target: RNG seeding.  The \"#\" can be replaced with a \n");
			std::printf("                 number of bits from 16 to 128, or nothing (which causes it to \n");
			std::printf("                 pick a value based upon RNG metadata).  This repeatedly seeds \n");
			std::printf("                 First, the RNG is seeded with a randomly chosen seed, which is\n");
			//			 12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::printf("                 pick a value based upon RNG metadata).  This repeatedly seeds \n");
			std::printf("                 the RNG over and over again, generating 8 bytes out output \n");
			std::printf("                 from each seed.  The output is concatenated and sent to the \n");
			std::printf("                 tests as normal.  Each seed has a hamming distance of 1 from \n");
			std::printf("                 the previous seed used.  All seeds are restricted in size to \n");
			std::printf("                 the number of bits chosen (16 bits means seeds 0 to 65535). \n");
			std::printf("                 This effectively tests the quality of the seeding algorithm. \n");
			std::printf("                 Note that if this mode is enabled then the results header \n");
			std::printf("                 element listing log2 of number of bytes tested will instead \n");
			std::printf("                 list the log2 of number of seeds tested.\n");
			//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::printf("  -ttep          Test target: Entropy pooling.  This should only be done on\n");
			std::printf("                 RNGs that support entropy pooling.  It is similar to \n");
			std::printf("                 -ttseed, but the entropy accumulation methods are used\n");
			std::printf("                 instead of simple seeding, and the amount of entropy used may\n");
			std::printf("                 be much larger.\n");
			std::printf("                 Note that if this mode is enabled then the results header \n");
			std::printf("                 element listing log2 of number of bytes tested will instead \n");
			std::printf("                 list the log2 of number of strings tested.\n");
			std::printf("test length options:\n");
			std::printf("  -tlmin LENGTH  sets the minimum test length to LENGTH.  The tests will run on\n");
			std::printf("                 that much data before it starts printing regular results.  A\n");
			std::printf("                 large minimum will prevent it from displaying results on any\n");
			std::printf("                 test lengths other than the maximum length (set by tlmax) and\n");
			std::printf("                 lengths that were explicitly requested (by tlshow).\n");
			std::printf("                 See notes on lengths for details on how to express the length\n");
			std::printf("                 you want.\n");
			std::printf("                 The default minimum is 1.5 seconds (-tlmin 1.5s).\n");
			std::printf("  -tlmax LENGTH  sets the maximum test length to LENGTH.  The tests will stop\n");
			std::printf("                 after that much data.  See notes on lengths for details on how\n");
			std::printf("                 to express the length you want.\n");
			std::printf("                 The default maximum is 32 terabytes (-tlmax 32TB).\n");
			std::printf("  -tlshow LENGTH sets an additional point at which to display interim results.\n");
			std::printf("                 You can set multiple such points if desired.\n");
			std::printf("                 These are in addition to the normal interim results points,\n");
			std::printf("                 which are at every amount of data that is a power of 2 after\n");
			std::printf("                 the minimum and before the maximum.\n");
			std::printf("                 See the notes on lengths for details on how to express the\n");
			std::printf("                 lengths you want.\n");
			std::printf("  -tlfail        Halts testing after interim results are displayed if those\n");
			std::printf("                 results include any failures. (default)\n");
			std::printf("  -tlmaxonly     The opposite of -tlfail\n");
			std::printf("results table format option:\n");
			std::printf("  -table_wide       wide results tables with additional columns\n");
			std::printf("  -table_narrow     uses no more than 79 characters per line (default)\n");
			std::printf("  -table_scoring    uses <= 79 characters per line, emphasizes scoring\n");
			//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::printf("  -table_alterntive uses <= 79 characters per line, but experimental format\n");
			std::printf("other options:\n");
			std::printf("  -multithreaded  enables multithreaded testing.  Typically up to ~6 cores can\n");
			std::printf("                  be usefully used at once.\n");
			std::printf("  -no_metatest    disables experimental meta-test option.  (default)\n");
			std::printf("  -do_metatest    enables experimental meta-test option.  After each results \n");
			std::printf("                  summary a metatest will be performed on it.\n");
			//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::printf("  -singlethreaded disables multithreaded testing.  (default)\n");
			std::printf("  -seed SEED      specifies a 128 bit integer to seed the tested RNG with.  If\n");
			std::printf("                  no seed is specified then a seed will be chosen randomly.  \n");
			std::printf("                  The value must be expressed in hexadecimal.  An '0x' prefix\n");
			std::printf("                  on the seed is acceptable but not necessary.  If the RNG in\n");
			std::printf("                  use is incapable of using a seed (such as stdin), then any \n");
			std::printf("                  arbitrary string is acceptable.\n");
			std::printf("notes on lengths:\n");
			//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::printf("  Each of the test length options requires a field named LENGTH.  These fields\n");
			std::printf("  can accept either an amount of time or an amount of data.  In either case, \n");
			std::printf("  several types of units are supported.  \n");
			std::printf("  A time should be expressed as a number postfixed with either s, m, h, or d, \n");
			std::printf("  to express a number of seconds, minutes, hours, or days.  \n");
			std::printf("  example: -tlmin 1.4s (sets the minimum test length to 1.4 seconds)\n");
			std::printf("  An amount of data can be expressed as a number with no postfix, in which case\n");
			std::printf("  the number will be treated as the log-based-2 of the amount of bytes to test\n");
			std::printf("  (in normal target mode) or the log-baed-2 of the number of seeds or strings \n");
			std::printf("  to test in alternate test target modes.\n");
			std::printf("  example: -tlmin 23 (sets the minimum test length to 8 million bytes or 8\n");
			std::printf("    million seeds, depending upon test target mode)\n");
			std::printf("  Alternatively, an amount of data can be expressed as a number followed by \n");
			std::printf("  KB, MB, GB, TB, or PB for kilobytes, megabytes, gigabytes, or petabytes.\n");
			std::printf("  example: -tlmin 14KB (sets the minimum test length to 14 kilobytes\n");
			std::printf("  If the B is omitted on KB, MB, GB, TB, or PB then it treat the metric\n");
			std::printf("  prefixes as refering to numbers of bytes in normal test target mode, or\n");
			std::printf("  numbers of seeds in seeding test target mode, or numbers of strings in \n");
			std::printf("  entropy pooling test target mode.\n");
			std::printf("  example: -tlmin 40M (sets the minimum test length to ~40 million bytes or ~40\n");
			std::printf("    million seeds, depending upon test target mode)\n");
			std::printf("  A minor detail: I use the de facto standard (in which K means 1024 when\n");
			std::printf("  dealing with quantities of binary information) not the official standard (in\n");
			std::printf("  which K means 1000 no matter what is being dealt with unless an 'i' follows\n");
			std::printf("  the 'K').\n");
			//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
			std::exit(0);
		}
		else if (!strcmp(cmd, "self_test") || !strcmp(cmd, "selftest")) {
			PractRand::self_test_PractRand();
			std::printf("PractRand self-test passed.\n");
			std::exit(0);
		}
		else if (!std::strcmp(cmd, "list_recommended_rngs")) {
			print_recommended_rng_list();
			std::exit(0);
		}
		else if (!std::strcmp(cmd, "list_reference_rngs")) {
			print_reference_rng_list();
			std::exit(0);
		}
		else {
			std::fprintf(stderr, "\"%s\" is not valid as the first command line parameter to %s\n", argv[1], name);
			std::fprintf(stderr, "Valid parameters include various RNG names, --version, --help, --self_test, \n");
			std::fprintf(stderr, "--list_recommended_rngs, and --list_reference_rngs.\n");
			std::fprintf(stderr, "aborting\n");
			std::exit(1);
			//                    12345678901234567890123456789012345678901234567890123456789012345678901234567890
		}
	}

	RNG_Factories::register_recommended_RNGs();
	RNG_Factories::register_nonrecommended_RNGs();
	RNG_Factories::register_input_RNGs();
	RNG_Factories::register_candidate_RNGs();
	RNG_Factories::RNG_factory_index["dummy_rng"] = RNG_Factories::_generic_notrecommended_RNG_factory<dummy_rng>;
	//Seeder_MetaRNG::register_name();
	//FastSeeder_MetaRNG::register_name();
	FastSeeder128_MetaRNG::register_name();
	RecursiveSeed64_MetaRNG::register_name();
	EntropyPool_MetaRNG::register_name();
	std::string errmsg;
	std::unique_ptr<RNGs::vRNG> rng(RNG_Factories::create_rng(argv[1], &errmsg));
	if (!rng) {
		if (errmsg.empty()) std::fprintf(stderr, "unrecognized RNG name \"%s\".  aborting.\n", argv[1]);
		else std::fprintf(stderr, "%s\n", errmsg.c_str());
		std::exit(1);
	}
	bool do_selftest = true;
	bool use_multithreading = false;
	bool end_on_failure = true;
	bool smart_thresholds = true;
	bool do_metatests = false;
	int table_format = 0;
	double threshold = 0.1;
	int folding = 1;//0 = no folding, 1 = standard folding, 2 = extra folding
	int test_set_index = lookup_te_value(0);
	int seeding_bits = 128;//when in ttseed mode, this determines how many bits the seeds are restricted to - no impact on the values used in seed_low and seed_high
	int mode = 0;//0 = normal, 1 = test seeding, 2 = test entropy pooling
	for (int i = 2; i < argc; i++) {
		int params_left = argc - i - 1;
		//-a
		//-e EXPECTED
		//-p THRESHOLD
		if (false) ;
		else if (!std::strcmp(argv[i], "-a")) {
			smart_thresholds = false;
			threshold = 1;
		}
		else if (!std::strcmp(argv[i], "-e")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			smart_thresholds = true;
			threshold = std::atof(argv[++i]);
			if (threshold < 0.000001 || threshold > 1000) {
				std::printf("invalid smart threshold: -e %s (must be between 0.000001 and 1000)\n", argv[i]);
				std::exit(0);
			}
		}
		else if (!std::strcmp(argv[i], "-p")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			smart_thresholds = false;
			threshold = std::atof(argv[++i]);
			if (threshold < 0.0000000001 || threshold > 1.0) {
				std::printf("invalid p-value threshold: -p %s (must be between 0.0000000001 and 1)\n", argv[i]);
				std::exit(0);
			}
		}
		//-tf FOLDING
		//-te EXPANDED
		else if (!std::strcmp(argv[i], "-tf")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			folding = std::atoi(argv[++i]);
			if (folding < 0 || folding > 2) {
				std::printf("invalid folding test set value: -tf %s\n", argv[i]);
				std::exit(0);
			}
		}
		else if (!std::strcmp(argv[i], "-te")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			int expanded = std::atoi(argv[++i]);
			test_set_index = lookup_te_value(expanded);//0 maps to 0, but other values may not map to themselves
			if (test_set_index == -1) {
				std::printf("invalid expanded test set value: -te %s\n", argv[i]);
				std::exit(0);
			}
		}
		//-ttnormal
		//-ttseed
		//-ttep
		else if (!std::strcmp(argv[i], "-ttnormal")) mode = 0;
		else if (!std::strncmp(argv[i], "-ttseed", 7)) {
			mode = 1;
			seeding_bits = atoi(argv[i]+7);
			if (seeding_bits < 0) issue_error("-ttseed### - bits negative?");
		}
		else if (!std::strcmp(argv[i], "-ttep"))     mode = 2;
		//-tlmin LENGTH
		//-tlmax LENGTH
		//-tlshow LENGTH
		else if (!std::strcmp(argv[i], "-tlmin")) i++;
		else if (!std::strcmp(argv[i], "-tlmax")) i++;
		else if (!std::strcmp(argv[i], "-tlshow")) i++;
		//-tlfail
		//-tlmaxonly
		else if (!std::strcmp(argv[i], "-tlfail")) end_on_failure = true;
		else if (!std::strcmp(argv[i], "-tlmaxonly")) end_on_failure = false;

		//-threads
		//-nothreads
		//-seed SEED
		else if (!std::strcmp(argv[i], "-multithreaded")) use_multithreading = true;
		else if (!std::strcmp(argv[i], "-singlethreaded")) use_multithreading = false;
		else if (!std::strcmp(argv[i], "-do_metatest")) do_metatests = true;
		else if (!std::strcmp(argv[i], "-no_metatest")) do_metatests = false;
		else if (!std::strcmp(argv[i], "-table_narrow")) table_format = 0;
		else if (!std::strcmp(argv[i], "-table_wide")) table_format = 1;
		else if (!std::strcmp(argv[i], "-table_alternative")) table_format = 2;
		else if (!std::strcmp(argv[i], "-table_scoring")) table_format = 3;
		else if (!std::strcmp(argv[i], "-seed")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			seed_str = argv[++i];
		}
		else if (!std::strcmp(argv[i], "-skip_selftest")) do_selftest = false;
		else {
			std::fprintf(stderr, "unrecognized testing option: %s\naborting\n", argv[i]);
			std::exit(1);
		}
	}
	if (do_selftest) PractRand::self_test_PractRand();

#if !defined MULTITHREADING_SUPPORTED
	if (use_multithreading) {
		std::printf("multithreading is not supported on this build.  If multithreading should be supported, try defining the MULTITHREADING_SUPPORTED preprocessor symbol during the build process.\n");
		std::exit(0);
	}
#endif
	enum {TL_MIN = 1, TL_MAX = 2, TL_SHOW = 3};//, TL_FLAG_UNSPECIFIED_UNITS = 16};
	std::map<double,int> show_times;
	std::map<Uint64,int> show_datas;
	double show_min = -2.0;
	double show_max = 1ull << 45;
	//walking parameters a second time to force the mode to be known prior to finding the test lengths
	for (int i = 2; i < argc; i++) {
		int params_left = argc - i - 1;
		if (false) ;
		else if (!std::strcmp(argv[i], "-tlmin")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			double length = interpret_length(argv[++i], !mode);
			if (length >= 0 && length < Tests::TestBlock::SIZE) { std::printf("invalid test length: %s\n", argv[i]); std::exit(0); }
			show_min = length;
		}
		else if (!std::strcmp(argv[i], "-tlmax")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			double length = interpret_length(argv[++i], !mode);
			if (length >= 0 && length < Tests::TestBlock::SIZE) { std::printf("invalid test length: %s\n", argv[i]); std::exit(0); }
			show_max = length;
		}
		else if (!std::strcmp(argv[i], "-tlshow")) {
			if (params_left < 1) {std::printf("command line option %s must be followed by a value\n", argv[i]); std::exit(0);}
			double length = interpret_length(argv[++i], !mode);
			if (length >= 0 && length < Tests::TestBlock::SIZE) { std::printf("invalid test length: %s\n", argv[i]); std::exit(0); }
			if (length < 0) show_times[-length] = TL_SHOW;
			else show_datas[Uint64(length) / Tests::TestBlock::SIZE] = TL_SHOW;
		}
	}
	if (show_min < 0) show_times[-show_min] = TL_MIN;
	else show_datas[Uint64(show_min) / Tests::TestBlock::SIZE] = TL_MIN;
	if (show_max < 0) show_times[-show_max] = TL_MAX;
	else show_datas[Uint64(show_max) / Tests::TestBlock::SIZE] = TL_MAX;

	std::time_t start_time = std::time(NULL);
	TimeUnit start_clock = get_time();

	seed_low = known_good.raw32();//128 bit space, as that's what the interface accepts, but 32 bit random value so that by default it's not too onerous to record/compare/whatever the value by hand
	seed_high = 0;
	if (seed_str && !(rng->get_flags() & PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED)) {
		if (!interpret_seed(seed_str, seed_low, seed_high)) {
			std::printf("\"%s\" is not a valid 128 bit hexadecimal seed\n", seed_str);
			std::exit(0);
		}
	}
	known_good.seed(seed_low, seed_high);//it's possible the RNG uses the same algorithm as the known good RNG, but if so it's a very good algorithm anyway, and this doesn't matter much

	//PractRand::RNGs::vRNG *testing_rng;
	if (mode == 0) {
		rng->seed(seed_low, seed_high);
		//testing_rng = rng;
	}
	else if (mode == 1) {
		//it would be nice to print a warning here for RNGs that use generic integer seeding
		//but that's a little difficult atm as there's no way to query whether an RNG does so

		//rng.reset(new Seeder_MetaRNG(rng.release()));
		if (!seeding_bits) {
			Uint64 low, high;
			rng->get_maximum_seed(low, high);
			if (!low && !high) {//no metadata available
				seeding_bits = 128;
				std::printf("warning - -ttseed enabled, seed size autodetected, no seed size metadata available, defaulting to 128 bit seeds\n");
			}
			else if (high) seeding_bits = 64 + std::floor(std::log2(1.0 + high) + 0.000001);
			else seeding_bits = 64 + std::floor(std::log2(1.0 + low) + 0.000001);//rounding up, but only when VERY close to the next bit
			// there are many cases where there are for example 2**47.99999999999999 valid seeds
			// and the cost of accidentally duplicating a seed once or twice isn't a big deal, so rounding up is desirable
			// BUT, there are also wide swaths between integer numbers of bits where rounding up would result in large percentages of seeds being duplicates
			// so the compromise is, round up in some cases, but only when absurdly close already
		}
		rng.reset(new FastSeeder128_MetaRNG(rng.release(), seeding_bits));
		//rng.reset(new FastSeeder_MetaRNG(rng.release(), seeding_bits));
		//testing_rng = new Seeder_MetaRNG(rng);
		//testing_rng->seed(seed);
		rng->seed(seed_low, seed_high);
	}
	else if (mode == 2) {
		if (!(rng->get_flags() & PractRand::RNGs::FLAG::SUPPORTS_ENTROPY_ACCUMULATION)) {
			std::printf("Entropy pooling is not supported by this RNG, so mode --ttep is invalid.\n");
			std::printf("aborting\n");
			std::exit(0);
		}
		rng->reset_entropy();
		Uint64 a = rng->raw64();
		rng->reset_entropy();
		Uint64 b = rng->raw64();
		if (a != b) {
			std::printf("entropy pooling RNG \"%s\" failed basic check 1.\naborting\n", rng->get_name().c_str());
			std::exit(0);
		}
		Uint64 s64 = known_good.raw64();
		rng->reset_entropy();
		rng->add_entropy64(s64);
		Uint64 c1 = rng->raw64();
		Uint64 c2 = rng->raw64();
		rng->reset_entropy();
		rng->add_entropy64(s64);
		Uint64 d = rng->raw64();
		rng->reset_entropy();
		rng->add_entropy64(s64+1);
		Uint64 e1 = rng->raw64();
		Uint64 e2 = rng->raw64();
		if (c1 != d) {
			std::printf("entropy pooling RNG \"%s\" failed basic check 2.\naborting\n", rng->get_name().c_str());
			std::exit(0);
		}
		if (c1 == e1 && c2 == e2) {
			std::printf("entropy pooling RNG \"%s\" probably failed basic check 3.\naborting\n", rng->get_name().c_str());
			std::exit(0);
		}
		//rng->seed(seed);
		//	I'd like to test varying length entropy strings, but known good EPs are failing eventually when varying length is allowed for some reason
		//rng.reset(new EntropyPool_MetaRNG(rng.release(), 48, 64));
		rng.reset(new EntropyPool_MetaRNG(rng.release(), 57, 57));
		//testing_rng->seed(seed);
		rng->seed(seed_low, seed_high);
	}
	else {
		std::printf("invalid mode, aborting\n");
		std::exit(1);
	}

	std::printf("RNG = %s, seed = ", rng->get_name().c_str());
	if (rng->get_flags() & PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED) {
		if (seed_str) std::printf("%s", seed_str);
		else std::printf("unknown");
	}
	else {
		print_128bit_hex(seed_low, seed_high);
	}
	const char *folding_names[3] = {"none", "standard", "extra"};
	std::printf("\ntest set = %s, folding = %s", test_sets[test_set_index].name, folding_names[folding]);
	if (folding == 1) {
		int native_bits = rng->get_native_output_size();
		if (native_bits > 0) std::printf(" (%d bit)", native_bits);
		else std::printf("(unknown format)");
	}

	std::printf("\n\n");

	Tests::ListOfTests tests( (Tests::TestBaseclass*)NULL);
	if (test_set_index == -1) { std::printf("internal error\n"); std::exit(1); }
	if (false) ;
	else if (folding == 0) tests = test_sets[test_set_index].callback();
	else if (folding == 1) tests = Tests::Batteries::apply_standard_foldings(rng.get(), test_sets[test_set_index].callback);
	else if (folding == 2) tests = Tests::Batteries::apply_extended_foldings(test_sets[test_set_index].callback);
	else { std::printf("internal error\n"); std::exit(1); }

//	Tests::ListOfTests tests = Tests::Batteries::get_expanded_standard_tests(rng);
#if defined MULTITHREADING_SUPPORTED
	std::unique_ptr<TestManager> tman;
	if (use_multithreading) tman.reset( new MultithreadedTestManager(&tests, &known_good));
	else tman.reset( new TestManager(&tests, &known_good));
#else
	std::unique_ptr<TestManager> tman(new TestManager(&tests, &known_good));
#endif
	tman->reset(rng.get());

	Uint64 blocks_tested = 0;
	bool already_shown = false;
	Uint64 next_power_of_2 = 1;
	bool showing_powers_of_2 = false;
	double time_passed = 0;
	while (true) {
		Uint64 blocks_to_test = next_power_of_2 - blocks_tested;
		enum {MAX_BLOCKS = 256 * 1024};
		if (blocks_to_test > MAX_BLOCKS) blocks_to_test = MAX_BLOCKS;
		while (!show_datas.empty()) {
			Uint64 data_checkpoint = show_datas.begin()->first - blocks_tested;
			if (data_checkpoint) {
				if (data_checkpoint < blocks_to_test) blocks_to_test = data_checkpoint;
				break;
			}
			int action = show_datas.begin()->second;
			if (action == TL_SHOW) {
				if (!already_shown) show_checkpoint(*tman, mode, time_passed, smart_thresholds, threshold, end_on_failure, do_metatests, table_format);
				already_shown = true;
			}
			else if (action == TL_MIN) {
				showing_powers_of_2 = true;
				if (!already_shown) show_checkpoint(*tman, mode, time_passed, smart_thresholds, threshold, end_on_failure, do_metatests, table_format);
				already_shown = true;
			}
			else if (action == TL_MAX) {
				if (!already_shown) show_checkpoint(*tman, mode, time_passed, smart_thresholds, threshold, end_on_failure, do_metatests, table_format);
				return 0;
			}
			else {std::printf("internal error: unrecognized test length code, aborting\n");std::exit(1);}
			show_datas.erase(show_datas.begin());
		}
		while (!show_times.empty()) {
			double time_checkpoint = show_times.begin()->first - time_passed;
			if (time_checkpoint > 0) break;
			int action = show_times.begin()->second;
			show_times.erase(show_times.begin());
			if (action == TL_SHOW) {
				if (!already_shown) show_checkpoint(*tman, mode, time_passed, smart_thresholds, threshold, end_on_failure, do_metatests, table_format);
				already_shown = true;
			}
			else if (action == TL_MIN) showing_powers_of_2 = true;
			else if (action == TL_MAX) {
				if (!already_shown) show_checkpoint(*tman, mode, time_passed, smart_thresholds, threshold, end_on_failure, do_metatests, table_format);
				return 0;
			}
			else {std::printf("internal error: unrecognized test length code, aborting\n");std::exit(1);}
		}

		if (blocks_tested == next_power_of_2) {
			if (showing_powers_of_2) {
				if (!already_shown) show_checkpoint(*tman, mode, time_passed, smart_thresholds, threshold, end_on_failure, do_metatests, table_format);
				already_shown = true;
			}
			next_power_of_2 <<= 1;
			continue;
		}
		tman->test(blocks_to_test);
		blocks_tested += blocks_to_test;
		already_shown = false;

		double clocks_passed = TimeUnit(get_time() - start_clock) * get_time_period();//may wrap too quickly
		int seconds_passed = std::time(NULL) - start_time;
		if (seconds_passed >= 1000 || seconds_passed > clocks_passed + 2.0) time_passed = seconds_passed;
		else time_passed = clocks_passed;
	}

	return 0;
}



