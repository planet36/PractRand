
#define MULTITHREADING_SUPPORTED

//master header, includes everything in PractRand for both
//  practical usage and research...
//  EXCEPT it does not include specific algorithms
#include "PractRand_full.h"

//specific algorithms: all recommended RNGs
#include "PractRand/RNGs/all.h"

//specific algorithms: non-recommended RNGs
#include "PractRand/RNGs/other/fibonacci.h"
#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/other/mult.h"
#include "PractRand/RNGs/other/simple.h"
#include "PractRand/RNGs/other/special.h"
#include "PractRand/RNGs/other/transform.h"

//for access to some functions the dummy PRNG might want to use:
#include "PractRand/rng_internals.h"

//tests used by the special and experimental test sets:
#include "PractRand/Tests/Birthday.h"
#include "PractRand/Tests/DistFreq4.h"
#include "PractRand/Tests/FPF.h"
#include "PractRand/Tests/FPMulti.h"
#include "PractRand/Tests/Gap16.h"

//helpers for the test programs, to deal with RNG names, test usage, etc
#include "RNG_from_name.h"
#include "parse_number.h"

#include "Candidate_RNGs.h"
#include "SeedingTester.h"
#include "TestManager.h"
#ifdef MULTITHREADING_SUPPORTED
#include "MultithreadedTestManager.h"
#endif

#include <bit>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <map>
#include <numbers>
#include <print>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

using namespace PractRand;

PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);

//using TimeUnit = std::chrono::system_clock::rep;
//TimeUnit get_time() { return std::chrono::system_clock::now().time_since_epoch().count(); }
//double get_time_period() { return std::chrono::system_clock::period::num / static_cast<double>(std::chrono::system_clock::period::den); }

/*
A minimal RNG implementation, just enough to make it usable.
Deliberately flawed, though still better than many platforms default RNGs
*/

class DummyRNG : public PractRand::RNGs::vRNG16 {
public:
	//declare state
	//Uint16 s1, s2, s3, s4;
	PractRand::RNGs::Polymorphic::NotRecommended::lcg16of64_varqual rng1;
	PractRand::RNGs::Polymorphic::NotRecommended::simpleB rng2;
	//and any helper methods you want:
	static Uint16 ddRot16(Uint16 value) { return std::rotr(value, (value >> (16 - 3)) << 1); }
	static Uint32 ddRot32(Uint32 value) { return std::rotr(value, (value >> (32 - 3)) << 2); }
	static Uint64 ddRot64(Uint64 value) { return std::rotr(value, (value >> (64 - 4)) << 2); }
	//constructor, if necessary
	DummyRNG() : rng1(8) {}
	//implement algorithm
	Uint16 raw16() override {
		Uint16 v1 = rng1.raw16();
		Uint16 v2 = rng2.raw16();
		Uint16 x = v1 ^ v2;
		return x;
		//Uint16 old = s4;
		//s4 += s3; s3 += s2; s2 += s1; s1 += s4;
		//s4 = s3; s3 = s2; s2 = s1; s1 += old;
		//Uint16 a = ddRot16(s1) - s2;
		//Uint16 b = ddRot16(s3) - s4;
		//Uint16 c = ddRot16(s2) - ddRot16(s4);
		//return ddRot16((s1 + s3) * 1) + ((s2 * 9) ^ std::rotl(s2 * 9, 5)) + (s4 * 3);
		//return (a & b) | (c & ~a);
		//return (a & b) | (b & c) | (a & c);
		//return a ^ b;
		//return (a * 1) + (b * 1) + (c * 1);
		//return (a * 11) + (b * 13) + (c * 15);
	}
	//allow PractRand to be aware of your internal state
	//uses include: default seeding mechanism (individual PRNGs can override), state serialization/deserialization, maybe eventually some avalanche testing tools
	void walk_state(PractRand::StateWalkingObject* walker) override {
		rng1.walk_state(walker);
		rng2.walk_state(walker);
		//walker->handle(s1);
		//walker->handle(s2);
		//walker->handle(s3);
		//walker->handle(s4);
	}
	//seeding from integers
	//not actually necessary, in the absence of such a method a default seeding-from-integer path will use walk_state to randomize the member variables
	//note that a separate path exists for seeding-from-another-PRNG
	void seed(Uint64 sv) override {
		rng1.seed(sv);
		rng2.seed(sv);
//		s1 = s2 = s3 = s4 = sv;
	}
	using vRNG::seed;
	//any name you want
	[[nodiscard]] std::string get_name() const override {return "DummyRNG";}
};
/*
	The above class is enough to create a PRNG compatible with PractRand.
	You can pass that to any non-template function in PractRand expecting a polymorphic RNG and it should work fine.
	(some template functions in PractRand require additional metadata or other weirdness)
	HOWEVER, that is not enough to allow this (or any other) command line tool to recognize the name of your RNG from the command line.
	For that, search for the line mentioning RNG_factory_index["dummy"] in the "main" function below here,
	that line allows it to recognize "dummy" on the command line as corresponding to this class.
*/

double print_result(const PractRand::TestResult& result, bool print_header = false) {
	if (print_header) std::println("  Test Name                         Raw       Processed     Evaluation");
	//                                     10        20        30        40        50        60        70        80
	std::print("  ");//2 characters
	//NAME
	if constexpr (true) {// 34 characters
		std::print("{}", result.name);
		int len = result.name.length();
		for (int i = len; i < 34; i++) std::print(" ");
	}

	//RAW TEST RESULT
	if constexpr (true) {// 10 characters?
		double raw = result.get_raw();
		if (raw > 99999.0) std::print("R>+99999  ");
		else if (raw < -99999.0) std::print("R<-99999  ");
		else if (std::abs(raw) < 999.95) std::print("R={:+6.1f}  ", raw);
		else std::print("R={:+6.0f}  ", raw);
		//if (std::abs(raw) < 99999.5) std::printf(" ");
		//if (std::abs(raw) < 999999.5) std::printf(" ");
		//if (std::abs(raw) < 9999999.5) std::printf(" ");
	}

	//RESULT AS A NUMERICAL "SUSPICION LEVEL" (log of distance from pvalue to closest extrema)
	if constexpr (false) {// 12 characters?
		bool printed = false;
		double susp = result.get_suspicion();
		if (result.type == PractRand::TestResult::TYPE_PASSFAIL) {
			std::print("  {}    ", result.get_pvalue() ? "\"pass\"" : "\"fail\"");
		}
		else if (result.type == PractRand::TestResult::TYPE_RAW) {
			std::print("            ");
		}
		else if (result.type == PractRand::TestResult::TYPE_BAD_P || result.type == PractRand::TestResult::TYPE_BAD_S || result.type == PractRand::TestResult::TYPE_RAW_NORMAL) {
			std::print("S=~{:+6.1f}   ", susp);
			printed = true;
		}
		else {
			std::print("S ={:+6.1f}   ", susp);
			printed = true;
		}
		if (printed) {
			if (std::abs(susp) < 9999.95) std::print(" ");
			if (std::abs(susp) < 999.95) std::print(" ");
		}
	}

	//RESULT AS A p-value
	if constexpr (true) {// 14 characters?
		if (result.type == PractRand::TestResult::TYPE_PASSFAIL) {
			std::print("  {}      ", result.get_pvalue() ? "\"pass\"" : "\"fail\"");
		}
		else if (result.type == PractRand::TestResult::TYPE_RAW) {
			std::print("              ");
		}
		else if (result.type == PractRand::TestResult::TYPE_BAD_P || result.type == PractRand::TestResult::TYPE_GOOD_P || result.type == PractRand::TestResult::TYPE_RAW_NORMAL) {
			double p = result.get_pvalue();
			double a = std::abs(p-0.5);
			std::print("{}", (result.type != PractRand::TestResult::TYPE_GOOD_P) ? "p~= " : "p = ");
			if (a > 0.49) {
				double s = result.get_suspicion();
				double ns = std::abs(s) + 1;
				double dec = ns / (std::numbers::ln10 / std::numbers::ln2);
				double dig = std::ceil(dec);
				double sig = std::floor(std::pow(0.1, dec - dig));
				if (dig > 999) { std::print(" {}        ", (s > 0) ? 1 : 0); }
				else {
					if (s > 0) std::print("1-{:1.0f}e-{:.0f}  ", sig, dig);
					else       std::print("  {:1.0f}e-{:.0f}  ", sig, dig);
					if (dig < 100) std::print(" ");
					if (dig < 10) std::print(" ");
				}
			}
			else if (result.type == PractRand::TestResult::TYPE_GOOD_P) { std::print("{:5.3f}     ", p); }
			else if (a >= 0.4) {                          std::print("{:4.2f}      ", p); }
			else {                                        std::print("{:3.1f}       ", p); }
		}
		else if (result.type == PractRand::TestResult::TYPE_BAD_S || result.type == PractRand::TestResult::TYPE_GOOD_S) {
			double s = result.get_suspicion();
			double p = result.get_pvalue();
			std::print("{}", (result.type == PractRand::TestResult::TYPE_BAD_S || result.type == PractRand::TestResult::TYPE_RAW_NORMAL) ? "p~=" : "p =");
			if (p >= 0.01 && p <= 0.99) { std::print(" {:.3f}     ", p); }
			else {
				double ns = std::abs(s) + 1;
				double dec = ns / (std::numbers::ln10 / std::numbers::ln2);
				double dig = std::ceil(dec);
				double sig = std::pow(0.1, dec - dig);
				sig = std::floor(sig * 10) * 0.1;
				if (dig > 9999) { std::print(" {}         ", (s > 0) ? 1 : 0); }
				else if (dig > 999) {
					sig = std::floor(sig);
					if (s > 0) std::print("1-{:1.0f}e-{:.0f}  ", sig, dig);
					else       std::print("  {:1.0f}e-{:.0f}  ", sig, dig);
				}
				else {
					if (s > 0) std::print("1-{:3.1f}e-{:.0f} ", sig, dig);
					else       std::print("  {:3.1f}e-{:.0f} ", sig, dig);
					if (dig < 100) std::print(" ");
					if (dig < 10) std::print(" ");
				}
			}
		}
	}

	double dec = std::numbers::ln2 / std::numbers::ln10;
	double as = (std::abs(result.get_suspicion()) + 1 - 1) * dec;// +1 for suspicion conversion, -1 to account for there being 2 failure regions (near-zero and near-1)
	double wmod = std::log(result.get_weight()) / std::log(0.5) * dec;
	double rs = as - wmod;
	//double ap = std::abs(0.5 - result.get_pvalue());
	//MESSAGE DESCRIBING RESULT IN ENGLISH
	if constexpr (true) {// 17 characters?
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
		if constexpr (false) ;
		else if (rs >999) std::print("  FAIL !!!!!!!!  ");
		else if (rs >325) std::print("  FAIL !!!!!!!   ");
		else if (rs >165) std::print("  FAIL !!!!!!    ");
		else if (rs > 85) std::print("  FAIL !!!!!     ");
		else if (rs > 45) std::print("  FAIL !!!!      ");
		else if (rs > 25) std::print("  FAIL !!!       ");
		else if (rs > 17) std::print("  FAIL !!        ");
		else if (rs > 12) std::print("  FAIL !         ");
		else if (rs >8.5) std::print("  FAIL           ");
		else if (rs >6.0) std::print(" VERY SUSPICIOUS ");
		else if (rs >4.0) std::print("very suspicious  ");
		else if (rs >3.0) std::print("suspicious       ");
		else if (rs >2.0) std::print("mildly suspicious");
		else if (rs >1.0) std::print("unusual          ");
		else if (rs >0.0) std::print("normalish       ");
		else              std::print("normal           ");
	}
	std::println("");
	return rs;
}

const char* seed_str = nullptr;

void show_checkpoint(TestManager* tman, int mode, Uint64 seed, double time, bool smart_thresholds, double threshold, bool end_on_failure) {
	std::print("rng={}", tman->get_rng()->get_name());

	std::print(", seed=");
	if (tman->get_rng()->get_flags() & PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED) {
		if (seed_str) std::print("{}", seed_str);
		else std::print("unknown");
	}
	else {
		if (seed >> 32) std::print("0x{:x}{:08x}", long(seed >> 32), long(seed & 0xFFffFFff));
		else std::print("0x{:x}", long(seed));
	}
	std::println("");

	std::print("length= ");
	Uint64 length = tman->get_blocks_so_far() * Tests::TestBlock::SIZE;
	double log2b = std::log(double(length)) / std::numbers::ln2;
	const char* unitstr[6] = {"kibibyte", "mebibyte", "gibibyte", "tebibyte", "pebibyte", "exbibyte"};
	int units = int(std::floor(log2b / 10)) - 1;
	if (units < 0 || units > 5) {std::println("internal error: length out of bounds?");std::exit(1);}
	if (length & (length-1))
		std::print("{:.3f} {}s", length * std::pow(0.5,units*10.0+10), unitstr[units] );
	else std::print("{:.0f} {}{}", length * std::pow(0.5,units*10.0+10), unitstr[units], length != (Uint64(1024)<<(units*10)) ? "s" : "" );
	if (length & (length-1)) std::print(" (2^{:.3f}", log2b - (mode?3:0)); else std::print(" (2^{:.0f}", log2b - (mode?3:0));
	const char* mode_unit_names[3] = {"bytes", "seeds", "entropy strings"};
	std::print(" {}), time= ", mode_unit_names[mode]);
	if (time < 99.95) std::println("{:.1f} seconds", time);
	else std::println("{:.0f} seconds", time);

	std::vector<PractRand::TestResult> results;
	tman->get_results(results);
	double total_weight = 0, min_weight = 9999999;
	for (const auto& result : results) {
		double weight = result.get_weight();
		total_weight += weight;
		if (weight < min_weight) min_weight = weight;
	}
	if (min_weight <= 0) {
		std::println("error: result weight too small");
		std::exit(1);
	}
	std::vector<int> marked;
	for (unsigned int i = 0; i < results.size(); i++) {
		results[i].set_weight(results[i].get_weight() / total_weight);
		if (!smart_thresholds) {
			if (std::abs(0.5 - results[i].get_pvalue()) < 0.5 - threshold) continue;
		}
		else {
			double T = threshold * results[i].get_weight() * 0.5;
			if (std::abs(0.5 - results[i].get_pvalue()) < (0.5 - T)) continue;
		}
		marked.push_back(i);
	}
	double biggest_decimal_suspicion = 0;
	for (unsigned int i = 0; i < marked.size(); i++) {
		double decimal_suspicion = print_result(results[marked[i]], i == 0);
		if (decimal_suspicion > biggest_decimal_suspicion) biggest_decimal_suspicion = decimal_suspicion;
	}
	if (marked.size() == results.size())
		;
	else if (marked.empty())
		std::println("  no anomalies in {} test result(s)", int(results.size()));
	else
		std::println("  ...and {} test result(s) without anomalies", int(results.size() - marked.size()));
	std::println("");
	(void)std::fflush(stdout);
	if (end_on_failure && biggest_decimal_suspicion > 8.5) {
		std::exit(0);
	}
}
double interpret_length(const std::string& lengthstr, bool normal_mode) {
	//(0-9)*[.(0-9)*][((K|M|G|T|P)[B])|(s|m|h|d)]
	int mode_factor = normal_mode ? 1 : 8;
	unsigned int pos = 0;
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
	default:
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
bool interpret_seed(const std::string& seedstr, Uint64& seed) {
	const char* first = seedstr.data();
	const char* last = first + seedstr.size();
	if (seedstr.starts_with("0x")) first += 2;
	auto [ptr, ec] = std::from_chars(first, last, seed, 16);
	return ec == std::errc() && ptr == last;
}

PractRand::Tests::ListOfTests testset_BirthdaySystematic() {
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayAlt(10), new PractRand::Tests::Birthday32());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdayAlt(22));
	return PractRand::Tests::ListOfTests(new PractRand::Tests::BirthdaySystematic128(26));
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Birthday32());
	//return PractRand::Tests::ListOfTests(new PractRand::Tests::Birthday64());
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
	return PractRand::Tests::ListOfTests(new Tests::FPF(4, 14, 6));
}
struct UnfoldedTestSet {
	int number;
	PractRand::Tests::ListOfTests(*callback)();
	const char* name;
};
UnfoldedTestSet test_sets[] = {
	{ .number=0, .callback=PractRand::Tests::Batteries::get_core_tests, .name="core" },//default value must come first
	{ .number=1, .callback=PractRand::Tests::Batteries::get_expanded_core_tests, .name="expanded" },
	{ .number=10, .callback=testset_BirthdaySystematic, .name="special (Birthday)" },
	{ .number=20, .callback=testset_experimental, .name="experimental" },
	{ .number=-1, .callback=nullptr, .name=nullptr }
};
int lookup_te_value(int te) {
	for (int i = 0; true; i++) {
		if (test_sets[i].number == te) return i;
		if (test_sets[i].number == -1) return -1;
	}
}

int main(int argc, char** argv) { // NOLINT(bugprone-exception-escape)
	PractRand::initialize_PractRand();
	PractRand::hook_error_handler(PractRand::print_err);
	std::println("RNG_test using PractRand version {}", PractRand::version_str);
	if (argc <= 1) {
		std::println("usage: {} RNG_name [options]  --  runs tests on RNG_name", argv[0]);
		std::println("or: {} -help  --  displays more instructions", argv[0]);
		std::println("or: {} -version  --  displays version information", argv[0]);
		std::println("RNG_name can be the name of any PractRand recommended RNG (example: sfc16) or");
		std::println("non-recommended RNG (example: mm32) or transformed RNG (exmple: SShrink(sfc16).");
		//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
		std::println("Alternatively, use stdin as an RNG name to read raw binay data piped in from an");
		std::println("external RNG.");
		std::println("options available include -a, -e, -p, -tf, -te, -ttnormal, -ttseed64, -ttep,");
		std::println("-tlmin, -tlmax, -tlshow, -multithreaded, -singlethreaded, and -seed.");
		std::println("For more information run: {} -help\n", argv[0]);
		std::exit(0);
	}
	if (!strcmp(argv[1], "-version") || !strcmp(argv[1], "--version") || !strcmp(argv[1], "-v")) {
		std::println("RNG_test version {}", PractRand::version_str);
		// arbitrarily declaring the version number of RNG_test to match the version number of PractRand
		std::println("A command line tool for testing RNGs with the PractRand library.");
		std::exit(0);
	}
	if (!strcmp(argv[1], "-help") || !strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
		std::println("syntax: {} RNG_name [options]", argv[0]);
		std::println("or: {} -help (to see this message)", argv[0]);
		std::println("or: {} -version (to see version number)", argv[0]);
		std::println("A command line tool for testing RNGs with the PractRand library.");
		std::println("RNG names:");
		std::println("  To use an external RNG, use stdin as an RNG name and pipe in the random");
		std::println("  numbers.  stdin8, stdin16, stdin32, and stdin64 also work, each interpretting");
		std::println("  the input in slightly different ways.  Use stdin if you're uncertain how many");
		//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
		std::println("  bits the RNG produces at a time, or if it's not one of those options.");
		std::println("  The lowest quality recommended RNGs are sfc16 and mt19937.");
		std::println("  The entropy pooling RNGs available are arbee and sha2_basd_pool.");
		std::println("  Small recommended RNGs include sfc16, sfc32, sfc64, jsf32, jsf64, .");
		std::println("threshold options:");
		std::println(" At most one threshold option should be specified.");
		//std::printf(" The default threshold setting is '-e 0.1', an alternative is '-p 0.001'\n");
		std::println(" The default threshold setting is '-e 0.1', alternatives are '-p 0.001' or '-a'");
		std::println("  -a             no threshold - display all test results.");
		std::println("  -e EXPECTED    sets intelligent p-value thesholds to display an expected");
		std::println("                 number of test results equal to EXPECTED.  If EXPECTED is zero");
		std::println("                 or less then intelligent p-value thresholds will be disabled");
		std::println("                 EXPECTED is a float with default value 0.1");
		std::println("  -p THRESHOLD   sets simple p-value thresholds to display any test results");
		std::println("                 within THRESHOLD of an extrema.  If THRESHOLD is zero or less");
		std::println("                 then simple p-value thresholds will be disabled");
		std::println("                 THRESHOLD is a float with recommended value 0.001");
		std::println("test set options:");
		std::println(" The default test set options are '-tf 1' and '-te 0'");
		std::println("  -tf FOLDING    FOLDING may be 0, 1, or 2.  0 means that the base tests are");
		std::println("                 run on only the raw test data.  1 means that the base tests");
		std::println("                 are run on the raw test data and also on a simple transform");
		std::println("                 that emphasizes the lowest bits.  2 means that the base tests");
		std::println("                 are run on a wider variety of transforms of the test data.");
		std::println("  -te EXPANDED   EXPANDED may be 0 or 1.  0 means that the base tests used are");
		std::println("                 the normal ones for PractRand, optimized for sensitivity per");
		std::println("                 time.  1 means that the expanded test set is used, optimized");
		std::println("                 for sensitivity per bit.");
		std::println("                 ... and now additional value(s) are supported.  Setting this");
		std::println("                 to 10 will use an systematically expanding Birthday Spacings");
		std::println("                 Test in place of a normal test set.  This test is separate");
		std::println("                 because it uses too much memory to run concurrently with other");
		std::println("                 tests");
		std::println("test target options:");
		std::println(" At most one test target option should be specified.");
		std::println(" The default test target option is '-ttnormal'");
		std::println("  -ttnormal      Test target: normal - the testing is done on the RNGs output.");
		std::println("  -ttseed64      Test target: RNG seeding from 64 bit integers.  First, the RNG");
		std::println("                 is seeded with a randomly chosen 64 bit integer.  Then 8 bytes");
		std::println("                 of output are taken from the RNG and given to the tests.  Then");
		std::println("                 another seed is chosen at a low hamming distance from the");
		std::println("                 prior seed and another 8 bytes of RNG output are given to the");
		std::println("                 tests.");
		std::println("                 This is repeated indefinitely, with care taken to minimize");
		std::println("                 the amount of duplicate seeds used.");
		//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
		std::println("  -ttep          Test target: Entropy pooling.  This should only be done on");
		std::println("                 RNGs that support entropy pooling.  It is similar to");
		std::println("                 -ttseed64, but the entropy accumulation methods are used");
		std::println("                 instead of simple seeding, and the amount of entropy used is");
		std::println("                 much larger.");
		std::println("  -walk_sequence  Some test-target modes will search seeds sequentially,");
		std::println("                  each subsequent seed 1 higher than the previous.");
		std::println("  -walk_greycode  Some test-target modes will search seeds in a simple");
		std::println("                  greycoded sequence, each subsequent seed at Hamming");
		std::println("                  distance 1 from the prior in a strict order.");
		std::println("  -walk_random    Some test-target modes will search seeds in a random walk,");
		std::println("                  each subsequent seed chosen at random from unused values");
		std::println("                  at Hamming distance 1 from the prior value.");
		std::println("  -walk_random_l  Some test-target modes will search seeds in a random walk,");
		std::println("                  each subsequent seed chosen at random from unused values");
		std::println("                  at Hamming distance 1 from the prior value, but lower bits");
		std::println("                  will be changed much more often than higher bits.");
		std::println("  -walk_random_h  Some test-target modes will search seeds in a random walk,");
		std::println("                  each subsequent seed chosen at random from unused values");
		std::println("                  at Hamming distance 1 from the prior value, but higher bits");
		std::println("                  will be changed much more often than lower bits.");
		std::println("test length options:");
		std::println("  -tlmin LENGTH  sets the minimum test length to LENGTH.  The tests will run on");
		std::println("                 that much data before it starts printing regular results.  A");
		std::println("                 large minimum will prevent it from displaying results on any");
		std::println("                 test lengths other than the maximum length (set by tlmax) and");
		std::println("                 lengths that were explicitly requested (by tlshow).");
		std::println("                 See notes on lengths for details on how to express the length");
		std::println("                 you want.");
		std::println("                 The default minimum is 1.5 seconds (-tlmin 1.5s).");
		std::println("  -tlmax LENGTH  sets the maximum test length to LENGTH.  The tests will stop");
		std::println("                 after that much data.  See notes on lengths for details on how");
		std::println("                 to express the length you want.");
		std::println("                 The default maximum is 32 tebibytes (-tlmax 32TB).");
		std::println("  -tlshow LENGTH sets an additional point at which to display interim results.");
		std::println("                 You can set multiple such points if desired.");
		std::println("                 These are in addition to the normal interim results points,");
		std::println("                 which are at every amount of data that is a power of 2 after");
		std::println("                 the minimum and before the maximum.");
		std::println("                 See the notes on lengths for details on how to express the");
		std::println("                 lengths you want.");
		std::println("  -tlfail        Halts testing after interim results are displayed if those");
		std::println("                 results include any failures. (default)");
		std::println("  -tlmaxonly     The opposite of -tlfail");
		std::println("other options:");
		std::println("  -multithreaded  enables multithreaded testing.  Typically up to 5 cores can");
		std::println("                  be used at once.");
		std::println("  -singlethreaded disables multithreaded testing.  (default)");
		std::println("  -seed SEED      specifies a 64 bit integer to seed the tested RNG with.  If");
		std::println("                  no seed is specified then a seed will be chosen randomly.");
		std::println("                  The value should be expressed in hexadecimal.  An '0x' prefix");
		std::println("                  on the seed is acceptable but not necessary.");
		std::println("notes on lengths:");
		//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
		std::println("  Each of the test length options requires a field named LENGTH.  These fields");
		std::println("  can accept either an amount of time or an amount of data.  In either case,");
		std::println("  several types of units are supported.");
		std::println("  A time should be expressed as a number postfixed with either s, m, h, or d,");
		std::println("  to express a number of seconds, minutes, hours, or days.");
		std::println("  example: -tlmin 1.4s (sets the minimum test length to 1.4 seconds)");
		std::println("  An amount of data can be expressed as a number with no postfix, in which case");
		std::println("  the number will be treated as the log-based-2 of the amount of bytes to test");
		std::println("  (in normal target mode) or the log-baed-2 of the number of seeds or strings");
		std::println("  to test in alternate test target modes.");
		std::println("  example: -tlmin 23 (sets the minimum test length to 8 mebibytes or 8");
		std::println("    million seeds, depending upon test target mode)");
		std::println("  Alternatively, an amount of data can be expressed as a number followed by");
		std::println("  KB, MB, GB, TB, or PB for kibibytes, mebibytes, gibibytes, tebibytes, or pebibytes.");
		std::println("  example: -tlmin 14KB (sets the minimum test length to 14 kibibytes");
		std::println("  If the B is omitted on KB, MB, GB, TB, or PB then it treats the metric");
		std::println("  prefixes as refering to numbers of bytes in normal test target mode, or");
		std::println("  numbers of seeds in seeding test target mode, or numbers of strings in");
		std::println("  entropy pooling test target mode.");
		std::println("  example: -tlmin 40M (sets the minimum test length to ~40 mebibytes or ~40");
		std::println("    million seeds, depending upon test target mode)");
		std::println("  A minor detail: I use binary prefixes (in which K means 1024 when");
		std::println("  dealing with quantities of binary information) not metric prefixes (in");
		std::println("  which K means 1000 no matter what is being dealt with, unless an 'i' follows");
		std::println("  the 'K').");
		//           12345678901234567890123456789012345678901234567890123456789012345678901234567890
		std::exit(0);
	}
	//PractRand::RNGs::vRNG *rng = get_rng(argv[1]);
	RNG_Factories::register_recommended_RNGs();
	RNG_Factories::register_nonrecommended_RNGs();
	RNG_Factories::register_input_RNGs();
	RNG_Factories::register_candidate_RNGs();
	RNG_Factories::RNG_factory_index["dummy"] = RNG_Factories::_generic_notrecommended_RNG_factory<DummyRNG>;
	Seeder_MetaRNG::register_name();
	EntropyPool_MetaRNG::register_name();
	std::string errmsg;
	RNGs::vRNG* rng = RNG_Factories::create_rng(argv[1], &errmsg);
	if (!rng) {
		if (errmsg.empty()) std::println(stderr, "unrecognized RNG name.  aborting.");
		else std::println(stderr, "{}", errmsg);
		std::exit(1);
	}

	bool do_self_test = true;
	bool use_multithreading = false;
	bool end_on_failure = true;
	bool smart_thresholds = true;
	double threshold = 0.1;
	int folding = 1;//0 = no folding, 1 = standard folding, 2 = extra folding
	int test_set_index = lookup_te_value(0);
	int mode = 0;//0 = normal, 1 = test seeding, 2 = test entropy pooling
	for (int i = 2; i < argc; i++) {
		int params_left = argc - i - 1;
		//-a
		//-e EXPECTED
		//-p THRESHOLD
		if constexpr (false) { ; }
		else if (!std::strcmp(argv[i], "-a")) {
			smart_thresholds = false;
			threshold = 1;
		}
		else if (!std::strcmp(argv[i], "-e")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			smart_thresholds = true;
			if (!parse_number(argv[++i], threshold) || threshold < 0.000001 || threshold > 1000) {
				std::println("invalid smart threshold: -e {} (must be between 0.000001 and 1000)", argv[i]);
				std::exit(0);
			}
		}
		else if (!std::strcmp(argv[i], "-p")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			smart_thresholds = false;
			if (!parse_number(argv[++i], threshold) || threshold < 0.0000000001 || threshold > 1.0) {
				std::println("invalid p-value threshold: -p {} (must be between 0.0000000001 and 1)", argv[i]);
				std::exit(0);
			}
		}
		//-tf FOLDING
		//-te EXPANDED
		else if (!std::strcmp(argv[i], "-tf")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			if (!parse_number(argv[++i], folding) || folding < 0 || folding > 2) {
				std::println("invalid folding test set value: -tf {}", argv[i]);
				std::exit(0);
			}
		}
		else if (!std::strcmp(argv[i], "-te")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			int expanded = 0;
			if (parse_number(argv[++i], expanded)) test_set_index = lookup_te_value(expanded);//0 maps to 0, but other values may not map to themselves
			else test_set_index = -1;
			if (test_set_index == -1) {
				std::println("invalid expanded test set value: -te {}", argv[i]);
				std::exit(0);
			}
		}
		//-ttnormal
		//-ttseed64
		//-ttep
		else if (!std::strcmp(argv[i], "-ttnormal")) { mode = 0; }
		else if (!std::strcmp(argv[i], "-ttseed64")) { mode = 1; }
		else if (!std::strcmp(argv[i], "-ttep"))     { mode = 2; }
		//-tlmin LENGTH
		//-tlmax LENGTH
		//-tlshow LENGTH
		else if (!std::strcmp(argv[i], "-tlmin")) { i++; }
		else if (!std::strcmp(argv[i], "-tlmax")) { i++; }
		else if (!std::strcmp(argv[i], "-tlshow")) { i++; }
		//-tlfail
		//-tlmaxonly
		else if (!std::strcmp(argv[i], "-tlfail")) { end_on_failure = true; }
		else if (!std::strcmp(argv[i], "-tlmaxonly")) { end_on_failure = false; }

		//-threads
		//-nothreads
		//-seed SEED
		else if (!std::strcmp(argv[i], "-multithreaded")) { use_multithreading = true; }
		else if (!std::strcmp(argv[i], "-singlethreaded")) { use_multithreading = false; }
		else if (!std::strcmp(argv[i], "-skip_selftest")) { do_self_test = false; }
		else if (!std::strcmp(argv[i], "-seed")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			seed_str = argv[++i];
		}
		else {
			std::println("unrecognized parameter: {}\naborting", argv[i]);
			std::exit(0);
		}
	}
#if !defined MULTITHREADING_SUPPORTED
	if (use_multithreading) {
		std::printf("multithreading is not supported on this build.  If multithreading should be supported, try defining the MULTITHREADING_SUPPORTED preprocessor symbol during the build process.\n");
		std::exit(0);
	}
#endif
	constexpr int TL_MIN = 1;
	constexpr int TL_MAX = 2;
	constexpr int TL_SHOW = 3;
	std::map<double,int> show_times;
	std::map<Uint64,int> show_datas;
	double show_min = -2.0;
	double show_max = 1ULL << 45;
	//walking parameters a second time to force the mode to be known prior to finding the test lengths
	for (int i = 2; i < argc; i++) {
		int params_left = argc - i - 1;
		if constexpr (false) { ; }
		else if (!std::strcmp(argv[i], "-tlmin")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			double length = interpret_length(argv[++i], !mode);
			if (!length) {std::println("invalid test length: {}", argv[i]);std::exit(0);}
			show_min = length;
		}
		else if (!std::strcmp(argv[i], "-tlmax")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			double length = interpret_length(argv[++i], !mode);
			if (!length) {std::println("invalid test length: {}", argv[i]);std::exit(0);}
			show_max = length;
		}
		else if (!std::strcmp(argv[i], "-tlshow")) {
			if (params_left < 1) {std::println("command line option {} must be followed by a value", argv[i]); std::exit(0);}
			double length = interpret_length(argv[++i], !mode);
			if (!length) {std::println("invalid test length: {}", argv[i]);std::exit(0);}
			if (length < 0) show_times[-length] = TL_SHOW;
			else show_datas[Uint64(length) / Tests::TestBlock::SIZE] = TL_SHOW;
		}
	}
	if (show_min < 0) show_times[-show_min] = TL_MIN;
	else show_datas[Uint64(show_min) / Tests::TestBlock::SIZE] = TL_MIN;
	if (show_max < 0) show_times[-show_max] = TL_MAX;
	else show_datas[Uint64(show_max) / Tests::TestBlock::SIZE] = TL_MAX;

	if (do_self_test) PractRand::self_test_PractRand();


	const auto start_time = std::chrono::steady_clock::now();

	Uint64 seed = known_good.raw32();//64 bit space, as that's what the interface accepts, but 32 bit random value so that by default it's not too onerous to record/compare/whatever the value by hand
	if (seed_str && !(rng->get_flags() & PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED)) {
		if (!interpret_seed(seed_str, seed)) {
			std::println("\"{}\" is not a valid 64 bit hexadecimal seed", seed_str);
			std::exit(0);
		}
	}
	known_good.seed(seed + 1);//the +1 is there just in case the RNG uses the same algorithm as the known good RNG

	PractRand::RNGs::vRNG* testing_rng = nullptr;
	if (mode == 0) {
		rng->seed(seed);
		testing_rng = rng;
	}
	else if (mode == 1) {
		//it would be nice to print a warning here for RNGs that use generic integer seeding
		//but that's a little difficult atm as there's no way to query whether an RNG does so
		testing_rng = new Seeder_MetaRNG(rng);
		testing_rng->seed(seed);
	}
	else if (mode == 2) {
		if (!(rng->get_flags() & PractRand::RNGs::FLAG::SUPPORTS_ENTROPY_ACCUMULATION)) {
			std::println("Entropy pooling is not supported by this RNG, so mode --ttep is invalid.");
			std::println("aborting");
			std::exit(0);
		}
		rng->reset_entropy();
		Uint64 a = rng->raw64();
		rng->reset_entropy();
		Uint64 b = rng->raw64();
		if (a != b) {
			std::println("entropy pooling RNG \"{}\" failed basic check 1.\naborting", rng->get_name());
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
			std::println("entropy pooling RNG \"{}\" failed basic check 2.\naborting", rng->get_name());
			std::exit(0);
		}
		if (c1 == e1 && c2 == e2) {
			std::println("entropy pooling RNG \"{}\" probably failed basic check 3.\naborting", rng->get_name());
			std::exit(0);
		}
		rng->seed(seed);
		//I'd like to test varying length entropy strings, but known good EPs are failing eventually when varying length is allowed for some reason
		testing_rng = new EntropyPool_MetaRNG(rng,48,64);
		testing_rng->seed(seed);
	}
	else {
		std::println("invalid mode, aborting");
		std::exit(1);
	}

	std::print("RNG = {}, seed = ", testing_rng->get_name());
	if (testing_rng->get_flags() & PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED) {
		if (seed_str) std::print("{}", seed_str);
		else std::print("unknown");
	}
	else {
		if (seed >> 32) std::print("0x{:x}{:08x}", long(seed >> 32), long((seed << 32) >> 32));
		else std::print("0x{:x}", long(seed));
	}
	const char* folding_names[3] = {"none", "standard", "extra"};
	std::print("\ntest set = {}, folding = {}", test_sets[test_set_index].name, folding_names[folding]);
	if (folding == 1) {
		int native_bits = testing_rng->get_native_output_size();
		if (native_bits > 0) std::print(" ({} bit)", native_bits);
		else std::print("(unknown format)");
	}

	std::println("\n");

	Tests::ListOfTests tests( static_cast<Tests::TestBaseclass*>(nullptr));
	if (test_set_index == -1) { std::println("internal error"); std::exit(1); }
	if constexpr (false) { ; }
	else if (folding == 0) { tests = test_sets[test_set_index].callback(); }
	else if (folding == 1) { tests = Tests::Batteries::apply_standard_foldings(testing_rng, test_sets[test_set_index].callback); }
	else if (folding == 2) { tests = Tests::Batteries::apply_extended_foldings(test_sets[test_set_index].callback); }
	else { std::println("internal error"); std::exit(1); }

//	Tests::ListOfTests tests = Tests::Batteries::get_expanded_standard_tests(rng);
#if defined MULTITHREADING_SUPPORTED
	TestManager* tman = nullptr;
	if (use_multithreading) tman = new MultithreadedTestManager(&tests, &known_good);
	else tman = new TestManager(&tests, &known_good);
#else
	TestManager* tman = new TestManager(&tests, &known_good);
#endif
	tman->reset(testing_rng);

	Uint64 blocks_tested = 0;
	bool already_shown = false;
	Uint64 next_power_of_2 = 1;
	bool showing_powers_of_2 = false;
	double time_passed = 0;
	while (true) {
		Uint64 blocks_to_test = next_power_of_2 - blocks_tested;
		constexpr int MAX_BLOCKS = 256 * 1024;
		if (blocks_to_test > MAX_BLOCKS) blocks_to_test = MAX_BLOCKS;
		while (!show_datas.empty()) {
			Uint64 data_checkpoint = show_datas.begin()->first - blocks_tested;
			if (data_checkpoint) {
				if (data_checkpoint < blocks_to_test) blocks_to_test = data_checkpoint;
				break;
			}
			int action = show_datas.begin()->second;
			if (action == TL_SHOW) {
				if (!already_shown) show_checkpoint(tman, mode, seed, time_passed, smart_thresholds, threshold, end_on_failure);
				already_shown = true;
			}
			else if (action == TL_MIN) {
				showing_powers_of_2 = true;
				if (!already_shown) show_checkpoint(tman, mode, seed, time_passed, smart_thresholds, threshold, end_on_failure);
				already_shown = true;
			}
			else if (action == TL_MAX) {
				if (!already_shown) show_checkpoint(tman, mode, seed, time_passed, smart_thresholds, threshold, end_on_failure);
				return 0;
			}
			else {std::println("internal error: unrecognized test length code, aborting");std::exit(1);}
			show_datas.erase(show_datas.begin());
		}
		while (!show_times.empty()) {
			double time_checkpoint = show_times.begin()->first - time_passed;
			if (time_checkpoint > 0) break;
			int action = show_times.begin()->second;
			show_times.erase(show_times.begin());
			if (action == TL_SHOW) {
				if (!already_shown) show_checkpoint(tman, mode, seed, time_passed, smart_thresholds, threshold, end_on_failure);
				already_shown = true;
			}
			else if (action == TL_MIN) { showing_powers_of_2 = true; }
			else if (action == TL_MAX) {
				if (!already_shown) show_checkpoint(tman, mode, seed, time_passed, smart_thresholds, threshold, end_on_failure);
				return 0;
			}
			else {std::println("internal error: unrecognized test length code, aborting");std::exit(1);}
		}

		if (blocks_tested == next_power_of_2) {
			if (showing_powers_of_2) {
				if (!already_shown) show_checkpoint(tman, mode, seed, time_passed, smart_thresholds, threshold, end_on_failure);
				already_shown = true;
			}
			next_power_of_2 <<= 1;
			continue;
		}
		tman->test(blocks_to_test);
		blocks_tested += blocks_to_test;
		already_shown = false;

		time_passed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
	}

	return 0;
}



