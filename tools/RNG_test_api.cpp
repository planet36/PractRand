
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifndef MULTITHREADING_SUPPORTED
#define MULTITHREADING_SUPPORTED 1
#endif

// master header, includes everything in PractRand for both
//   practical usage and research...
//   EXCEPT it does not include specific algorithms
#include "PractRand_full.h"

// specific algorithms: all recommended RNGs
#include "PractRand/RNGs/all.h"

#include "SeedingTester.h"

#include "RNG_test_api.h"

using namespace PractRand;

#include "TestManager.h"
#ifdef MULTITHREADING_SUPPORTED
#include "MultithreadedTestManager.h"
#endif

typedef std::chrono::system_clock::rep TimeUnit;
TimeUnit get_time() {
    return std::chrono::system_clock::now().time_since_epoch().count();
}

double get_time_period() {
    return std::chrono::system_clock::period().num /
           static_cast<double>(std::chrono::system_clock::period().den);
}

class RNG_custom8 : public PractRand::RNGs::vRNG8 {
private:
    uint8_t (*raw8_callback)();
    void (*seed64_callback)(uint64_t seed);

public:
    RNG_custom8(uint8_t (*raw8_callback)())
        : raw8_callback(raw8_callback), seed64_callback(nullptr) {}

    RNG_custom8(uint8_t (*raw8_callback)(),
                void (*seed64_callback)(uint64_t seed))
        : raw8_callback(raw8_callback), seed64_callback(seed64_callback) {}

    PractRand::Uint8 raw8() override { return raw8_callback(); }
    PractRand::Uint64 get_flags() const override {
        return PractRand::RNGs::FLAG::STATE_UNAVAILABLE;
    }
    std::string get_name() const override { return "RNG_custom8"; }
    int get_native_output_size() const override { return -1; }
    void walk_state(PractRand::StateWalkingObject *) override {}

    void seed(Uint64 seed) override {
        seed64_callback ? seed64_callback(seed) : (void)seed;
    }
};

double print_result(const PractRand::TestResult &result,
                    bool print_header = false) {
    if (print_header)
        std::printf("  Test Name                         Raw       Processed   "
                    "  Evaluation\n");
    //                                     10        20        30        40 50
    //                                     60        70        80
    std::printf("  "); // 2 characters
    // NAME
    if constexpr (true) { // 34 characters
        std::printf("%s", result.name.c_str());
        int len = result.name.length();
        for (int i = len; i < 34; i++)
            std::printf(" ");
    }

    // RAW TEST RESULT
    if constexpr (true) { // 10 characters?
        double raw = result.get_raw();
        if (raw > 99999.0)
            std::printf("R>+99999  ");
        else if (raw < -99999.0)
            std::printf("R<-99999  ");
        else if (std::abs(raw) < 999.95)
            std::printf("R=%+6.1f  ", raw);
        else
            std::printf("R=%+6.0f  ", raw);
        // if (std::abs(raw) < 99999.5) std::printf(" ");
        // if (std::abs(raw) < 999999.5) std::printf(" ");
        // if (std::abs(raw) < 9999999.5) std::printf(" ");
    }

    // RESULT AS A NUMERICAL "SUSPICION LEVEL" (log of distance from pvalue to
    // closest extrema)
    if constexpr (false) { // 12 characters?
        bool printed = false;
        double susp = result.get_suspicion();
        if (result.type == result.TYPE_PASSFAIL)
            std::printf("  %s    ",
                        result.get_pvalue() ? "\"pass\"" : "\"fail\"");
        else if (result.type == result.TYPE_RAW)
            std::printf("            ");
        else if (result.type == result.TYPE_BAD_P ||
                 result.type == result.TYPE_BAD_S ||
                 result.type == result.TYPE_RAW_NORMAL) {
            std::printf("S=~%+6.1f   ", susp);
            printed = true;
        } else {
            std::printf("S =%+6.1f   ", susp);
            printed = true;
        }
        if (printed) {
            if (std::abs(susp) < 9999.95)
                std::printf(" ");
            if (std::abs(susp) < 999.95)
                std::printf(" ");
        }
    }

    // RESULT AS A p-value
    if constexpr (true) { // 14 characters?
        if (result.type == result.TYPE_PASSFAIL)
            std::printf("  %s      ",
                        result.get_pvalue() ? "\"pass\"" : "\"fail\"");
        else if (result.type == result.TYPE_RAW)
            std::printf("              ");
        else if (result.type == result.TYPE_BAD_P ||
                 result.type == result.TYPE_GOOD_P ||
                 result.type == result.TYPE_RAW_NORMAL) {
            double p = result.get_pvalue();
            double a = std::abs(p - 0.5);
            std::printf((result.type != result.TYPE_GOOD_P) ? "p~= " : "p = ");
            if (a > 0.49) {
                double s = result.get_suspicion();
                double ns = std::abs(s) + 1;
                double dec = ns / (std::log(10.0) / std::log(2.0));
                double dig = std::ceil(dec);
                double sig = std::floor(std::pow(0.1, dec - dig));
                if (dig > 999)
                    std::printf(" %d        ", (s > 0) ? 1 : 0);
                else {
                    if (s > 0)
                        std::printf("1-%1.0fe-%.0f  ", sig, dig);
                    else
                        std::printf("  %1.0fe-%.0f  ", sig, dig);
                    if (dig < 100)
                        std::printf(" ");
                    if (dig < 10)
                        std::printf(" ");
                }
            } else if (result.type == result.TYPE_GOOD_P)
                std::printf("%5.3f     ", p);
            else if (a >= 0.4)
                std::printf("%4.2f      ", p);
            else
                std::printf("%3.1f       ", p);
        } else if (result.type == result.TYPE_BAD_S ||
                   result.type == result.TYPE_GOOD_S) {
            double s = result.get_suspicion();
            double p = result.get_pvalue();
            std::printf((result.type == result.TYPE_BAD_S ||
                         result.type == result.TYPE_RAW_NORMAL)
                            ? "p~="
                            : "p =");
            if (p >= 0.01 && p <= 0.99)
                std::printf(" %.3f     ", p);
            else {
                double ns = std::abs(s) + 1;
                double dec = ns / (std::log(10.0) / std::log(2.0));
                double dig = std::ceil(dec);
                double sig = std::pow(0.1, dec - dig);
                sig = std::floor(sig * 10) * 0.1;
                if (dig > 9999)
                    std::printf(" %d         ", (s > 0) ? 1 : 0);
                else if (dig > 999) {
                    sig = std::floor(sig);
                    if (s > 0)
                        std::printf("1-%1.0fe-%.0f  ", sig, dig);
                    else
                        std::printf("  %1.0fe-%.0f  ", sig, dig);
                } else {
                    if (s > 0)
                        std::printf("1-%3.1fe-%.0f ", sig, dig);
                    else
                        std::printf("  %3.1fe-%.0f ", sig, dig);
                    if (dig < 100)
                        std::printf(" ");
                    if (dig < 10)
                        std::printf(" ");
                }
            }
        }
    }

    double dec = std::log(2.) / std::log(10.0);
    double as = (std::abs(result.get_suspicion()) + 1 - 1) *
                dec; // +1 for suspicion conversion, -1 to account for there
                     // being 2 failure regions (near-zero and near-1)
    double wmod = std::log(result.get_weight()) / std::log(0.5) * dec;
    double rs = as - wmod;
    // double ap = std::abs(0.5 - result.get_pvalue());
    // MESSAGE DESCRIBING RESULT IN ENGLISH
    if constexpr (true) { // 17 characters?
        /*
            Threshold Values:
            The idea is to assign a suspicioun level based not just upon the
            p-value but also the number of p-values and their relative
           importance. If there are a million p-values then we probably don't
           care about anything less extreme than a one in ten million event. But
           if there's one important p-value and a million unimportant ones then
            the important one doesn't have to be that extreme to rouse our
           suspicion.

            Output Format:
            unambiguous failures are indented 2 spaces to make them easier to
           spot probable failures with (barely) enough room for ambiguity are
           indented 1 space the most extreme failures get a sequence of
           exclamation marks to distinguish them
        */
        if constexpr (false)
            ;
        else if (rs > 999)
            std::printf("  FAIL !!!!!!!!  ");
        else if (rs > 325)
            std::printf("  FAIL !!!!!!!   ");
        else if (rs > 165)
            std::printf("  FAIL !!!!!!    ");
        else if (rs > 85)
            std::printf("  FAIL !!!!!     ");
        else if (rs > 45)
            std::printf("  FAIL !!!!      ");
        else if (rs > 25)
            std::printf("  FAIL !!!       ");
        else if (rs > 17)
            std::printf("  FAIL !!        ");
        else if (rs > 12)
            std::printf("  FAIL !         ");
        else if (rs > 8.5)
            std::printf("  FAIL           ");
        else if (rs > 6.0)
            std::printf(" VERY SUSPICIOUS ");
        else if (rs > 4.0)
            std::printf("very suspicious  ");
        else if (rs > 3.0)
            std::printf("suspicious       ");
        else if (rs > 2.0)
            std::printf("mildly suspicious");
        else if (rs > 1.0)
            std::printf("unusual          ");
        else if (rs > 0.0)
            std::printf("normalish       ");
        else
            std::printf("normal           ");
    }
    std::printf("\n");
    return rs;
}

extern "C" {
double show_checkpoint(TestManager *tman, int mode, uint64_t seed, double time,
                       bool smart_thresholds, double threshold) {

    std::printf("length= ");
    Uint64 length = tman->get_blocks_so_far() * Tests::TestBlock::SIZE;
    double log2b = std::log(double(length)) / std::log(2.0);
    const char *unitstr[6] = {"kilobyte", "megabyte", "gigabyte",
                              "terabyte", "petabyte", "exabyte"};
    int units = int(std::floor(log2b / 10)) - 1;
    if (units < 0 || units > 5) {
        std::printf("internal error: units=%d length out of bounds?\n", units);
        std::exit(1);
    }
    if (length & (length - 1))
        std::printf("%.3f %ss", length * std::pow(0.5, units * 10.0 + 10),
                    unitstr[units]);
    else
        std::printf("%.0f %s%s", length * std::pow(0.5, units * 10.0 + 10),
                    unitstr[units],
                    length != (Uint64(1024) << (units * 10)) ? "s" : "");
    if (length & (length - 1))
        std::printf(" (2^%.3f", log2b - (mode ? 3 : 0));
    else
        std::printf(" (2^%.0f", log2b - (mode ? 3 : 0));
    const char *mode_unit_names[3] = {"bytes", "seeds", "entropy strings"};
    std::printf(" %s), time= ", mode_unit_names[mode]);
    if (time < 99.95)
        std::printf("%.1f seconds\n", time);
    else
        std::printf("%.0f seconds\n", time);

    std::vector<PractRand::TestResult> results;
    tman->get_results(results);
    double total_weight = 0, min_weight = 9999999;
    for (unsigned int i = 0; i < results.size(); i++) {
        double weight = results[i].get_weight();
        total_weight += weight;
        if (weight < min_weight)
            min_weight = weight;
    }
    if (min_weight <= 0) {
        std::printf("error: result weight too small\n");
        std::exit(1);
    }
    std::vector<int> marked;
    for (unsigned int i = 0; i < results.size(); i++) {
        results[i].set_weight(results[i].get_weight() / total_weight);
        if (!smart_thresholds) {
            if (std::abs(0.5 - results[i].get_pvalue()) < 0.5 - threshold)
                continue;
        } else {
            double T = threshold * results[i].get_weight() * 0.5;
            if (std::abs(0.5 - results[i].get_pvalue()) < (0.5 - T))
                continue;
        }
        marked.push_back(i);
    }
    double biggest_decimal_suspicion = 0;
    for (unsigned int i = 0; i < marked.size(); i++) {
        double decimal_suspicion = print_result(results[marked[i]], i == 0);
        if (decimal_suspicion > biggest_decimal_suspicion)
            biggest_decimal_suspicion = decimal_suspicion;
    }
    if (marked.size() == results.size())
        ;
    else if (marked.size() == 0)
        std::printf("  no anomalies in %d test result(s)\n",
                    int(results.size()));
    else
        std::printf("  ...and %d test result(s) without anomalies\n",
                    int(results.size() - marked.size()));
    std::printf("\n");
    std::fflush(stdout);

    return biggest_decimal_suspicion;
}
} // extern "C"

double interpret_length(const std::string &lengthstr, bool normal_mode) {
    //(0-9)*[.(0-9)*][((K|M|G|T|P)[B])|(s|m|h|d)]
    int mode_factor = normal_mode ? 1 : 8;
    unsigned int pos = 0;
    double value = 0;
    for (; pos < lengthstr.size(); pos++) {
        char c = lengthstr[pos];
        if (c < '0')
            break;
        if (c > '9')
            break;
        value = value * 10 + (c - '0');
    }
    if (!pos)
        return 0;
    if (pos == lengthstr.size())
        return std::pow(2.0, value) * mode_factor;
    if (lengthstr[pos] == '.') {
        pos++;
        double sig = 0.1;
        for (; pos < lengthstr.size(); pos++, sig *= 0.1) {
            char c = lengthstr[pos];
            if (c < '0')
                break;
            if (c > '9')
                break;
            value += (c - '0') * sig;
        }
        if (pos == lengthstr.size())
            return std::pow(2.0, value) * mode_factor;
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
        scale = -1; // one second
        expect_B = false;
        break;
    case 'm':
        scale = -60; // one minute
        expect_B = false;
        break;
    case 'h':
        scale = -3600; // one hour
        expect_B = false;
        break;
    case 'd':
        scale = -86400; // one day
        expect_B = false;
        break;
    default:
        break;
    }
    pos++;
    if (pos == lengthstr.size()) {
        if (scale < 0) {
            if (value < 0.05)
                value = 0.05;
            return value * scale;
        }
        return value * scale * mode_factor;
    }
    if (!expect_B)
        return 0;
    if (lengthstr[pos++] != 'B')
        return 0;
    if (pos != lengthstr.size())
        return 0;
    return value * scale;
}
bool interpret_seed(const std::string &seedstr, Uint64 &seed) {
    // would prefer strtol, but that is insufficiently portable when it has to
    // handle 64 bit values
    Uint64 value = 0;
    unsigned int position = 0;
    if (seedstr.length() >= 3 && seedstr[0] == '0' && seedstr[1] == 'x')
        position = 2;
    while (position < seedstr.length()) {
        int c = seedstr[position++];
        if (value >> 60)
            return false; // too long
        value *= 16;
        if (c >= '0' && c <= '9')
            value += (c - '0');
        else if (c >= 'a' && c <= 'f')
            value += (c - 'a') + 10;
        else if (c >= 'A' && c <= 'F')
            value += (c - 'A') + 10;
        else
            return false; // invalid character
    }
    seed = value;
    return true;
}

#include "PractRand/Tests/Birthday.h"
#include "PractRand/Tests/DistFreq4.h"
#include "PractRand/Tests/FPF.h"
#include "PractRand/Tests/FPMulti.h"
#include "PractRand/Tests/Gap16.h"

static uint64_t state = 5342;
static uint64_t inc = 521;

uint32_t pcg32() {
    uint64_t oldstate = state;
    // Advance internal state
    state = oldstate * 6364136223846793005ULL + (inc | 1);
    // Calculate output function (XSH RR), uses old state for
    // maximum instruction level parallelism
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

uint8_t raw8_callback() { return pcg32(); }
// uint8_t raw8_callback() { return 0; }

extern "C" {

int main() { run_tests(raw8_callback, 21, true, true, true); }

int run_tests(uint8_t (*raw8_callback)(), uint32_t log2_max_test_blocks,
              bool smart_thresholds, bool end_on_failure,
              bool use_multithreading) {
    PractRand::initialize_PractRand();
    PractRand::hook_error_handler(PractRand::print_err);

    uint64_t max_test_blocks = 1ull << log2_max_test_blocks;

    std::string errmsg;
    RNGs::vRNG *testing_rng = new RNG_custom8(raw8_callback);

    double threshold = 0.1;
    int folding = 1; // 0 = no folding, 1 = standard folding, 2 = extra folding
    int mode = 0;    // 0 = normal, 1 = test seeding, 2 = test entropy pooling

    enum {
        TL_MIN = 1,
        TL_MAX = 2,
        TL_SHOW = 3
    }; //, TL_FLAG_UNSPECIFIED_UNITS = 16};
    std::map<double, int> show_times;
    std::map<Uint64, int> show_datas;
    double show_min = -2.0;
    double show_max = 1ull << 45;

    if (show_min < 0)
        show_times[-show_min] = TL_MIN;
    else
        show_datas[Uint64(show_min) / Tests::TestBlock::SIZE] = TL_MIN;
    if (show_max < 0)
        show_times[-show_max] = TL_MAX;
    else
        show_datas[Uint64(show_max) / Tests::TestBlock::SIZE] = TL_MAX;

    // PractRand::self_test_PractRand();

    std::time_t start_time = std::time(nullptr);
    TimeUnit start_clock = get_time();

    PractRand::RNGs::Polymorphic::hc256 known_good(PractRand::SEED_AUTO);

    Uint64 seed =
        known_good
            .raw32(); // 64 bit space, as that's what the interface accepts, but
                      // 32 bit random value so that by default it's not too
                      // onerous to record/compare/whatever the value by hand
    known_good.seed(seed + 1); // the +1 is there just in case the RNG uses the
                               // same algorithm as the known good RNG

    //(mode == 0)
    // testing_rng->seed(seed);

    Tests::ListOfTests tests = Tests::Batteries::get_core_tests();

    printf("running core tests\n");

#if defined MULTITHREADING_SUPPORTED
    printf("Using multithreading\n");
    TestManager *tman;
    if (use_multithreading)
        tman = new MultithreadedTestManager(&tests, &known_good);
    else
        tman = new TestManager(&tests, &known_good);
#else
    TestManager *tman = new TestManager(&tests, &known_good);
#endif
    tman->reset(testing_rng);

    Uint64 blocks_tested = 0;
    bool already_shown = true;
    Uint64 next_power_of_2 = 1;
    bool showing_powers_of_2 = false;
    double time_passed = 0;

    double BDS_THRESHOLD = 8.5;
    uint64_t MAX_BLOCKS = 256 * 1024;
    double bds = 0; //  biggest decimal suspicion
    while (true) {
        Uint64 blocks_to_test = next_power_of_2 - blocks_tested;
        if (blocks_to_test > MAX_BLOCKS)
            blocks_to_test = MAX_BLOCKS;
        while (!show_datas.empty()) {
            Uint64 data_checkpoint = show_datas.begin()->first - blocks_tested;
            if (data_checkpoint) {
                if (data_checkpoint < blocks_to_test)
                    blocks_to_test = data_checkpoint;
                break;
            }
            int action = show_datas.begin()->second;

            if (action == TL_SHOW) {
            } else if (action == TL_MIN) {
                showing_powers_of_2 = true;
            } else if (action == TL_MAX) {
                printf("action == TL_MAX\n");
                return 0;
            } else {
                std::printf("internal error: unrecognized test length code, "
                            "aborting\n");
                return -1;
            }

            if (!already_shown) {
                bds = show_checkpoint(tman, mode, seed, time_passed,
                                      smart_thresholds, threshold);
                already_shown = true;
            }
            if (end_on_failure && bds > BDS_THRESHOLD) {
                std::printf("biggest decimal suspicion = %f > %f, aborting\n",
                            bds, BDS_THRESHOLD);
                return -1;
            }
            show_datas.erase(show_datas.begin());
        }
        while (!show_times.empty()) {
            double time_checkpoint = show_times.begin()->first - time_passed;
            if (time_checkpoint > 0)
                break;
            int action = show_times.begin()->second;
            show_times.erase(show_times.begin());
            if (action == TL_SHOW) {
            } else if (action == TL_MIN) {
                showing_powers_of_2 = true;
                // already_shown = true;
            } else if (action == TL_MAX) {
                return 0;
            } else {
                std::printf("internal error: unrecognized test length code, "
                            "aborting\n");
                return -1;
            }
            if (!already_shown) {
                bds = show_checkpoint(tman, mode, seed, time_passed,
                                      smart_thresholds, threshold);
                already_shown = true;
            }

            if (end_on_failure && bds > BDS_THRESHOLD) {
                std::printf(
                    "xx biggest decimal suspicion = %f > %f, aborting\n", bds,
                    BDS_THRESHOLD);
                return -2;
            }
        }

        if (!showing_powers_of_2 && (blocks_tested == next_power_of_2)) {
            already_shown = true;
        }
        if (!already_shown) {
            bds = show_checkpoint(tman, mode, seed, time_passed,
                                  smart_thresholds, threshold);

            if (end_on_failure && bds > BDS_THRESHOLD) {
                std::printf(
                    "xx biggest decimal suspicion = %f > %f, aborting\n", bds,
                    BDS_THRESHOLD);
                return -2;
            }
            already_shown = true;
        }
        if (blocks_tested > max_test_blocks)
            break;

        if (blocks_tested == next_power_of_2) {
            next_power_of_2 <<= 1;
            continue;
        }
        tman->test(blocks_to_test);
        blocks_tested += blocks_to_test;
        already_shown = false;

        double clocks_passed = (get_time() - start_clock) *
                               get_time_period(); // may wrap too quickly
        int seconds_passed = std::time(nullptr) - start_time;
        if (seconds_passed >= 1000 || seconds_passed > clocks_passed + 2.0)
            time_passed = seconds_passed;
        else
            time_passed = clocks_passed;
    }
    printf("done\n");

    return 0;
}
} // extern "C"
