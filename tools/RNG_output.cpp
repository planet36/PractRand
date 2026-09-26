#include <charconv>
#include <cmath>
#include <csignal>     /* signal, sig_atomic_t */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <print>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#ifdef _WIN32 // needed to allow binary stdout on windows
#include <fcntl.h>
#include <io.h>
#endif

//master header, includes everything in PractRand for both
//  practical usage and research...
//  EXCEPT it does not include specific algorithms
//  also it does not include PractRand/RNG_adaptors.h as that includes lots of templated stuff
#include "PractRand_full.h"
//the full version is needed because non-recommended RNGs are supported

//specific RNG algorithms, to produce (pseudo-)random numbers
#include "PractRand/RNGs/all.h"

#include "PractRand/RNGs/other/fibonacci.h"
#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/other/mult.h"
#include "PractRand/RNGs/other/simple.h"
#include "PractRand/RNGs/other/special.h"
#include "PractRand/RNGs/other/transform.h"

//not actually part of the library headers, just some inline code for an abstract factory for PractRand RNG name -> in
#include "RNG_from_name.h"
#include "parse_number.h"

using namespace PractRand;
#include "Candidate_RNGs.h"


bool interpret_seed(const std::string& seedstr, Uint64& seed) {
	const char* first = seedstr.data();
	const char* last = first + seedstr.size();
	if (seedstr.starts_with("0x")) first += 2;
	auto [ptr, ec] = std::from_chars(first, last, seed, 16);
	return ec == std::errc() && ptr == last;
}
void print_usage(const char* program_name) {
	std::cerr << "usage:\n\t" << program_name << " RNG_name bytes_to_output [64bit_hexadecimal_seed]\n";
	std::cerr << "  example:\n\t" << program_name << " jsf32 16\n";
	std::cerr << "  prints 16 bytes using the jsf32 RNG with a randomly chosen seed, with an \n";
	std::cerr << "    error message if fewer than 16 bytes were successfully outputted.\n";
	std::cerr << "usage:\n\t" << program_name << " RNG_name inf [64bit_hexadecimal_seed]\n";
	std::cerr << "  as above, but it prints indefinitely, with no error message when aborted\n";
	std::cerr << "usage:\n\t" << program_name << " RNG_name name\n";
	std::cerr << "  It prints the result of (RNG)->get_name()\n";
	std::cerr << "  Which is often the same as the RNG_name parameter, but not always.\n";
	exit(0);
}

sig_atomic_t signaled = 0;

void signal_handler(int param)
{
	signaled = param;
}

#include "SeedingTester.h"

int main(int argc, char** argv) { // NOLINT(bugprone-exception-escape)
#ifdef _WIN32
	_setmode( _fileno(stdout), _O_BINARY); // needed to allow binary stdout on windows
#endif
	if (argc < 3 || argc > 4) print_usage(argv[0]);
	PractRand::initialize_PractRand();
	PractRand::hook_error_handler(PractRand::print_err);

	RNG_Factories::register_recommended_RNGs();
	RNG_Factories::register_nonrecommended_RNGs();
	RNG_Factories::register_candidate_RNGs();
	Seeder_MetaRNG::register_name();
	EntropyPool_MetaRNG::register_name();
	std::string errmsg;
	RNGs::vRNG* rng = RNG_Factories::create_rng(argv[1], &errmsg);

	if (!rng) {
		if (errmsg.empty()) { std::println(stderr, "RNG_output ERROR: unrecognized RNG name"); print_usage(argv[0]); }
		else { std::println(stderr, "RNG_output ERROR: RNG_Factories returned error message:\n{}", errmsg); exit(1); }
	}

	double _n = 0;//stays 0 if argv[2] is not a number, such as "name"
	parse_number(argv[2], _n);
	Uint64 n = 0;
	if (_n <= 0 || _n >= 18446744073709551616.0) {
		if (!strcmp(argv[2], "name")) {
			std::println("{}", rng->get_name());
			exit(0);
		}
		else if (!strcmp(argv[2], "inf")) {
			_n = 0;
			n = 0xFFFFffffFFFFffffULL;
		}
		else {
			std::println(stderr, "RNG_output ERROR: invalid number of output bytes"); print_usage(argv[0]);
		}
	}
	else { n = Uint64(_n); }

	if (argc == 3) { rng->autoseed(); }
	else {
		Uint64 seed = 0;
		if (!interpret_seed(argv[3],seed)) {std::println(stderr, "RNG_output ERROR: \"{}\" is not a valid 64 bit hexadecimal seed", argv[3]); std::exit(0);}
		rng->seed(seed);
	}

	void(*prev_handler)(int) = nullptr;
	prev_handler = signal(SIGINT, signal_handler);  if (prev_handler == SIG_ERR) { std::cerr << "WARNING: Setting signal handler for SIGINT has failed.\n"; }
	prev_handler = signal(SIGTERM, signal_handler); if (prev_handler == SIG_ERR) { std::cerr << "WARNING: Setting signal handler for SIGTERM has failed.\n"; }
#ifdef __linux__
	prev_handler = signal(SIGPIPE, signal_handler); if (prev_handler == SIG_ERR) { std::cerr << "WARNING: Setting signal handler for SIGPIPE has failed.\n"; }
#endif

	constexpr int BUFFER_SIZE = 8;
	//Uint64 buffer[BUFFER_SIZE];
	PractRand::Tests::TestBlock buffer[BUFFER_SIZE];
	while (n && !signaled) {
		//for (int i = 0; i < BUFFER_SIZE; i++) buffer[i] = rng->raw64();
		std::size_t bytes_this_loop = n > (BUFFER_SIZE * PractRand::Tests::TestBlock::SIZE) ? (BUFFER_SIZE * PractRand::Tests::TestBlock::SIZE) : n;
		buffer[0].fill(rng, (bytes_this_loop + PractRand::Tests::TestBlock::SIZE - 1) >> PractRand::Tests::TestBlock::SIZE_L2);
		size_t bytes_written = std::fwrite(&buffer[0], 1, bytes_this_loop, stdout);
		n -= bytes_written;
		if ( bytes_written != bytes_this_loop ) {
			//if (std::ferror(stdout)) std::perror("I/O error when writing to standard output"); // this was generating spurious error messages on windows
			break;
		}
	}
	(void)std::fflush(stdout);
	if (signaled) {
		//std::cerr << "WARNING: Received signal " << signaled << ". Closing the application." << std::endl; // this was generating spurious error messages on linux
	}
	if (n && _n) {
		std::cerr << "RNG_output ERROR: " << Uint64(_n) << " bytes were requested, but only " << (Uint64(_n) - n) << " bytes were written.\n";
	}
	return 0;
}
