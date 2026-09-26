
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>



/*
Two functions are performed that may need to use platform-specific functionality:
	1. obtaining entropy for seeding
		preferably by asking the system cryptographic random number generator for the entropy (eg /dev/urandom)
		if that is not possible, then use the time (as measured with libc calls) instead
			to do: add a fallback case that blocks while listening to timing noise for entropy
	2. obtaining a unique 64 bit number that should never duplicate, and should be thread-safe.
		preferably with an atomic increment
		if that is not possible then with a malloc(1) that is never freed
			(this should get used at most once per thread creation)
*/

using namespace PractRand;

bool PractRand::Internals::add_entropy_automatically( PractRand::RNGs::vRNG* entropy_pool, [[maybe_unused]] int milliseconds ) {
	//the intention is for "millisecond" to be an amount of time that this function is permitted to spend on obtaining entropy
	//but currently nothing that spends time in a controlled fashion is implemented, so it's meaningless

	constexpr int DESIRED_BITS = 256;
	//constexpr int N32 = (DESIRED_BITS+31)/32;
	constexpr int N64 = (DESIRED_BITS+63)/64;

#if 1
	{//unix (linux/bsd/osx/etc, all flavors supposedly)
		//mostly safe to use even on platforms where it won't work
		std::FILE* f = std::fopen("/dev/urandom", "rb");
		uint64_t buf[N64];
		if (f) {
			if (std::fread(buf,N64*sizeof(buf[0]),1,f) == 1) {
				for (const auto i : buf) entropy_pool->add_entropy64(i);
				(void)std::fclose(f);
				entropy_pool->flush_buffers();
				std::memset(buf, 0, sizeof(buf));
				return true;
			}
			(void)std::fclose(f);
		}
	}
#endif
#if 0 //DISABLED
	{
	//disabled to avoid the possibility of blocking
		if (millseconds && f = std::fopen("/dev/random", "rb")) {
			//skip this if a good source was already found, because this can block
			uint64_t buf[N64];
			if(std::fread(buf,N64*sizeof(buf[0]),1,f)) {
				for (int i = 0; i < N64; i++) entropy_pool->add_entropy64(buf[i]);
				entropy_pool->flush_buffers();
				std::memset(buf, 0, sizeof(buf));
				return true;
			}
		}
	}
#endif
#if 1
	{//libc
		//not much entropy, but we take what we can get
		entropy_pool->add_entropy64(static_cast<uint64_t>(std::time(nullptr)));
		entropy_pool->add_entropy64(static_cast<uint64_t>(std::clock()));
		entropy_pool->add_entropy64(reinterpret_cast<uint64_t>(entropy_pool));
		auto* p = static_cast<uint64_t*>(std::malloc(sizeof(uint64_t)));
		entropy_pool->add_entropy64(reinterpret_cast<uint64_t>(p));
		//entropy_pool->add_entropy64(*p);//commented to avoid issues with memory debuggers
		free(p);
	}
#endif

	entropy_pool->flush_buffers();
	return false;
}

uint64_t PractRand::Internals::issue_unique_identifier ( ) {
#if 1
	static std::atomic<uint64_t> count = 0;
	return count++;
#elif defined __GNUC__
	static volatile uint64_t count = 0;
	return __sync_fetch_and_add(&count, uint64_t(1) );
#else
	//ugly, but without more knowledge of the target system or more dependencies there's not much more that can be done
	return (uint64_t)std::malloc(1);
#endif
}

/*
	//don't care about the units since it's only used as an entropy source
	//however, rdtscp is to be avoided since not enough CPUs support it
uint64_t PractRand::Internals::high_resolution_time() {
#if defined __GNUC__ && ( defined(__i386__) || defined(__x86_64__) )
	//from wikipedia
	uint32_t low, high;
	__asm__ __volatile__("rdtsc" : "=a"(low), "=d"(high) :: "ecx" );
	//__asm__ __volatile__("rdtscp" : "=a"(low), "=d"(high) :: "ecx" );
	return (uint64_t(high) << 32) | low;
#elif defined __GNUC__ && defined(__powerpc__)
	//from http://www.mcs.anl.gov/~kazutomo/rdtsc.html
	uint64_t result = 0;
	unsigned long int upper, lower,tmp;
	__asm__ volatile(
		"0:                  \n"
		"\tmftbu   %0           \n"
		"\tmftb    %1           \n"
		"\tmftbu   %2           \n"
		"\tcmpw    %2,%0        \n"
		"\tbne     0b         \n"
		: "=r"(upper),"=r"(lower),"=r"(tmp)
	);
	result = upper;
	result = result<<32;
	result = result|lower;

	return(result);
#else
	//to do: figure out the appropriate preprocessor defines to check for gettimeofday

	//very poor resolution, but we don't have a lot of alternatives at this point:
	return (uint64_t)std::clock();
#endif
}
*/
