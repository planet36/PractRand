#pragma once

#include <cstdint>

/*
Things to configure in this file:
1.  Integer types (compiler and CPU dependent)
2.  SIMD (compiler and CPU dependent, only used by ChaCha)
*/

/*
1.  Integer type sizes - reconfigure if needed
	If you know C/C++ programing this one should be self-explanatory.  The
default configuration should work for most 32 bit and 64 bit platforms.
If you are on a 16 bit or 8 bit platform then you will have to
reconfigure it.  PractRand does not use the C99 standard for this because
that is not yet standard for C++ (soon...), and some compilers I use do
not yet include the C99 integer size headers.
	If your compiler does not support 64 bit integers then PractRand can
not be built at all - many many things in PractRand require them.
	If this section is misconfigured then PractRand should fail to
compile with error messages about _compiletime_assert_integer_size.
*/
namespace PractRand {
	//typical for 32 and 64 bit compiler targets
	using Uint8 = std::uint8_t;
	using Uint16 = std::uint16_t;
	using Uint32 = std::uint32_t;
	using Uint64 = std::uint64_t;
	using Sint8 = std::int8_t;
	using Sint16 = std::int16_t;
	using Sint32 = std::int32_t;
	using Sint64 = std::int64_t;

	//double-checking those integer sizes:
	static_assert(sizeof(Uint8 )==1);
	static_assert(sizeof(Uint16)==2);
	static_assert(sizeof(Uint32)==4);
	static_assert(sizeof(Uint64)==8);
	static_assert(sizeof(Sint8 )==1);
	static_assert(sizeof(Sint16)==2);
	static_assert(sizeof(Sint32)==4);
	static_assert(sizeof(Sint64)==8);
}

/*
2. SIMD stuff can greatly speed up the ChaCha RNG
(not used by much yet, only ChaCha so far and only on MSVC)

SIMD strongly NOT recommended at this time, as I haven't come up with a
clean way to force a properly aligned malloc to be used when appropriate.
(in theory the spec may guarantee that all mallocs are aligned
appropriately, but in practice that does not seem to be the case)
I may have to have the Raw::ChaCha class contain pointers at memory blocks
instead of memory blocks themselves, there doesn't seem to be any other
way to allow a sane interface.
*/

#define PRACTRAND_NO_SIMD

#if defined _MSC_VER
#define PRACTRAND_ALIGN_128 __declspec(align(16))
#endif//defined _MSC_VER && _M_IX86_FP > 0
