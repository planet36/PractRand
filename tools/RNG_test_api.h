#ifndef _RNG_TEST_API_H_
#define _RNG_TEST_API_H_

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#else
#ifndef bool
#define bool uint8_t
#endif // ifndef bool
#endif // __cplusplus

int run_tests(uint8_t (*raw8_callback)(), uint32_t log2_max_test_blocks,
              bool smart_thresholds, bool end_on_failure,
              bool use_multithreading);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // _RNG_TEST_API_H_