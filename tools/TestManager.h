#pragma once

#include "PractRand/RNGs/hc256.h"
#include "PractRand/rng_basics.h"
#include "PractRand/test_batteries.h"
#include "PractRand/tests.h"

#include <cstring>
#include <vector>

class TestManager {
protected:
	std::vector<PractRand::Tests::TestBlock> buffer;
	PractRand::RNGs::vRNG* rng{nullptr};
	PractRand::RNGs::vRNG* known_good;
	PractRand::Tests::ListOfTests* tests;
	unsigned int max_buffer_amount;
	int prefix_blocks;
	int main_blocks;
	int blocks_to_repeat{};
	uint64_t blocks_so_far;
	bool freshly_created;
	int prep_blocks(uint64_t& blocks);
public:
	[[nodiscard]] const PractRand::RNGs::vRNG* get_rng() const {return rng;}//RNG being tested
	[[nodiscard]] uint64_t get_blocks_so_far() const {return blocks_so_far;}//number of blocks tested

	explicit TestManager(PractRand::Tests::ListOfTests* tests_, PractRand::RNGs::vRNG* known_good_=nullptr, int max_buffer_amount_ = 1 << (25-10));
	//rng_ = RNG to test
	//tests_ = list of tests use on the RNG
	//known_good_ = sometimes the tests or test manager need good random numbers for some reason
	//max_buffer_amount_ = size in kibibytes of the maximum amount of random data to keep buffered up at on time

	virtual ~TestManager();//destructor (destroys the tests in the ListOfTests)

	virtual void reset(PractRand::RNGs::vRNG* rng_);//resets contents for starting a new test run ; if rng is NULL then it will reuse the current RNG

	virtual void test(uint64_t blocks);//does testing... the number of blocks is ADDITIONAL blocks to test, not total blocks to test

	virtual void finish_testing() {}//waits until every block passed to test() is tested

	virtual void get_results( std::vector<PractRand::TestResult>& result_vec ) final;//gets the results
};

TestManager::TestManager(PractRand::Tests::ListOfTests* tests_, PractRand::RNGs::vRNG* known_good_, int max_buffer_amount_) : known_good(known_good_), tests(tests_) {
	if (!known_good) known_good = new PractRand::RNGs::Polymorphic::hc256(PractRand::SEED_AUTO);
	blocks_so_far = 0;
	max_buffer_amount = max_buffer_amount_;
	prefix_blocks = 0;
	main_blocks = 0;
	for (auto& test : tests->tests) test->init(known_good);
	freshly_created = true;
}
TestManager::~TestManager() {
	for (auto& test : tests->tests) test->deinit();
	for (auto& test : tests->tests) delete test;
}
void TestManager::reset(PractRand::RNGs::vRNG* rng_) {
	if (!freshly_created) for (auto& test : tests->tests) test->deinit();
	freshly_created = false;
	for (auto& test : tests->tests) test->init(known_good);
	blocks_to_repeat = 0;
	for (auto& test : tests->tests) {
		int rb = test->get_blocks_to_repeat();
		if (blocks_to_repeat < rb) blocks_to_repeat = rb;
	}
	buffer.resize(max_buffer_amount + blocks_to_repeat);
	if (rng_) rng = rng_;
	main_blocks = 0;
	prefix_blocks = 0;
	blocks_so_far = 0;
}
int TestManager::prep_blocks(uint64_t& blocks) {
	uint64_t _delta_blocks = blocks;
	if (_delta_blocks > max_buffer_amount) _delta_blocks = max_buffer_amount;
	int delta_blocks = int(_delta_blocks);
	blocks -= delta_blocks;
	size_t repeat_region_start = 0, repeat_region_size = 0;
	if (prefix_blocks + main_blocks >= blocks_to_repeat) {
		repeat_region_start = prefix_blocks + main_blocks - blocks_to_repeat;
		repeat_region_size = blocks_to_repeat;
	}
	else {
		repeat_region_start = 0;
		repeat_region_size = prefix_blocks + main_blocks;
	}
	if (repeat_region_start != 0)
		std::memmove(buffer.data(), &buffer[repeat_region_start], repeat_region_size * PractRand::Tests::TestBlock::SIZE);
	prefix_blocks = repeat_region_size;
	main_blocks = delta_blocks;
	buffer[prefix_blocks].fill(rng, main_blocks);
	blocks_so_far += delta_blocks;
	return delta_blocks;
}
void TestManager::test(uint64_t num_blocks) {
	while (num_blocks) {
		prep_blocks(num_blocks);
		for (auto& test : tests->tests)
			test->test_blocks(&buffer[prefix_blocks], main_blocks);
	}
}
void TestManager::get_results( std::vector<PractRand::TestResult>& result_vec ) {
	finish_testing();
	for (auto& test : tests->tests) {
		test->get_results(result_vec);
	}
}
