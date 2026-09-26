#pragma once

#include <thread>

class MultithreadedTestManager : public TestManager {
	static void run_test(PractRand::Tests::TestBaseclass* test, PractRand::Tests::TestBlock* base_block, Uint64 numblocks) {
		constexpr int MAX_BLOCKS_PER_CALL = 1ULL << 18;
		while (numblocks > MAX_BLOCKS_PER_CALL) {
			test->test_blocks(base_block,MAX_BLOCKS_PER_CALL);
			numblocks -= MAX_BLOCKS_PER_CALL;
			base_block += MAX_BLOCKS_PER_CALL;
		}
		if (numblocks) test->test_blocks(base_block,static_cast<int>(numblocks));
	}
	std::vector<std::jthread> threads;
	void wait_on_threads() {
		for (auto& thread : threads) thread.join();
		threads.clear();
	}


public:
	//std::vector<PractRand::Tests::TestBlock> buffer1, buffer2;
	std::vector<PractRand::Tests::TestBlock> alt_buffer;
	//RNGs::vRNG *rng;
	//RNGs::vRNG *known_good;
	//Tests::ListOfTests *tests;
	//unsigned int max_buffer_amount;
	//int prefix_blocks;
	//int main_blocks;
	//Uint64 blocks_so_far;
	void multithreaded_prep_blocks(Uint64 num_blocks) {
		int new_prefix_blocks = blocks_to_repeat;
		if (new_prefix_blocks > prefix_blocks + main_blocks)
			new_prefix_blocks = prefix_blocks + main_blocks;
		if (new_prefix_blocks) {
			std::memcpy(
				buffer.data(),
				&alt_buffer[prefix_blocks + main_blocks - new_prefix_blocks],
				PractRand::Tests::TestBlock::SIZE * new_prefix_blocks
			);
		}
		prefix_blocks = new_prefix_blocks;
		main_blocks = (num_blocks > max_buffer_amount) ? max_buffer_amount : Uint32(num_blocks);
		buffer[prefix_blocks].fill(rng, main_blocks);
		blocks_so_far += main_blocks;
	}

	MultithreadedTestManager(PractRand::Tests::ListOfTests* tests_, PractRand::RNGs::vRNG* known_good_, int max_buffer_amount_ = 1 << (27-10)) : TestManager(tests_, known_good_, max_buffer_amount_) {
		//buffer1.resize(max_buffer_amount + Tests::TestBaseclass::REPEATED_BLOCKS);
		for (auto& test : tests->tests) test->init(known_good);
	}
	void test(Uint64 num_blocks) override {
		while (num_blocks) {
			multithreaded_prep_blocks(num_blocks);
			num_blocks -= main_blocks;
			wait_on_threads();
			alt_buffer.swap(buffer);
			for (auto& test : tests->tests) {
				threads.emplace_back( run_test, test, &alt_buffer[prefix_blocks], main_blocks );
			}
		}
		wait_on_threads();
	}
	void reset(PractRand::RNGs::vRNG* rng_) override {//resets contents for starting a new test run ; if rng is NULL then it will reuse the current RNG
		if (!freshly_created) for (auto& test : tests->tests) test->deinit();
		freshly_created = false;
		for (auto& test : tests->tests) test->init(known_good);
		blocks_to_repeat = 0;
		for (auto& test : tests->tests) {
			int rb = test->get_blocks_to_repeat();
			if (blocks_to_repeat < rb) blocks_to_repeat = rb;
		}
		buffer.resize(max_buffer_amount + blocks_to_repeat);
		alt_buffer.resize(max_buffer_amount + blocks_to_repeat);
		if (rng_) rng = rng_;
		if (!rng) issue_error();
		main_blocks = 0;
		prefix_blocks = 0;
		blocks_so_far = 0;
	}
};
