namespace PractRand {
	namespace Tests {
		class Coup16 final : public TestBaseclass {
		protected:
			enum { S = 65536 / 32 };
			Uint32 flags[S];
			FixedSizeCount<Uint8, 65536> counts;
			//to do: also measure sequential correlation, or do something like BCFN on whether or not each sample is >= the expected value
		public:
			Coup16() {}
			virtual void init(PractRand::RNGs::vRNG *known_good) override;
			virtual std::string get_name() const override;
			virtual void get_results(std::vector<TestResult> &results) override;
			virtual void test_blocks(TestBlock *data, int numblocks) override;
		};
	}//Tests
}//PractRand
