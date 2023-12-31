#pragma once

#include "PractRand/test_batteries.h"

namespace PractRand::Tests::Transforms {
			class multiplex : public TestBaseclass {
			protected:
				ListOfTests subtests;
				std::string name;
				Uint64 blocks_already;//# of blocks outputed to subtests
			public:
				multiplex ( const char *name_, const ListOfTests &testlist );
				void init( RNGs::vRNG *known_good ) override;

				//virtual double get_result();
				void get_results ( std::vector<TestResult> &results ) override;

				virtual bool recommend_subtest_tree_descent() const;
				void deinit() override;
				~multiplex ( ) override;
				std::string get_name() const override;
				void test_blocks(TestBlock *data, int numblocks) override;
				int get_blocks_to_repeat() const override;

				virtual int get_num_children() const;
				virtual TestBaseclass *get_child(int index) const;
				const ListOfTests &_get_children() const {return subtests;}
				virtual Uint64 get_blocks_passed_through(int index) const;// {return blocks_already;}
				//virtual std::string get_child_name  (int index) const;
				//virtual double      get_child_result(int index);
			};
			class switching final : public multiplex {
				std::vector<Uint64> lengths;
				std::vector<Uint64> blocks_already_per;
				Uint64 total_length;
				Uint64 phase;
				unsigned int which;
			public:
				switching( const char *name_, const ListOfTests &testlist,
					std::vector<Uint64> lengths_ );
				switching( const char *name_, const ListOfTests &testlist,
					Uint64 length );
				void init( RNGs::vRNG *known_good ) override;
				void test_blocks(TestBlock *data, int numblocks) override;
			//	virtual double get_result();
				Uint64 get_blocks_passed_through(int index) const override;
			};
			class Transform_Baseclass : public multiplex {
			private:
				int flush_size;
			protected:
				std::vector<TestBlock> buffered;
				int leftovers;//if non-zero then the last buffered block is treated as unready
			public:
				Transform_Baseclass( const char *name_, const ListOfTests &testlist, int flush_size_ = 4*1024 ) :
					multiplex(name_, testlist), flush_size(flush_size_)
				{}
				void init( RNGs::vRNG *known_good ) override;
				void test_blocks(TestBlock *data, int numblocks) override = 0;
				//virtual double get_result() {flush(true);return Transforms::multiplex::get_result();}
				void get_results ( std::vector<TestResult> &results ) override {flush(true); Transforms::multiplex::get_results(results);}
				virtual void flush(bool aggressive = true) final;
			};
			class shrink : public Transform_Baseclass {
			protected:
				Uint8 sbox[65536];
				int pattern;//
				int outbitsL;//1,2,4, or 8
				int unitsL;//0= 8 bit input words, 1 = 16 bit input words
			public:
				shrink ( const char *name_, const ListOfTests &testlist, int pattern_, int outbitsL_=1, int unitsL_=0 )
				:
					Transform_Baseclass(name_, testlist),
					pattern(pattern_),
					outbitsL(outbitsL_),
					unitsL(unitsL_)
				{}
				void init( RNGs::vRNG *known_good ) override;
				std::string get_name() const override;
				void test_blocks(TestBlock *data, int numblocks) override;
			};
			class lowbits final : public Transform_Baseclass {
			protected:
				int lowbitsL;//0 = 1 bit of output per input word, 1 = 2 bits of output per input word, 2=4,3=8,4=15,5=32
				int unitsL;//-1= 4 bit input words, 0= 8 bit input words, 1 = 16 bit input words, 2 = 32, 3 = 64
			public:
				lowbits ( const char *name_, const ListOfTests &testlist, int lowbitsL_=1, int unitsL_=0 )
				:
					Transform_Baseclass(name_, testlist),
					lowbitsL(lowbitsL_),
					unitsL(unitsL_)
				{}
				void init( RNGs::vRNG *known_good ) override;
				std::string get_name() const override;
				void test_blocks(TestBlock *data, int numblocks) override;
			};
			/*class bitsyr8x8 : public Transform_Baseclass {
			public:
				bitsyr8x8 ( const char *name_, ListOfTests &testlist )
					: Transform_Baseclass(name_, testlist) {}

				virtual ~bitsyr8x8 ( );
				virtual std::string get_name() const;
				virtual void test_blocks(TestBlock *data, int numblocks);
			};*/
			class FirstNofM final : public Transform_Baseclass {
			protected:
				int bytes_used;
				int bytes_stride;
				int input_phase;
			public:
				FirstNofM ( const char *name_, const ListOfTests &testlist, int bytes_used_, int bytes_stride_ )
				:
					Transform_Baseclass(name_, testlist),
					bytes_used(bytes_used_),
					bytes_stride(bytes_stride_)
				{}
				void init( RNGs::vRNG *known_good ) override;
				std::string get_name() const override;
				void test_blocks(TestBlock *data, int numblocks) override;
			};
}//PractRand
