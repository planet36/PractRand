#pragma once

#include "PractRand/test_helpers.h"
#include "PractRand/tests.h"

namespace PractRand::Tests {
		class DistC6 : public TestBaseclass {
		public:
			explicit DistC6 ( int length_ = 9, int unitsL_ = 0,
				int bits_clipped_0_ = 1,
				int bits_clipped_1_ = 0,
				int bits_clipped_2_ = 0 );
			void init( PractRand::RNGs::vRNG* known_good ) override;
			[[nodiscard]] std::string get_name() const override;
			void get_results ( std::vector<TestResult>& results ) override;

			void test_blocks(TestBlock* data, int numblocks) override;

		protected:
			static constexpr int ENABLE_REORDER = 1;
			static constexpr int ENABLE_8_BIT_BYPASS = 1;
			//configuration:
			int length;
			int unitsL;
			int bits_clipped_0;
			int bits_clipped_1;
			int bits_clipped_2;
			//precalcs:
			int bits_per_sample;
			int size;
			uint32_t mask_pre{uint32_t(-1)};
			uint32_t lookup_table[ENABLE_8_BIT_BYPASS ? 256 : 65]{};//reorder_bits(reorder_codes(transform_bitcount( X ))) or reorder_bits(reorder_codes(transform_bitcount(count_bits8( X ))))
			uint8_t reorder_codes_table[64]{};
			uint8_t unreorder_codes_table[64]{};
			//state:
			int warmup{};
			uint32_t last_index{};
			VariableSizeCount<uint8_t> counts;
			//internal helpers:
			[[nodiscard]] virtual int transform_bitcount ( int bit_count ) const final;
			static uint32_t _reorder_bits ( uint32_t rcode, int bits_per_sample, int length ) ;
			[[nodiscard]] uint32_t reorder_bits( uint32_t rcode ) const {return _reorder_bits(rcode, bits_per_sample, length);}
			[[nodiscard]] uint32_t unreorder_bits ( uint32_t rbcode ) const {return _reorder_bits(rbcode, length, bits_per_sample);}
			[[nodiscard]] uint32_t reorder_codes   ( uint32_t code ) const {return reorder_codes_table[code];}
			[[nodiscard]] uint32_t unreorder_codes ( uint32_t rcode ) const {return unreorder_codes_table[rcode];}
			void generate_reorder_codes ();
			[[nodiscard]] uint32_t _advance_index ( uint32_t index, int rbcode ) const {
				if constexpr (ENABLE_REORDER)
					return ((index & mask_pre) << 1) | rbcode ;
				else
					return ((index & mask_pre) << bits_per_sample) | rbcode ;
			}
			void advance_index ( int bit_count ) {
				last_index = _advance_index(last_index, lookup_table[bit_count]);}
		};
		class DistC7 final : public DistC6 {
		public:
			explicit DistC7(int length_ = 9, int unitsL_ = 0,
				int bits_clipped_0_ = 1,
				int bits_clipped_1_ = 0,
				int bits_clipped_2_ = 0
				);
			void init(PractRand::RNGs::vRNG* known_good) override;
			[[nodiscard]] std::string get_name() const override;
			void get_results(std::vector<TestResult>& results) override;

			void test_blocks(TestBlock* data, int numblocks) override;
		protected:
			//configuration:
			//precalcs:
			//state:
			bool odd{};
			VariableSizeCount<uint8_t> odd_counts;
			//internal helpers:
		};
}//PractRand
