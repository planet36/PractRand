#pragma once

#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"

#if 1
namespace PractRand::RNGs::Adaptors {
			template<class base_rng> class NORMALIZE_OUTPUT;
			template<class base_rng> class NORMALIZE_SEEDING;
			template<class base_rng> class NORMALIZE_DISTRIBUTIONS;


			template<class base_rng> class NORMALIZE;
			template<class base_rng> class NORMALIZE_OUTPUT_TYPE;
			//template<class base_rng> class _TEMPLATIZE_OUTPUT;
			template<class base_rng> class NORMALIZE_SEEDING_TYPE;
			template<class base_rng> class NORMALIZE_DISTRIBUTIONS_TYPE;

			namespace Internal {
				template<class base_rng, int bits> class ADAPT_OUTPUT_1_TO_ALL;

				template<class base_rng, bool needs_int_seeding> class ADAPT_SEEDING;
				template<class base_rng> class ADAPT_SEEDING<base_rng,true> : public base_rng {
				public:
					static constexpr int FLAGS = base_rng::FLAGS & ~ RNGs::FLAG::NEEDS_GENERIC_SEEDING;
					//static constexpr int RNG_WRAPPER_LEVEL = base_rng::RNG_WRAPPER_LEVEL+1;
					using base_rng_type = base_rng;
					void seed     (uint64_t seed) {StateWalkingObject* walker = int_to_rng_seeder(seed); this->walk_state(walker); delete walker;}
					void seed     (vRNG* seeder){StateWalkingObject* walker = vrng_to_rng_seeder(seeder); this->walk_state(walker); delete walker;}
					void autoseed ()            {StateWalkingObject* walker = get_autoseeder(this); this->walk_state(walker); delete walker;}
				};
				template<class base_rng> class ADAPT_SEEDING<base_rng,false> : public base_rng {
				public:
					static constexpr int FLAGS = base_rng::FLAGS & ~ RNGs::FLAG::NEEDS_GENERIC_SEEDING;
					//static constexpr int RNG_WRAPPER_LEVEL = base_rng::RNG_WRAPPER_LEVEL+1;
					using base_rng_type = base_rng;
					using base_rng :: seed;
					void seed     (vRNG* seeder){StateWalkingObject* walker = vrng_to_rng_seeder(seeder); this->walk_state(walker); delete walker;}
					void autoseed ()            {StateWalkingObject* walker = get_autoseeder(this); this->walk_state(walker); delete walker;}
				};

				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<base_rng, 8> : public base_rng {
				public:
					static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;
					//static constexpr int RNG_WRAPPER_LEVEL = base_rng::RNG_WRAPPER_LEVEL+1;
					using base_rng_type = base_rng;
					uint16_t raw16() {return this->raw8()  + (static_cast<uint16_t>(this->raw8()) <<  8);}
					uint32_t raw32() {return raw16() + (static_cast<uint32_t>(raw16()) << 16);}
					uint64_t raw64() {return raw32() + (static_cast<uint64_t>(raw32()) << 32);}
				};
				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<base_rng, 16> : public base_rng {
				public:
					static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;
					//static constexpr int RNG_WRAPPER_LEVEL = base_rng::RNG_WRAPPER_LEVEL+1;
					using base_rng_type = base_rng;
					uint8_t  raw8()  {return static_cast<uint8_t>(this->raw16());}
					uint32_t raw32() {return this->raw16() + (static_cast<uint32_t>(this->raw16()) << 16);}
					uint64_t raw64() {return raw32() + (static_cast<uint64_t>(raw32()) << 32);}
				};
				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<base_rng, 32> : public base_rng {
				public:
					static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;
					//static constexpr int RNG_WRAPPER_LEVEL = base_rng::RNG_WRAPPER_LEVEL+1;
					using base_rng_type = base_rng;
					uint8_t  raw8()  {return static_cast<uint8_t>(this->raw32());}
					uint16_t raw16() {return static_cast<uint16_t>(this->raw32());}
					uint64_t raw64() {return this->raw32() + (static_cast<uint64_t>(this->raw32()) << 32);}
				};
				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<base_rng, 64> : public base_rng {
				public:
					static constexpr int OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL;
					//static constexpr int RNG_WRAPPER_LEVEL = base_rng::RNG_WRAPPER_LEVEL+1;
					using base_rng_type = base_rng;
					uint8_t  raw8()  {return static_cast<uint8_t>(this->raw64());}
					uint16_t raw16() {return static_cast<uint16_t>(this->raw64());}
					uint32_t raw32() {return static_cast<uint32_t>(this->raw64());}
				};


				template<class base_rng, int output_type, int output_bits> class NORMALIZE_OUTPUT_HELPER;
				template<class base_rng, int output_bits> class NORMALIZE_OUTPUT_HELPER<base_rng, OUTPUT_TYPES::NORMAL_1, output_bits> {
					public:using t = ADAPT_OUTPUT_1_TO_ALL<base_rng, output_bits>;
				};
				template<class base_rng, int output_bits> class NORMALIZE_OUTPUT_HELPER<base_rng, OUTPUT_TYPES::NORMAL_ALL, output_bits> {
					public:using t = base_rng;
				};

				template<class base_rng> class ADD_DISTRIBUTIONS : public base_rng {
				public:
					//static constexpr int DISTRIBUTIONS_TYPE = DISTRIBUTIONS_TYPE__NORMAL;
					uint32_t randi ( uint32_t max ) {
						uint32_t mask = 0, tmp = 0;
						max -= 1;
						mask = max;
						mask |= mask >> 1; mask |= mask >>  2; mask |= mask >> 4;
						mask |= mask >> 8; mask |= mask >> 16;
						while (true) {
							tmp = this->raw32() & mask;
							if (tmp <= max) return tmp;
						}
					}
					uint32_t randi ( uint32_t min, uint32_t max ) {return randi(max-min) + min;}

					uint64_t randli ( uint64_t max ) {
						uint64_t mask = 0, tmp = 0;
						max -= 1;
						mask = max;
						mask |= mask >> 1; mask |= mask >>  2; mask |= mask >>  4;
						mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;
						while (true) {
							tmp = this->raw64() & mask;
							if (tmp <= max) return tmp;
						}
					}
					uint64_t randli ( uint64_t min, uint64_t max ) {return randli(max-min) + min;}

					uint32_t randi_fast ( uint32_t max ) {return randi_fast_implementation(this->raw32(), max);}
					uint32_t randi_fast ( uint32_t min, uint32_t max ) {return randi_fast(max-min) + min;}

					//random floating point numbers:
					float randf ( ) { return float(this->raw32() * (1.0 / 4294967296.0)); }
					float randf ( float m ) { return randf() * m; }
					float randf ( float min, float max ) { return randf() * (max-min) + min; }

					double randlf ( ) { return static_cast<double>(this->raw64()) * (1.0 / 18446744073709551616.0); }
					double randlf ( double m ) { return randlf() * m; }
					double randlf ( double min, double max ) { return randlf() * (max-min) + min; }

					//Boost / C++0x TR1 compatibility:
#if defined PRACTRAND_BOOST_COMPATIBILITY
					typedef uint64_t result_type;
					result_type operator()() {return this->raw64();}
					static constexpr bool has_fixed_value = true;
					static constexpr result_type min_value = 0;
					static constexpr result_type max_value = ~(result_type)0;
					result_type min() const {return min_value;}
					result_type max() const {return max_value;}
#endif
				};
				template<class base_rng, bool needs_distributions_added> class NORMALIZE_DISTRIBUTIONS_HELPER;
				template<class base_rng> class NORMALIZE_DISTRIBUTIONS_HELPER<base_rng,true> {
				public:using t = ADD_DISTRIBUTIONS<base_rng>;
				};
				template<class base_rng> class NORMALIZE_DISTRIBUTIONS_HELPER<base_rng,false> {
				public:using t = base_rng;
				};

			}//namespace Internal

			template<class base_rng> class NORMALIZE_SEEDING_TYPE {
			public: using t =
				Internal::ADAPT_SEEDING<
					base_rng, static_cast<bool>(base_rng::FLAGS & RNGs::FLAG::NEEDS_GENERIC_SEEDING)
				>;
			//public:typedef typename base_rng t;
			};
			template<class base_rng> class NORMALIZE_SEEDING : public NORMALIZE_SEEDING_TYPE<base_rng>::t {};

			template<class base_rng> class NORMALIZE_OUTPUT_TYPE {
				public:using t = Internal::NORMALIZE_OUTPUT_HELPER<base_rng,base_rng::OUTPUT_TYPE, base_rng::OUTPUT_BITS>::t;
			};
			template<class base_rng> class NORMALIZE_OUTPUT : public NORMALIZE_OUTPUT_TYPE<base_rng>::t {};

			template<class base_rng> class NORMALIZE_DISTRIBUTIONS_TYPE {
				//public:typedef typename Internal::_NORMALIZE_DISTRUBTIONS_HELPER<base_rng,bool(base_rng::DISTRUBTIONS_TYPE & DISTRIBUTIONS_TYPE__NORMAL)>::t t;
				public:using t = Internal::ADD_DISTRIBUTIONS<base_rng>;
			};
			template<class base_rng> class NORMALIZE_DISTRIBUTIONS : public NORMALIZE_DISTRIBUTIONS_TYPE<base_rng>::t {};

			template<class base_rng> class NORMALIZE {
				public:using t = NORMALIZE_SEEDING_TYPE<
					typename NORMALIZE_DISTRIBUTIONS_TYPE<
						typename NORMALIZE_OUTPUT_TYPE<base_rng>::t
					>::t
				>::t;
			};
			template<class base_rng> class RAW_TO_LIGHT_WEIGHT_RNG : public NORMALIZE<base_rng>::t {
			public:
				explicit RAW_TO_LIGHT_WEIGHT_RNG(SEED_AUTO_TYPE) {this->autoseed();}
				explicit RAW_TO_LIGHT_WEIGHT_RNG(SEED_NONE_TYPE) {}
				explicit RAW_TO_LIGHT_WEIGHT_RNG(uint64_t s) {this->seed(s);}
				explicit RAW_TO_LIGHT_WEIGHT_RNG(vRNG* seeder) {this->seed(seeder);}
			};
			//to do:
			//template<class base_rng> class RAW_TO_POLYMORPHIC_RNG;
			//template<class base_rng> class NORMAL_TO_POLYMORPHIC_RNG;
}//namespace PractRand
#endif
