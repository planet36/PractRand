#ifndef __PRACTRAND_RNG_ADAPTORS_H__
#define __PRACTRAND_RNG_ADAPTORS_H__

#if 1
namespace PractRand {
	namespace RNGs {
		namespace Adaptors {
			template<class base_rng> class RAW_TO_LIGHT_WEIGHT_RNG;


			namespace _Internal {
				template<int bits, class base_rng> class ADAPT_OUTPUT_1_TO_ALL;//adds support for all raw output formats (8, 16, 32, and 64 bit output), requires support for one of those already
				template<bool needs_integer_seeding, class base_rng> class ADAPT_SEEDING;//adds support for autoseeding & the like, optionally adds a crude form of seeding from integers
				template<class base_rng> class ADD_DISTRIBUTIONS;//adds support for uniform interger & floating-point distributions, requires support for all raw output bits

				template<bool needs_seeding, class base_rng> class NORMALIZE_SEEDING;//many light-weight RNGs already support seeding
				template<int output_type, int output_bits, class base_rng> class NORMALIZE_OUTPUT;//different RNGs need different output formats added, some don't even need any added

				template<class base_rng> class ADAPT_SEEDING<true, base_rng> : public base_rng {
				public:
					enum { FLAGS = (base_rng::FLAGS & ~RNGs::FLAG::NEEDS_GENERIC_SEEDING) };
					typedef base_rng base_rng_type;
					void seed(Uint64 seed_low, Uint64 seed_high) {StateWalkingObject *walker = Internals::int_to_rng_seeder(seed_low, seed_high); this->walk_state(walker); delete walker;}
					void seed(vRNG *seeder){ StateWalkingObject *walker = Internals::vrng_to_rng_seeder(seeder); this->walk_state(walker); delete walker; }
					void autoseed()        { StateWalkingObject *walker = Internals::get_autoseeder(this); this->walk_state(walker); delete walker; }
				};
				template<class base_rng> class ADAPT_SEEDING<false, base_rng> : public base_rng {
				public:
					enum { FLAGS = (base_rng::FLAGS & ~RNGs::FLAG::NEEDS_GENERIC_SEEDING) };
					typedef base_rng base_rng_type;
					using base_rng_type::seed;
					void seed(vRNG *seeder){ StateWalkingObject *walker = Internals::vrng_to_rng_seeder(seeder); this->walk_state(walker); delete walker; }
					void autoseed()        { StateWalkingObject *walker = Internals::get_autoseeder(this); this->walk_state(walker); delete walker; }
				};

				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<8, base_rng> : public base_rng {
				public:
					enum {OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL};
					typedef base_rng base_rng_type;
					Uint16 raw16() {return this->raw8()  + (((Uint16) this->raw8()) <<  8);}
					Uint32 raw32() {return raw16() + (((Uint32)raw16()) << 16);}
					Uint64 raw64() {return raw32() + (((Uint64)raw32()) << 32);}
				};
				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<16, base_rng> : public base_rng {
				public:
					enum {OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL};
					typedef base_rng base_rng_type;
					Uint8  raw8()  {return (Uint8)this->raw16();}
					Uint32 raw32() {return this->raw16() + (((Uint32)this->raw16()) << 16);}
					Uint64 raw64() {return raw32() + (((Uint64)raw32()) << 32);}
				};
				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<32, base_rng> : public base_rng{
				public:
					enum {OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL};
					typedef base_rng base_rng_type;
					Uint8  raw8()  {return (Uint8)this->raw32();}
					Uint16 raw16() {return (Uint16)this->raw32();}
					Uint64 raw64() {return this->raw32() + (((Uint64)this->raw32()) << 32);}
				};
				template<class base_rng> class ADAPT_OUTPUT_1_TO_ALL<64, base_rng> : public base_rng{
				public:
					enum {OUTPUT_TYPE = OUTPUT_TYPES::NORMAL_ALL};
					typedef base_rng base_rng_type;
					Uint8  raw8()  {return (Uint8)this->raw64();}
					Uint16 raw16() {return (Uint16)this->raw64();}
					Uint32 raw32() {return (Uint32)this->raw64();}
				};

				template<int output_bits, class base_rng> class NORMALIZE_OUTPUT<OUTPUT_TYPES::NORMAL_1, output_bits, base_rng> {
					public:typedef ADAPT_OUTPUT_1_TO_ALL<output_bits, base_rng> t;
				};
				template<int output_bits, class base_rng> class NORMALIZE_OUTPUT<OUTPUT_TYPES::NORMAL_ALL, output_bits, base_rng> {
					public:typedef base_rng t;
				};

				template<class base_rng> class ADD_DISTRIBUTIONS : public base_rng {
				public:
//						enum {DISTRIBUTIONS_TYPE = DISTRIBUTIONS_TYPE__NORMAL};
					Uint32 rand_i32 ( Uint32 max ) {
						Uint32 mask, tmp;
						max -= 1;
						mask = max;
						mask |= mask >> 1; mask |= mask >>  2; mask |= mask >> 4;
						mask |= mask >> 8; mask |= mask >> 16;
						while (1) {
							tmp = this->raw32() & mask;
							if (tmp <= max) return tmp;
						}
					}

					Uint64 rand_i64 ( Uint64 max ) {
						Uint64 mask, tmp;
						max -= 1;
						mask = max;
						mask |= mask >> 1; mask |= mask >>  2; mask |= mask >>  4;
						mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;
						while (1) {
							tmp = this->raw64() & mask;
							if (tmp <= max) return tmp;
						}
					}
					Uint32 rand_i32_fast ( Uint32 max ) {return randi32_fast_implementation(this->raw32(), max);}

					//random floating point numbers:
					float rand_float() { return rand_float_implementation(this->raw32()); }
					float rand_float( float m ) { return rand_float() * m; }
					double rand_double() { return rand_double_implementation(this->raw64()); }
					double rand_double ( double m ) { return rand_double() * m; }

					//Boost / C++0x TR1 compatibility:
#if defined PRACTRAND_BOOST_COMPATIBILITY
					typedef Uint64 result_type;
					result_type operator()() {return this->raw64();}
					static const bool has_fixed_value = true;
					static const result_type min_value = 0;
					static const result_type max_value = ~(result_type)0;
					result_type min() const {return min_value;}
					result_type max() const {return max_value;}
#endif
				};

				template<class base_rng> class NORMALIZE_ALL {
				public:
					typedef ADAPT_SEEDING< (base_rng::FLAGS & PractRand::RNGs::FLAG::NEEDS_GENERIC_SEEDING),
						ADD_DISTRIBUTIONS<
							typename NORMALIZE_OUTPUT<base_rng::OUTPUT_TYPE, base_rng::OUTPUT_BITS, base_rng>::t 
						> 
					> t;
				};

//*/
			}//namespace _Internal


			template<class base_rng> class RAW_TO_LIGHT_WEIGHT_RNG : public _Internal::NORMALIZE_ALL<base_rng>::t {
			public:
				RAW_TO_LIGHT_WEIGHT_RNG(SEED_AUTO_TYPE) { this->autoseed(); }
				RAW_TO_LIGHT_WEIGHT_RNG(SEED_NONE_TYPE) {}
				RAW_TO_LIGHT_WEIGHT_RNG(Uint64 seed_low, Uint64 seed_high=0) {this->seed(seed_low, seed_high);}
				RAW_TO_LIGHT_WEIGHT_RNG(vRNG *seeder) {this->seed(seeder);}
			};
			//to do:
			//template<class base_rng> class RAW_TO_POLYMORPHIC_RNG;
			//template<class base_rng> class NORMAL_TO_POLYMORPHIC_RNG;
		}//namespace Adaptors
	}//namespace RNGs
}//namespace PractRand
#endif

#endif //__PRACTRAND_RNG_ADAPTORS_H__
