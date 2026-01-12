/*
A minimal RNG implementation, just enough to make it usable.
	PractRand requires three methods from PRNGs:
		one that returns raw 8, 16, 32, or 64 bit output
		one that returns the name of the PRNG
		and one that allows PractRand to interface with the PRNGs internal state
			if that's not possible for some reason, you can put in an empty implementation that does nothing
			but then default seeding methods won't work, and neither will state serialization/deserialization and anything that relies upon them
			as long as you're ready for those to fail, you can just give walk_state() an empty implementation
			but it's usually easier to give it a real implementation than to have to provide proper seeding functions, let alone the rest
Deliberately flawed, though still better than many platform's default RNGs

Actually, I change this all the time to whatever I'm looking at the moment - I see some PRNG code online, I paste it here to work with it briefly.  

*/

#if defined __GNUC__
#include <x86intrin.h>
#endif

class dummy_rng : public PractRand::RNGs::vRNG32 {
public:
	//declare state
	//Uint64 mix, old_mix, old_rot, old_output, weyl;
	Uint32 a, b, c, d;

	//and any helper methods you want, most of these aren't actually used by any particular algorithm:
	static Uint8 rotate8(Uint8 v, int bits) { return (v << bits) | (v >> (8 - bits)); }
	static Uint16 rotate16(Uint16 v, int bits) { return (v << bits) | (v >> (16 - bits)); }
	static Uint32 rotate32(Uint32 v, int bits) { return (v << bits) | (v >> (32 - bits)); }
	static Uint64 rotate64(Uint64 v, int bits) { return (v << bits) | (v >> (64 - bits)); }
	static Uint32 rotate_right16(Uint32 v, int bits) { return (v >> bits) | (v << (16 - bits)); }
	static Uint32 rotate_right32(Uint32 v, int bits) { return (v >> bits) | (v << (32 - bits)); }
	static Uint64 rotate_right64(Uint64 v, int bits) { return (v >> bits) | (v << (64 - bits)); }
	static void ADC(Uint8 &destination, Uint8 source, Uint8 &carry) { Uint16 x = destination; x += source + carry; destination = x; carry = (x >> (8 * sizeof(destination))) & 1; }
	static void ADC(Uint16 &destination, Uint16 source, Uint16 &carry) { Uint32 x = destination; x += source + carry; destination = x; carry = (x >> (8 * sizeof(destination))) & 1; }
	static void ADC(Uint32 &destination, Uint32 source, Uint32 &carry) { Uint64 x = destination; x += source + carry; destination = x; carry = (x >> (8 * sizeof(destination))) & 1; }
	static Uint64 byteswap64(Uint64 v) { return (Uint64(Uint8(v >> 0)) << 56) | (Uint64(Uint8(v >> 8)) << 48) | (Uint64(Uint8(v >> 16)) << 40) | (Uint64(Uint8(v >> 24)) << 32) | (Uint64(Uint8(v >> 32)) << 24) | (Uint64(Uint8(v >> 40)) << 16) | (Uint64(Uint8(v >> 48)) << 8) | (Uint64(Uint8(v >> 56)) << 0); }
	static Uint32 byteswap32(Uint32 v) { return (Uint32(Uint8(v >> 0)) << 24) | (Uint32(Uint8(v >> 8)) << 16) | (Uint32(Uint8(v >> 16)) << 8) | (Uint32(Uint8(v >> 24)) << 0); }
	static Uint16 byteswap16(Uint16 v) { return rotate16(v, 8); }
	static Uint16 pseudo_bswap16(Uint16 v) { return (Uint16((v >> 0) & 15) << 12) | (Uint16((v >> 4) & 15) << 8) | (Uint16((v >> 8) & 15) << 4) | (Uint16((v >> 12) & 15) << 0); }
	static Uint8 pseudo_bswap8(Uint8 v) { return (((v >> 0) & 3) << 6) | (((v >> 2) & 3) << 4) | (((v >> 4) & 3) << 2) | (((v >> 6) & 3) << 0); }
	static Uint8 add_carry_u16(Uint8 carry_in, Uint16 input_value_A, Uint16 input_value_B, Uint16 *output_address) {
		// gcc doesn't seem to have _addcarry_u16, I don't actually care about performance here atm and it's easier to write a C wrapper than figure out something else
		Uint32 sum = input_value_A;
		sum += input_value_B;
		sum += carry_in;
		*output_address = Uint16(sum);
		return sum >> 16;
	}

	//constructor, if necessary
	dummy_rng() {}
	dummy_rng(PractRand::SEED_AUTO_TYPE) { autoseed(); }
	dummy_rng(PractRand::SEED_NONE_TYPE) {}

	//implement algorithm
	Uint32 raw32() {
		/*Uint64 next_mix = old_rot + old_output;
		old_output = 0x9e3779b97f4a7c15ull * mix;
		old_rot = rotate64(old_mix, 18);
		old_mix = mix ^ weyl;
		mix = next_mix;
		weyl += 0x9e3779b97f4a7c15ull;
		return old_output;//*/

		/*a += b;
		b ^= c;
		c += d;
		d ^= e * 9;
		e += a;
		a = rotate64(a, 12);
		c = rotate64(c, 19);
		return b;//*/
		/*a += b;
		b ^= c;
		c += d;
		d ^= a;
		c = rotate64(c, 19);
		return a;//*/

		/*Uint64 tmp = a + c;
		c = c * 9 + 149;
		//c *= 0x9e3779b9;
		a += b;
		b ^= tmp;
		return tmp;//*/

		/*Uint64 old = a + b;
		a = b * 0x9e3779b97f4a7c15ull;//0x7f4a7c15;
		b = rotate64(c, 27); // ~failed at 38 on sh32, 41+ on sh27
		c = old;
		return a;//*/

		/*Uint32 old = a + c;
		a = b * 0x7f4a7c15; //	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
		b = rotate32(b, 31);//	11	16	18	24	26	28	31	37	37	37	37	36	37	36	37	36	34	37	36	36	36	38*	37	36	36	25	24	24	21	21	16	16
		c++;
		b ^= old;
		return a;//*/

		/*Uint32 old = a + b;
		a = b * 0x7f4a7c15;		//	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
		b = rotate32(old, 15);	//	10	16	16	21	21	21	21	28	35	35	36	37	37	37	37	>34
		return a;//*/

		/*Uint32 old = a + b + c++;
		a = b * 0x7f4a7c15;   //	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
		b = rotate32(old, 8); //	13	16	16	21	21	21	21	29	>42		>42
		return a;//*/

		/*Uint32 old = a + b;
		a = b * 0x7f4a7c15;   //	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	
		b = rotate32(old, 11); //	11	16	16	21	21	21	21	28	34	35	37	36	36	37	37	36	35	36	37	36	36	37	35	34	33	25	22	20	20	20	19	18
		return a;//*/

		/*Uint32 old = a + b;
		a = b * 0x7f4a7c15 ^ c++;  //	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	
		b = rotate32(old, 8);      //	12	16	16	24	26	29	31	37	>45
		return a;//*/

		Uint64 tmp = a * Uint64(0x9e3779b97f4a7c15ull);
		Uint32 tmp_low = Uint32(tmp), tmp_high = Uint32(tmp >> 32);
		a = b ^ tmp_high;
		b = tmp_low;
		return b;

		//a = byteswap64(a * 0x9e3779b97f4a7c15ull);
		//return a;

		/*Uint64 high_result;
		Uint64 old1 = a << 1;
		Uint64 old2 = b;
		a = _umul128(0x9e3779b97f4a7c15ull, a, &high_result);
		b = b * 0x9e3779b97f4a7c15ull + high_result + old1;
		unsigned char carry;
		carry = _addcarry_u64(0, a, 1, &a);
		_addcarry_u64(carry, b, 0, &b);
		return old2;//*/

		/* // this does pretty well on tests for a 48 bit PRNG, failing DC9 at 2^40 bytes, passing everything else
		a = rotate16(a, 5) ^ (b & c);
		b = rotate16(b, 7) ^ c;
		c = rotate16(c, 10) ^ (a + b);
		return b;//*/
	}
	//allow PractRand to be aware of your internal state
	//uses include: default seeding mechanism (individual PRNGs can override), state serialization/deserialization, maybe eventually some avalanche testing tools
	void walk_state(PractRand::StateWalkingObject *walker) {
		//walker->handle(mix); walker->handle(old_rot); walker->handle(old_output); walker->handle(old_mix); walker->handle(weyl);
		walker->handle(a); walker->handle(b); walker->handle(c);
		walker->handle(d);
	}
	//seeding from integers
	//not strictly necessary, in the absence of such a method a default seeding-from-integer path will use walk_state to set the member variables to a hash of the seed
	//note that a separate path exists for seeding-from-output-of-another-PRNG
	/*void seed(Uint64 sv) {
		mix = old_rot = old_output = old_mix = 0;
		weyl = sv;
		for (int i = 0; i < 10; i++) raw16();
	}*/

	//any name you want
	std::string get_name() const { return "dummy_rng"; }
};

/*
The above class is enough to create a PRNG compatible with PractRand.
You can pass that to any non-template function in PractRand expecting a polymorphic RNG and it should work fine.
(some template functions in PractRand require additional metadata or other weirdness)
HOWEVER, that is not enough to allow this (or any other) command line tool to recognize the name of your RNG from the command line.
For that, search for the line mentioning RNG_factory_index["dummy_rng"] in the "main" function,
that line allows it to recognize "dummy" on the command line as corresponding to this class.
*/
