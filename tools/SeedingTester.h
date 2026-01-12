/*class Seeder_MetaRNG : public PractRand::RNGs::vRNG64 {
public:
	PractRand::RNGs::Polymorphic::hc256 known_good;
	PractRand::RNGs::vRNG *base_rng;
	Uint64 current_seed;
	Uint64 seed_mask;
	int seed_bits;

	std::set<Uint64> unordered_history;
	std::deque<Uint64> history;
	unsigned int history_limit;

	void clear_history() {
		unordered_history.clear();
		history.clear();
	}
	Seeder_MetaRNG(PractRand::RNGs::vRNG *base_rng_, int seed_bits_ = 64) : known_good(PractRand::SEED_NONE), base_rng(base_rng_), seed_bits(seed_bits_), history_limit(1024) {
		if (seed_bits == 64) { seed_mask = ~0ull; }
		else seed_mask = 1ull << seed_bits;
		//current_seed = known_good.raw64();
		//record_seed(current_seed);
	}
	~Seeder_MetaRNG() { delete base_rng; }
	void autoseed() {
		clear_history();
		known_good.autoseed();
		current_seed = known_good.raw64() & seed_mask;
		record_seed(current_seed);
	}
	void seed(PractRand::Uint64 s) {
		clear_history();
		known_good.seed(s);
		current_seed = known_good.raw64() & seed_mask;
		record_seed(current_seed);
	}
	bool record_seed(Uint64 new_seed) {
		if (!unordered_history.insert(new_seed).second) return false;
		current_seed = new_seed;
		history.push_back(current_seed);
		if (history.size() > history_limit) {
			unordered_history.erase(history.front());
			history.pop_front();
		}
		return true;
	}
	void evolve_seed() {
		Uint64 bits_to_try = seed_mask;
		while (true) {
			Uint64 bit = 1ull << known_good.randi32(seed_bits);
			Uint64 new_seed = current_seed ^ bit;
			if (record_seed(new_seed)) return;
			bits_to_try &= ~bit;
			if (!bits_to_try) {// no possible seed exists at hamming distance 1, try higher hamming distances
				while (true) {
					bit = 1ull << known_good.randi32(seed_bits);
					new_seed ^= bit;
					if (record_seed(new_seed)) return;
				}
			}
		}
	}
	Uint64 raw64() {
		base_rng->seed(current_seed);
		Uint64 rv = base_rng->raw64();
		evolve_seed();
		return rv;
	}
	std::string get_name() const {
		std::ostringstream tmp;
		tmp << "SeedingTester(" << base_rng->get_name();
		if (seed_bits != 64) tmp << "," << seed_bits;
		tmp << ")";
		return tmp.str();
	}
	void walk_state(StateWalkingObject *) {}
	static PractRand::RNGs::vRNG *_factory(std::vector<std::string> &params) {
		if (params.size() < 1 || params.size() > 2) { params.push_back("wrong number of parameters - should be SeedingTester(rng) or SeedingTester(rng,seed_bits)"); return NULL; }
		PractRand::RNGs::vRNG *rng = RNG_Factories::create_rng(params[0]);
		int seeding_bits = 64;
		if (params.size() == 2) seeding_bits = atoi(params[1].c_str());
		if (seeding_bits < 32 || seeding_bits > 64) { params.push_back("SeedingTester(rng,seed_bits) - seed_bits should be at least 32 and no more than 64"); return NULL; }
		if (!rng) return NULL;
		return new Seeder_MetaRNG(rng, seeding_bits);
	}
	static void register_name() {
		RNG_Factories::RNG_factory_index["SeedingTester"] = _factory;
	}
};
class FastSeeder_MetaRNG : public PractRand::RNGs::vRNG64 {
public:
	//	Record keeping might have been a limiting factor, so I tried replacing it with a lighter-weight alternative here
	//	Not sure that will actually help though - sometimes seeding the PRNG is the bigger performance issue
	//	If this does well I might replace make this the default behavior of -ttseed64.  
	//	It is doing well.  
	//	despite the name, this calls the regular seed function, NOT seed_fast
	PractRand::RNGs::Polymorphic::hc256 known_good;
	PractRand::RNGs::vRNG *base_rng;
	Uint64 current_seed;
	Uint64 seed_mask;
	int seed_bits;

	typedef Uint16 HashType;
	enum {
		HISTORY_SIZE_L2 = 6,
		HISTORY_SIZE = 1 << HISTORY_SIZE_L2,
		HASHTABLE_BITS = 16, // must be at least 6, and must also be greater than HISTORY_SIZE_L2, but no more than the bits in HashType
		HASHTABLE_SIZE = 1 << (HASHTABLE_BITS - 6),
		MAX_HASHTABLE_BITS = 8 * sizeof(HashType)
	};
	typedef Uint8 _compile_time_assertion[HASHTABLE_BITS <= MAX_HASHTABLE_BITS ? 1 : -1];//if this is an error, then HashType needs to be a bigger integer type (about 10 lines up from here)
	Uint32 warmup;
	Uint32 index;
	HashType ordered_history[HISTORY_SIZE];//history is kept in hashed format - we don't really need full details
	Uint64 unordered_history[HASHTABLE_SIZE];//and unordered history is the actual hashtable

	void clear_history() {
		index = 0;
		warmup = HISTORY_SIZE;
		//for (int i = 0; i < HISTORY_SIZE; i++) ordered_history[i] = 0;
		for (int i = 0; i < HASHTABLE_SIZE; i++) unordered_history[i] = 0;
	}
	FastSeeder_MetaRNG(PractRand::RNGs::vRNG *base_rng_, int seed_bits_ = 64) : known_good(PractRand::SEED_NONE), base_rng(base_rng_), seed_bits(seed_bits_) {
		if (seed_bits == 64) { seed_mask = ~0ull;  }
		else seed_mask = (1ull << seed_bits) - 1;
	}
	~FastSeeder_MetaRNG() { delete base_rng; }
	void autoseed() {
		clear_history();
		known_good.autoseed();
		current_seed = known_good.raw64() & seed_mask;
		record_seed(current_seed);
	}
	void seed(Uint64 seed_low, Uint64 seed_high) {
		clear_history();
		known_good.seed(seed_low, seed_high);
		current_seed = known_good.raw64() & seed_mask;
		record_seed(current_seed);
	}
	HashType calculate_hash(Uint64 seed) {
		seed *= 0x9e3779b97f4a7c15ull; // hope the system has fast 64 bit multiplication, otherwise this isn't much of an optimization
		seed = (seed >> 32) | (seed << 32);
		seed *= 0x9e3779b97f4a7c15ull;
		seed ^= seed >> 32;
		seed *= 0x9e3779b97f4a7c15ull;
		return HashType(seed >> (64 - HASHTABLE_BITS));
	}
	bool record_seed(Uint64 new_seed) {
		HashType hash = calculate_hash(new_seed);
		Uint64 &packed_word = unordered_history[hash >> 6];
		Uint64 bit = 1ull << (hash & 63);
		if (packed_word & bit) return false; // checked if the new hash is available
		packed_word |= bit; // reserve the new hash
		current_seed = new_seed;

		Uint32 masked_index = index++ & (HISTORY_SIZE - 1);
		HashType old_hash = ordered_history[masked_index];
		ordered_history[masked_index] = hash;
		if (warmup) { // check if there is an old hash that needs to be removed
			warmup -= 1;
		}
		else {
			Uint64 &old_packed_word = unordered_history[old_hash >> 6];
			Uint64 old_bit = 1ull << (old_hash & 63);
			old_packed_word &= ~old_bit; // the old hash is now removed
		}
		return true;
	}
	void evolve_seed() {
		Uint64 bits_to_try = seed_mask;
		while (true) {
			Uint64 bit = 1ull << known_good.randi32(seed_bits);
			Uint64 new_seed = current_seed ^ bit;
			if (record_seed(new_seed)) return;
			bits_to_try &= ~bit;
			if (!bits_to_try) {// no possible seed exists at hamming distance 1, try higher hamming distances
				while (true) {
					bit = 1ull << known_good.randi32(seed_bits);
					new_seed ^= bit;
					if (record_seed(new_seed)) return;
				}
			}
		}
	}
	Uint64 raw64() {
		base_rng->seed(current_seed);
		Uint64 rv = base_rng->raw64();
		evolve_seed();
		return rv;
	}
	std::string get_name() const {
		std::ostringstream tmp;
		tmp << "SeedingTester2(" << base_rng->get_name();
		if (seed_bits != 64) tmp << "," << seed_bits;
		tmp << ")";
		return tmp.str();
	}
	void walk_state(StateWalkingObject *) {}
	static PractRand::RNGs::vRNG *_factory(std::vector<std::string> &params) {
		if (params.size() < 1 || params.size() > 2) { params.push_back("wrong number of parameters - should be SeedingTester2(rng) or SeedingTester2(rng,seed_bits)"); return NULL; }
		PractRand::RNGs::vRNG *rng = RNG_Factories::create_rng(params[0]);
		int seeding_bits = 64;
		if (params.size() == 2) seeding_bits = atoi(params[1].c_str());
		if (seeding_bits < 15 || seeding_bits > 64) { params.push_back("SeedingTester2(rng,seed_bits) - seed_bits should be at least 15 and no more than 64 (default is 64)"); return NULL; }
		if (!rng) return NULL;
		return new FastSeeder_MetaRNG(rng, seeding_bits);
	}
	static void register_name() {
		RNG_Factories::RNG_factory_index["SeedingTester2"] = _factory;
	}
};*/
class FastSeeder128_MetaRNG : public PractRand::RNGs::vRNG64 {
public:
	//	Same thing, but 128 bit seed support this time
	PractRand::RNGs::Polymorphic::hc256 known_good;
	PractRand::RNGs::vRNG *base_rng;
	//Uint64 current_seed_low, current_seed_high;
	//Uint64 seed_mask_low, seed_mask_high;
	int seed_bits;

	typedef Uint16 HashType;
	enum {
		HISTORY_SIZE_L2 = 6,
		HISTORY_SIZE = 1 << HISTORY_SIZE_L2,
		HASHTABLE_BITS = 16, // must be at least 6, and must also be greater than HISTORY_SIZE_L2, but no more than the bits in HashType
		HASHTABLE_SIZE = 1 << (HASHTABLE_BITS - 6),
		MAX_HASHTABLE_BITS = 8 * sizeof(HashType)
	};
	typedef Uint8 _compile_time_assertion[HASHTABLE_BITS <= MAX_HASHTABLE_BITS ? 1 : -1];//if this is an error, then HashType needs to be a bigger integer type (about 10 lines up from here)
	Uint32 warmup;
	Uint32 index;
	struct SeedType {
		Uint64 low, high;
		SeedType &operator^=(const SeedType &other) { low ^= other.low; high ^= other.high; return *this; }
		SeedType &operator&=(const SeedType &other) { low &= other.low; high &= other.high; return *this; }
		//SeedType operator~() { SeedType rv; rv.low = ~low; rv.high = ~high; return rv; }
		operator bool() const { return low || high; }
		void set_zero() { low = high = 0; }
	};
	SeedType current_seed, seed_mask;
	HashType ordered_history[HISTORY_SIZE];//history is kept in hashed format - we don't really need full details
	Uint64 unordered_history[HASHTABLE_SIZE];//and unordered history is the actual hashtable

	void clear_history() {
		index = 0;
		warmup = HISTORY_SIZE;
		//for (int i = 0; i < HISTORY_SIZE; i++) ordered_history[i] = 0; // doesn't need to be initialized, because warmup says there's nothing there
		for (int i = 0; i < HASHTABLE_SIZE; i++) unordered_history[i] = 0;
	}
	FastSeeder128_MetaRNG(PractRand::RNGs::vRNG *base_rng_, int seed_bits_ = 128) : known_good(PractRand::SEED_NONE), base_rng(base_rng_), seed_bits(seed_bits_) {
		if (seed_bits == 128) { seed_mask.low = seed_mask.high = ~0ull; }
		else if (seed_bits >= 64) {
			seed_mask.high = (1ull << (seed_bits - 64)) - 1;
			seed_mask.low = ~0ull;
		}
		else {
			seed_mask.high = 0;
			seed_mask.low = (1ull << seed_bits) - 1;
		}
	}
	~FastSeeder128_MetaRNG() { delete base_rng; }
	void autoseed() {
		clear_history();
		known_good.autoseed();
		current_seed.low = known_good.raw64();
		current_seed.high = known_good.raw64();
		current_seed &= seed_mask;
		record_seed(current_seed);
	}
	void seed(Uint64 seed_low, Uint64 seed_high) {
		clear_history();
		known_good.seed(seed_low, seed_high);
		current_seed.low = known_good.raw64();
		current_seed.high = known_good.raw64();
		current_seed &= seed_mask;
		record_seed(current_seed);
	}
	HashType calculate_hash(SeedType seed) {
		Uint64 hash2 = seed.high * 0x854DF6C09CF0E321ull; // hope the system has fast 64 bit multiplication, otherwise this sucks
		Uint64 hash1 = seed.low * 0x9e3779b97f4a7c15ull;
		hash1 ^= hash2;
		hash1 ^= (hash2 << 32) | (hash2 >> 32);
		hash1 *= 0x9e3779b97f4a7c15ull;
		hash1 ^= hash1 >> 32;
		hash1 *= 0x9e3779b97f4a7c15ull;
		return HashType(hash1 >> (64 - HASHTABLE_BITS));
	}
	bool record_seed(SeedType new_seed) {
		HashType hash = calculate_hash(new_seed);
		Uint64 &packed_word = unordered_history[hash >> 6];//hashtable
		Uint64 bit = 1ull << (hash & 63);
		if (packed_word & bit) return false; // checked if the new hash is available
		packed_word |= bit; // reserve the new hash
		current_seed = new_seed;

		Uint32 masked_index = index++ & (HISTORY_SIZE - 1);
		HashType old_hash = ordered_history[masked_index];
		ordered_history[masked_index] = hash;
		if (warmup) { // check if there is an old hash that needs to be removed
			warmup -= 1;
		}
		else {
			Uint64 &old_packed_word = unordered_history[old_hash >> 6];
			Uint64 old_bit = 1ull << (old_hash & 63);
			old_packed_word &= ~old_bit; // the old hash is now removed
		}
		return true;
	}
	SeedType pick_random_bit(int seed_bits) {
		int bit = known_good.rand_i32(seed_bits);
		SeedType rv;
		if (bit < 64) {
			rv.high = 0;
			rv.low = 1ull << bit;
		}
		else {
			rv.low = 0;
			rv.high = 1ull << (bit - 64);
		}
		return rv;
	}
	void evolve_seed() {
		SeedType bits_untried = seed_mask;
		int num_bits_left_to_try = seed_bits;
		while (true) {
			SeedType mutation = pick_random_bit(seed_bits);
			SeedType new_seed = current_seed;
			new_seed ^= mutation;
			if (record_seed(new_seed)) return;
			mutation &= bits_untried;
			if (mutation) {
				num_bits_left_to_try -= 1;
				bits_untried ^= mutation;
			}
			if (!num_bits_left_to_try) {// no possible seed exists at hamming distance 1, try higher hamming distances
				while (true) {
					mutation = pick_random_bit(seed_bits);
					new_seed ^= mutation;
					if (record_seed(new_seed)) return;
				}
			}
		}
	}
	Uint64 raw64() {
		base_rng->seed(current_seed.low, current_seed.high);
		Uint64 rv = base_rng->raw64();
		evolve_seed();
		return rv;
	}
	std::string get_name() const {
		std::ostringstream tmp;
		tmp << "SeedingTester(" << base_rng->get_name();
		if (seed_bits != 128) tmp << "," << seed_bits;
		tmp << ")";
		return tmp.str();
	}
	void walk_state(StateWalkingObject *) {}
	static PractRand::RNGs::vRNG *_factory(std::vector<std::string> &params) {
		if (params.size() < 1 || params.size() > 2) { params.push_back("wrong number of parameters - should be SeedingTester128(rng) or SeedingTester128(rng,seed_bits)"); return NULL; }
		PractRand::RNGs::vRNG *rng = RNG_Factories::create_rng(params[0]);
		int seeding_bits = 128;
		if (params.size() == 2) seeding_bits = atoi(params[1].c_str());
		if (seeding_bits < 16 || seeding_bits > 128) { params.push_back("SeedingTester(rng,seed_bits) - seed_bits should be at least 16 and no more than 128 (default is 128)"); return NULL; }
		if (!rng) return NULL;
		return new FastSeeder128_MetaRNG(rng, seeding_bits);
	}
	static void register_name() {
		RNG_Factories::RNG_factory_index["SeedingTester"] = _factory;
	}
};
class RecursiveSeed64_MetaRNG : public PractRand::RNGs::vRNG64 {
public:
	PractRand::RNGs::vRNG *base_rng;
	PractRand::RNGs::vRNG *alternate_instance;

	RecursiveSeed64_MetaRNG(PractRand::RNGs::vRNG *base_rng_) : base_rng(base_rng_) {}
	~RecursiveSeed64_MetaRNG() { delete base_rng; }
	void autoseed() {
		base_rng->autoseed();
	}
	void seed(Uint64 seed_low, Uint64 seed_high) {
		base_rng->seed(seed_low, 0);
	}
	Uint64 raw64() {
		Uint64 rv = base_rng->raw64();
		base_rng->seed(rv, 0);
		return rv;
	}
	std::string get_name() const {
		std::ostringstream tmp;
		tmp << "RecursiveSeed64(" << base_rng->get_name() << ")";
		return tmp.str();
	}
	void walk_state(StateWalkingObject *) {}
	static PractRand::RNGs::vRNG *_factory(std::vector<std::string> &params) {
		if (params.size() != 1) { params.push_back("wrong number of parameters - should be RecursiveSeed64(rng)"); return NULL; }
		PractRand::RNGs::vRNG *rng = RNG_Factories::create_rng(params[0]);
		int seeding_bits = 64;
		if (!rng) return NULL;
		return new RecursiveSeed64_MetaRNG(rng);
	}
	static void register_name() {
		RNG_Factories::RNG_factory_index["RecursiveSeed64"] = _factory;
	}
};

class EntropyPool_MetaRNG : public PractRand::RNGs::vRNG64 {
public:
	typedef PractRand::Uint64 Transform;
	PractRand::RNGs::Polymorphic::hc256 known_good;
	PractRand::RNGs::vRNG *base_entropy_pool;
	unsigned int min_length, max_length;
	unsigned int history_length;
	std::vector<Uint8> current_seed;
	Transform last_transform;
	std::multiset<Uint64> unordered_history;//hashes only
	std::deque<std::pair<std::multiset<Uint64>::iterator, Transform> > history;//hashes first, then transform applied - newest at front

	//what it should be: (but the current version is good enough)
	//std::map<Uint64,Uint64> unordered_history;// hashes -> positions;
	//std::deque<std::pair<std::map<Uint64,Uint64>::iterator,Transform> > history;//hashes first, then transform applied - newest at front
	//Uint64 position;//starts at 0, incremented after every entropy string

	EntropyPool_MetaRNG(PractRand::RNGs::vRNG *base_entropy_pool_, int min_length_, int max_length_) : known_good(PractRand::SEED_NONE), base_entropy_pool(base_entropy_pool_), min_length(min_length_), max_length(max_length_), history_length(16384) {
		//int len = (min_length + max_length) / 2;
		current_seed.reserve(max_length);
		//current_seed.resize(len);
		//for (int i = 0; i < len; i++) current_seed[i] = known_good.raw8();
		//last_transform = ???;
	}
	~EntropyPool_MetaRNG() { delete base_entropy_pool; }
	void autoseed() {
		known_good.autoseed();
		int len = (min_length + max_length) / 2;
		current_seed.resize(len);
		for (int i = 0; i < len; i++) current_seed[i] = known_good.raw8();
	}
	void seed(Uint64 seed_low, Uint64 seed_high) {
		known_good.seed(seed_low, seed_high);
		int len = (min_length + max_length) / 2;
		current_seed.resize(len);
		for (int i = 0; i < len; i++) current_seed[i] = known_good.raw8();
	}
	Transform pick_random_transform(const std::vector<Uint8> &message) {
		while (true) {
			if (known_good.rand_float() < 0.999) {//toggle bit
				return known_good.rand_i32(message.size() * 8) + (Uint64(0) << 56);
			}
			if (known_good.rand_float() < 0.50) {//insertion
				if (message.size() >= max_length) continue;
				//low 8 bits = value to insert ; next 28 bits = position to insert at ; top 8 bits = action type
				Uint64 position = known_good.rand_i32(message.size() + 1);
				return known_good.raw8() + (position << 8) + (Uint64(1) << 56);
			}
			else {//deletion
				if (message.size() <= min_length) continue;
				Uint64 position = known_good.rand_i32(message.size());
				return message[position] + (position << 8) + (Uint64(2) << 56);
			}
		}
	}
	void apply_transform(std::vector<Uint8> &message, Transform transform) {
		//toggle bit, add byte at end, add byte at begining, remove byte at end, remove byte at begining
		//adds or removals must include the data added or removed in addition to the action
		switch (transform >> 56) {
		case 0://toggle bit
			message[transform >> 3] ^= 1 << (transform & 7);
			break;
		case 1://insert byte
		{
				   int position = (transform >> 8) & ((1ull << 28) - 1);
				   int value = transform & 255;
				   int old_size = message.size();
				   if (position > old_size) { std::printf("internal error - invalid EntropyPool_MetaRNG transform (insert)\n"); std::exit(1); }
				   message.resize(old_size + 1);
				   if (position < old_size) std::memmove(&message[position + 1], &message[position], old_size - position);
				   message[position] = value;
		}
			break;
		case 2://delete byte
		{
				   int position = (transform >> 8) & ((1ull << 28) - 1);
				   int value = transform & 255;
				   int old_size = message.size();
				   if (message[position] != value || position >= old_size) { std::printf("internal error - invalid EntropyPool_MetaRNG transform (deletion)\n"); std::exit(1); }
				   if (position != old_size - 1) std::memmove(&message[position], &message[position + 1], old_size - 1 - position);
				   message.resize(old_size - 1);
		}
			break;
		default:
			std::printf("internal error - unrecognized EntropyPool_MetaRNG transform\n");
			std::exit(1);
		}
	}
	void apply_inverse_transform(std::vector<Uint8> &message, Transform transform) {
		switch (transform >> 56) {
		case 0://reverse a toggle bit by doing the same thing
			apply_transform(message, transform);
			break;
		case 1://reverse an insertion with a deletion
			apply_transform(message, transform + (Uint64(1) << 56));
			break;
		case 2://reverse a deletion with an insertion
			apply_transform(message, transform - (Uint64(1) << 56));
			break;
		}
	}
	PractRand::Uint64 hash_message(const std::vector<Uint8> &message) {
		base_entropy_pool->reset_entropy();
		if (!message.empty())
			base_entropy_pool->add_entropy_N(&message[0], message.size());
		//base_entropy_pool->add_entropy64(0);
		base_entropy_pool->flush_buffers();
		return base_entropy_pool->raw64();
	}
	void check_history_length() {
		while (history.size() > history_length) {
			unordered_history.erase(history.back().first);
			history.pop_back();
		}
	}
	int hamming_distance(const Uint8 *message1, const Uint8 *message2, int n) {
		Uint32 sum = 0;
		for (int i = 0; i < n; i++) sum += PractRand::Internals::count_ones8(message1[i] ^ message2[i]);
		return sum;
	}
	bool check_conflict(const std::vector<Uint8> &message) {
		std::vector<Uint8> rewound = current_seed;
		for (std::deque<std::pair<std::multiset<Uint64>::iterator, Transform> >::iterator it = history.begin(); it != history.end(); it++) {
			//if (message.size() == rewound.size() && !std::memcmp(&message[0], &rewound[0], message.size())) {
			if (message.size() == rewound.size() && !hamming_distance(&message[0], &rewound[0], message.size())) {
				return true;
			}
			apply_inverse_transform(rewound, it->second);
		}
		return false;
	}
	void evolve_seed() {
		while (true) {
			Transform t = pick_random_transform(current_seed);
			std::vector<Uint8> new_seed = current_seed;
			apply_transform(new_seed, t);
			Uint64 hash = hash_message(new_seed);
			std::pair<std::multiset<Uint64>::iterator, std::multiset<Uint64>::iterator> sitr = unordered_history.equal_range(hash);
			if (sitr.first == sitr.second || !check_conflict(new_seed)) {//no conflicts
				current_seed.swap(new_seed);
				std::multiset<Uint64>::iterator it;
				if (sitr.first == unordered_history.end()) it = unordered_history.insert(hash);
				else it = unordered_history.insert(--sitr.first, hash);
				history.push_front(std::pair<std::multiset<Uint64>::iterator, Transform>(it, t));
				check_history_length();
				return;
			}
		}
	}
	PractRand::Uint64 raw64() {
		PractRand::Uint64 rv = hash_message(current_seed);
		evolve_seed();
		return rv;
	}
	std::string get_name() const {
		std::ostringstream tmp;
		tmp << "EntropyPoolingTester(" << base_entropy_pool->get_name() << "," << min_length << "to" << max_length << ")";
		return tmp.str();
	}
	void walk_state(PractRand::StateWalkingObject *) {}
	static PractRand::RNGs::vRNG *_factory(std::vector<std::string> &params) {
		if (params.size() != 3) { params.push_back("wrong number of parameters - should be EntropyPoolingTester(rng,minlength,maxlength)"); return NULL; }
		int minlength = std::atoi(params[1].c_str());
		int maxlength = std::atoi(params[1].c_str());
		if (minlength < 1 || maxlength < 1 || minlength > maxlength || maxlength > 500) { params.push_back("EntropyPoolingTester parameters out of range - 0 < minlength <= maxlength < 500"); return NULL; }
		PractRand::RNGs::vRNG *rng = RNG_Factories::create_rng(params[0]);
		if (!rng) return NULL;
		return new EntropyPool_MetaRNG(rng, minlength, maxlength);
	}
	static void register_name() {
		RNG_Factories::RNG_factory_index["EntropyPoolingTester"] = _factory;
	}
};