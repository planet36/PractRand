#pragma once

//template<typename RNG> double measure_RNG_performance();
//returns mebibytes per second, of calls to raw64() on 64 bit RNGs or raw32 on other RNGs

template<typename RNG> double measure_RNG_performance_16(RNG *rng) {
	constexpr int NUM_CLOCKS_TO_TEST = int(CLOCKS_PER_SEC * .5) + 1;
	//RAW_RNG rng(PractRand::SEED_AUTO);
	Uint16 buffy[1024];
	long clock0 = clock();
	long clock1, clock2;
	while ((clock1 = clock()) == clock0) ;
	int j = 0;
	do {
		for (auto & i : buffy) i = rng->raw16();
		j++;
	} while ((clock2=clock())-clock1 < NUM_CLOCKS_TO_TEST);

	double delta = (clock2 - clock1) / double(CLOCKS_PER_SEC);//seconds
	double amount = j / 256.0;//mebibytes
	double rate = amount / delta;

	//just to make very sure that some smart compiler won't optimize everything away:
	Uint16 a = 0;
	for (const auto i : buffy) a |= i;
	if (a == 0) std::printf("unlikely!");

	return rate;
}
template<typename RNG> double measure_RNG_performance_32(RNG *rng) {
	constexpr int NUM_CLOCKS_TO_TEST = int(CLOCKS_PER_SEC * 0.5) + 1;
	//RAW_RNG rng(PractRand::SEED_AUTO);
	Uint32 buffy[1024] = {0};
	long clock0 = clock();
	long clock1, clock2;
	while ((clock1 = clock()) == clock0) ;
	int j = 0;
	do {
		for (auto & i : buffy) i = rng->raw32();
		j++;
	} while ((clock2=clock())-clock1 < NUM_CLOCKS_TO_TEST);

	double delta = (clock2 - clock1) / double(CLOCKS_PER_SEC);//seconds
	double amount = j / 256.0;//mebibytes
	double rate = amount / delta;

	//just to make very sure that some smart compiler won't optimize everything away:
	Uint32 a = 0;
	for (const auto i : buffy) a |= i;
	if (a == 0) std::printf("unlikely!");

	return rate;
}
template<typename RNG> double measure_RNG_performance_64(RNG *rng) {
	constexpr int NUM_CLOCKS_TO_TEST = int(CLOCKS_PER_SEC * 0.5) + 1;
	//RAW_RNG rng(PractRand::SEED_AUTO);
	Uint64 buffy[1024] = {0};
	long clock0 = clock();
	long clock1, clock2;
	while ((clock1 = clock()) == clock0) ;
	int j = 0;
	do {
		for (auto & i : buffy) i = rng->raw64();
		j++;
	} while ((clock2=clock())-clock1 < NUM_CLOCKS_TO_TEST);

	double delta = (clock2 - clock1) / double(CLOCKS_PER_SEC);//seconds
	double amount = j / 128.0;//mebibytes
	double rate = amount / delta;

	//just to make very sure that some smart compiler won't optimize everything away:
	Uint64 a = 0;
	for (const auto i : buffy) a |= i;
	if (a == 0) std::printf("unlikely!");

	return rate;
}
template<typename RNG> double measure_RNG_performance() {
	RNG rng(PractRand::SEED_AUTO);
	if (RNG::OUTPUT_BITS == 64) return measure_RNG_performance_64(&rng);
	else return measure_RNG_performance_32(&rng);
}
