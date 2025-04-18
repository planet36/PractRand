#pragma once

#include "PractRand/RNGs/other/fibonacci.h"
#include "PractRand/RNGs/other/indirection.h"
#include "PractRand/RNGs/other/mult.h"
#include "PractRand/RNGs/other/simple.h"
#include "PractRand/RNGs/other/transform.h"
#include "PractRand/rng_basics.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace Special_RNGs {
	template<typename Word>
	class _stdin_reader {
		static void read_failed() {
			(void)std::fprintf(stderr, "error reading standard input\n");
			std::exit(0);
		}
		static constexpr int BUFF_SIZE = 32768 / sizeof(Word);
		Word *pos, *end;
		bool ended{false};
		Word buffer[BUFF_SIZE];
		void refill() {
			std::size_t n = std::fread(&buffer[0], sizeof(Word), BUFF_SIZE, stdin);
			if (n < BUFF_SIZE) ended = true;
			if (!n) read_failed();
			pos = &buffer[0];
			end = &buffer[n];
		}
	public:
		_stdin_reader() { refill(); }
		Word read() { if (pos == end) refill(); return *(pos++); }
	};
	class RNG_stdin : public PractRand::RNGs::vRNG8 {
		_stdin_reader<PractRand::Uint8> source;
		PractRand::Uint8 raw8() override {
			return source.read();
		}
		PractRand::Uint64 get_flags() const override { return PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED | PractRand::RNGs::FLAG::STATE_UNAVAILABLE; }
		std::string get_name() const override { return "RNG_stdin"; }
		int get_native_output_size() const override {return -1;}
		void walk_state(PractRand::StateWalkingObject *) override {}
	};
	class RNG_stdin8 : public PractRand::RNGs::vRNG8 {
		_stdin_reader<PractRand::Uint8> source;
		PractRand::Uint8 raw8() override { return source.read(); }
		PractRand::Uint64 get_flags() const override { return PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED | PractRand::RNGs::FLAG::STATE_UNAVAILABLE; }
		std::string get_name() const override { return "RNG_stdin8"; }
		void walk_state(PractRand::StateWalkingObject *) override {}
	};
	class RNG_stdin16 : public PractRand::RNGs::vRNG16 {
		_stdin_reader<PractRand::Uint16> source;
		PractRand::Uint16 raw16() override { return source.read(); }
		PractRand::Uint64 get_flags() const override { return PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED | PractRand::RNGs::FLAG::STATE_UNAVAILABLE; }
		std::string get_name() const override { return "RNG_stdin16"; }
		void walk_state(PractRand::StateWalkingObject *) override {}
	};
	class RNG_stdin32 : public PractRand::RNGs::vRNG32 {
		_stdin_reader<PractRand::Uint32> source;
		PractRand::Uint32 raw32() override { return source.read(); }
		PractRand::Uint64 get_flags() const override { return PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED | PractRand::RNGs::FLAG::STATE_UNAVAILABLE; }
		std::string get_name() const override { return "RNG_stdin32"; }
		void walk_state(PractRand::StateWalkingObject *) override {}
	};
	class RNG_stdin64 : public PractRand::RNGs::vRNG64 {
		_stdin_reader<PractRand::Uint64> source;
		PractRand::Uint64 raw64() override { return source.read(); }
		PractRand::Uint64 get_flags() const override { return PractRand::RNGs::FLAG::SEEDING_UNSUPPORTED | PractRand::RNGs::FLAG::STATE_UNAVAILABLE; }
		std::string get_name() const override { return "RNG_stdin64"; }
		void walk_state(PractRand::StateWalkingObject *) override {}
	};
}

namespace PractRand::RNG_Sets {
			const char *recommended_rngs[] = {
				"hc256", "trivium",
				"efiix64x384", "efiix32x384", "efiix16x384", "efiix8x384",
				"isaac64x256", "isaac32x256",
				"chacha(8)", "salsa(8)",
				"arbee",
				"xsm64", "xsm32",
				"jsf64", "jsf32", "sfc64", "sfc32", "sfc16",
				"rarns16", "rarns32", "rarns64",
				"mt19937",
				nullptr
			};
			const int num_recommended_rngs = sizeof(recommended_rngs) / sizeof(recommended_rngs[0]) - 1;
			const char *nonrecommended_simple[] = {
				"xorshift32", "xorshift64", "xorshift32of128", "xoroshiro128plus",
				"xorwow32of96", "xorwow32x6", "xsalta16x3", "xsaltb16x3", "xsaltc16x3",
				"sapparot", "sap16of48", "sap32of96",
				"flea32x1", "jsf16",
				"sfc_v1_16", "sfc_v1_32", "sfc_v2_16", "sfc_v2_32", "sfc_v3_16", "sfc_v3_32",
				"simpleA", "simpleB", "simpleC", "simpleD", "simpleE", "simpleF", "simpleG",
				"trivium_weakenedA", "trivium_weakenedB",
				"mo_Lesr32", "mo_ResrRers32", "mo_Rers32of64", "mo_Resdra32of64",
				"murmlacish",
				"gjishA", "gjishB", "gjishC", "gjishD",
				"ara16", "ara32", "arx16", "arx32", "hara16", "harx16", "learx16", "hlearx16", "alearx16", "arac16", "arxc16",
				nullptr
			};
			const int num_nonrecommended_simple = sizeof(nonrecommended_simple) / sizeof(nonrecommended_simple[0]) - 1;
			const char *nonrecommended_nonlcg[] = {
				"garthy16", "garthy32", "binarymult16", "binarymult32", "rxmult16", "multish3x32", "multish4x16",
				"mwrca16", "mwrca32", "mwrcc16", "mwrcc32", "mwrcca16", "mwrcca32",
				"old_mwlac16", "mwlac_varA", "mwlac_varB", "mwlac_varC", "mwlac_varD", "mwlac_varE",
				"mo_Cmfr", "mo_Cmr32of64", "mulcr16", "mulcr32", "mmr16", "mmr32",
				nullptr
			};
			const int num_nonrecommended_nonlcg = sizeof(nonrecommended_nonlcg) / sizeof(nonrecommended_nonlcg[0]) - 1;
			const char *nonrecommended_lcgish[] = {
				"lcg(16,32)", "lcg(16,40)", "lcg(16,48)", "lcg(16,56)", "lcg(16,64)", "lcg(16,72)", "lcg(16,80)", "lcg(16,96)", "lcg(16,112)",
				"xlcg(16,32)", "xlcg(16,40)", "xlcg(16,48)", "xlcg(16,56)", "xlcg(16,64)", "xlcg(16,72)",
				"clcg(16,64)", "clcg(16,68)", "clcg(16,72)", "clcg(16,76)", "clcg(16,80)", "clcg(16,84)",
				"cxlcg(16,64)", "cxlcg(16,68)", "cxlcg(16,72)",
				"bblcg(32,160,32)", "bblcg(32,192,32)", "bblcg(32,224,32)", "bblcg(32,256,32)", "bblcg(32,288,32)",
				"pcg32_norot", "pcg32", "cmrg32of192", "xsh_lcg_bad",
				nullptr
			};
			const int num_nonrecommended_lcgish = sizeof(nonrecommended_lcgish) / sizeof(nonrecommended_lcgish[0]) - 1;
			const char *nonrecommended_cbuf[] = {
				"mm32", "mm32_awc", "mm16of32", "mm16of32_awc", "mm4691",
				"cbuf_accum", "cbuf_accum_big", "cbuf_2accum_small", "cbuf_2accum", "dual_cbuf_small", "dual_cbuf", "dual_cbufa_small", "dual_cbuf_accum",
				"fibmul16of32", "fibmul32of64", "fibmulmix16", "ranrot32small", "ranrot32", "ranrot32big", "ranrot3tap32small", "ranrot3tap32", "ranrot3tap32big", "ranrot32hetsmall", "ranrot32het", "ranrot32hetbig",
				"mt19937_unhashed", "salsa(3)", "chacha(3)", "salsa(4)", "chacha(4)",
				nullptr
			};
			const int num_nonrecommended_cbuf = sizeof(nonrecommended_cbuf) / sizeof(nonrecommended_cbuf[0]) - 1;
			const char *nonrecommended_indirect[] = {
				"ibaa8(1)", "ibaa8(2)", "ibaa8(3)", "ibaa16(1)", "ibaa16(2)", "ibaa16(3)", "ibaa32(1)", "ibaa32(2)", "ibaa32(3)",
				"rc4_weakenedA", "rc4_weakenedB", "rc4_weakenedC", "rc4_weakenedD", "rc4",
				"isaac16(2)", "isaac16(3)", "isaac16(4)", "isaac32(2)", "isaac32(3)", "isaac32(4)",
				"genindA(5)", "genindA(7)", "genindA(9)", "genindB(1)", "genindB(2)", "genindB(3)", "genindB(4)",
				"genindC(2)", "genindC(3)", "genindC(4)", "genindC(5)", "genindD(6)", "genindD(9)",
				"genindE(1)", "genindE(2)", "genindE(3)", "genindF(2)", "genindF(3)", "genindF(4)", "genindF(5)",
				nullptr
			};
			const int num_nonrecommended_indirect = sizeof(nonrecommended_indirect) / sizeof(nonrecommended_indirect[0]) - 1;
}


namespace RNG_Factories {
	std::map<std::string, PractRand::RNGs::vRNG *(*)(std::vector<std::string> &params)> RNG_factory_index;
/*	struct RNG_Description {
		const char *name;
		const char *full_name;
		PractRand::RNGs::vRNG *(*factory)(std::vector<std::string> &);
		const char *description;
		const char *param_description;
		int quality;
		int quality_subscores;//decimal digits for each quality subscore - empirical, cycle lengths, statespace, analysis
		int speed;
		int word;
		int size;
		Uint64 flags;
	};*/

	bool parse_argument_list(const std::string &raw, std::string &name, std::vector<std::string> &parameters) {
		// "BDS(BDS(chacha(8),4),8)" -> "BDS", "BDS(chacha(8),4)", "8"
		// that is, function name, parameter 1, parameter 2, ...
		parameters.clear();
		const char *base = raw.c_str();
		const char *tmp = strpbrk(base, "(,)");
		if (!tmp) {
			name = raw;
			return true;
		}
		if (*tmp != '(') return false;
		if (tmp == base) return false;//zero length name not allowed
		std::size_t checkpoint_pos = tmp-base+1;
		name = raw.substr(0, checkpoint_pos-1);
		int level = 1;
		while (level) {
			tmp = strpbrk(tmp+1, "(,)");
			if (!tmp) return false;
			else if (*tmp == '(') level++;
			else if (*tmp == ')' || *tmp == ',') {
				if (level == 1) {
					std::size_t new_checkpoint = tmp-base+1;
					if (new_checkpoint == checkpoint_pos + 1) return false;//parameter = empty string
					parameters.push_back(raw.substr(checkpoint_pos,new_checkpoint-checkpoint_pos-1));
					checkpoint_pos = new_checkpoint;
				}
				if (*tmp == ')') level--;
			}
		}
		return checkpoint_pos == raw.size();
	}
	PractRand::RNGs::vRNG *create_rng(const std::string &raw, std::string *error_message = nullptr) {
		std::string rng_name;
		std::vector<std::string> parameters;
		if (!parse_argument_list(raw, rng_name, parameters)) return nullptr;
		std::map<std::string, PractRand::RNGs::vRNG *(*)(std::vector<std::string> &params)>::iterator it = RNG_factory_index.find(rng_name);
		if (it == RNG_factory_index.end()) return nullptr;
		unsigned long os = parameters.size();
		PractRand::RNGs::vRNG *rng = it->second(parameters);
		if (error_message && !rng && os < parameters.size()) *error_message = parameters[os];
		return rng;
	}

	PractRand::RNGs::vRNG *rngset_lookup_recommended(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1 || index > num_recommended_rngs) { params.push_back("rngset lookup index out of range"); return nullptr; }
		return create_rng(recommended_rngs[index - 1]);
	}
	PractRand::RNGs::vRNG *rngset_lookup_nonrecommended_simple(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1 || index > num_nonrecommended_simple) { params.push_back("rngset lookup index out of range"); return nullptr; }
		return create_rng(nonrecommended_simple[index - 1]);
	}
	PractRand::RNGs::vRNG *rngset_lookup_nonrecommended_nonlcg(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1 || index > num_nonrecommended_nonlcg) { params.push_back("rngset lookup index out of range"); return nullptr; }
		return create_rng(nonrecommended_nonlcg[index - 1]);
	}
	PractRand::RNGs::vRNG *rngset_lookup_nonrecommended_lcgish(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1 || index > num_nonrecommended_lcgish) { params.push_back("rngset lookup index out of range"); return nullptr; }
		return create_rng(nonrecommended_lcgish[index - 1]);
	}
	PractRand::RNGs::vRNG *rngset_lookup_nonrecommended_cbuf(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1 || index > num_nonrecommended_cbuf) { params.push_back("rngset lookup index out of range"); return nullptr; }
		return create_rng(nonrecommended_cbuf[index - 1]);
	}
	PractRand::RNGs::vRNG *rngset_lookup_nonrecommended_indirect(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1 || index > num_nonrecommended_indirect) { params.push_back("rngset lookup index out of range"); return nullptr; }
		return create_rng(nonrecommended_indirect[index - 1]);
	}
	PractRand::RNGs::vRNG *rngset_lookup_nonrecommended(std::vector<std::string> &params) {
		if (params.size() < 1) { params.push_back("rngset lookup requires index"); return nullptr; }
		int index = atoi(params.front().c_str());
		using namespace PractRand::RNG_Sets;
		if (index < 1) {// || index > num_nonrecommended_simple + num_nonrecommended_nonlcg + num_nonrecommended_lcgish + num_nonrecommended_cbuf + num_nonrecommended_indirect) {
			params.push_back("rngset lookup index out of range");
			return nullptr;
		}
		if (index <= num_nonrecommended_simple) return create_rng(nonrecommended_simple[index - 1]);
		index -= num_nonrecommended_simple;
		if (index <= num_nonrecommended_nonlcg) return create_rng(nonrecommended_nonlcg[index - 1]);
		index -= num_nonrecommended_nonlcg;
		if (index <= num_nonrecommended_lcgish) return create_rng(nonrecommended_lcgish[index - 1]);
		index -= num_nonrecommended_lcgish;
		if (index <= num_nonrecommended_cbuf) return create_rng(nonrecommended_cbuf[index - 1]);
		index -= num_nonrecommended_cbuf;
		if (index <= num_nonrecommended_indirect) return create_rng(nonrecommended_indirect[index - 1]);
		params.push_back("rngset lookup index out of range");
		return nullptr;
	}
	template<class RNG>
	PractRand::RNGs::vRNG *_generic_recommended_RNG_factory(std::vector<std::string> &params) {
		if (params.size() != 0) return nullptr;
		return new RNG(PractRand::SEED_NONE);
	}
	template<class RNG>
	PractRand::RNGs::vRNG *_generic_notrecommended_RNG_factory(std::vector<std::string> &params) {
		if (params.size() != 0) return nullptr;
		return new RNG();
	}
	PractRand::RNGs::vRNG *lcg_factory(std::vector<std::string> &params) {
		if (params.size() != 2) {params.push_back("wrong number of parameters to lcg - should be lcg(out_bits,total_bits)");return nullptr;}
		int out_bits = atoi(params[0].c_str());
		int total_bits = atoi(params[1].c_str());
		if (out_bits!=8 && out_bits != 16 && out_bits != 32) {params.push_back("lcg out_bits must be 8, 16, or 32 bits");return nullptr;}
		if (total_bits < out_bits || total_bits > 128) {params.push_back("lcg total_bits invalid: must be out_bits <= total_bits <= 128");return nullptr;}
		if (out_bits == 8) {
			if (total_bits <= 64)
				return  new PractRand::RNGs::Polymorphic::NotRecommended::lcg8of64_varqual(total_bits - out_bits);
			else return new PractRand::RNGs::Polymorphic::NotRecommended::lcg8of128_varqual(total_bits - out_bits);
		}
		else if (out_bits == 16) {
			if (total_bits <= 64)
				return  new PractRand::RNGs::Polymorphic::NotRecommended::lcg16of64_varqual( total_bits - out_bits);
			else return new PractRand::RNGs::Polymorphic::NotRecommended::lcg16of128_varqual(total_bits - out_bits);
		}
		else /*if (out_bits == 32)*/ {
			if (total_bits <= 64)
				return  new PractRand::RNGs::Polymorphic::NotRecommended::lcg32of64_varqual( total_bits - out_bits);
			else return new PractRand::RNGs::Polymorphic::NotRecommended::lcg32of128_varqual(total_bits - out_bits);
		}
	}
	PractRand::RNGs::vRNG *xlcg_factory(std::vector<std::string> &params) {
		if (params.size() != 2) {params.push_back("wrong number of parameters to xlcg - should be xlcg(out_bits,total_bits)");return nullptr;}
		int out_bits = atoi(params[0].c_str());
		int total_bits = atoi(params[1].c_str());
		if (out_bits!=8 && out_bits != 16 && out_bits != 32) {params.push_back("xlcg out_bits must be 8, 16, or 32 bits");return nullptr;}
		if (total_bits < out_bits || total_bits > 128) {params.push_back("xlcg total_bits invalid: must be out_bits <= total_bits <= 64");return nullptr;}
		if (out_bits == 8) {
			if (total_bits <= 64)
				return new PractRand::RNGs::Polymorphic::NotRecommended::xlcg8of64_varqual(total_bits - out_bits);
			else return new PractRand::RNGs::Polymorphic::NotRecommended::xlcg8of128_varqual(total_bits - out_bits);
		}
		else if (out_bits == 16) {
			if (total_bits <= 64)
				return new PractRand::RNGs::Polymorphic::NotRecommended::xlcg16of64_varqual(total_bits - out_bits);
			else return new PractRand::RNGs::Polymorphic::NotRecommended::xlcg16of128_varqual(total_bits - out_bits);
		}
		else if (out_bits == 32) {
			if (total_bits <= 64)
				return new PractRand::RNGs::Polymorphic::NotRecommended::xlcg32of64_varqual(total_bits - out_bits);
			else return new PractRand::RNGs::Polymorphic::NotRecommended::xlcg32of128_varqual(total_bits - out_bits);
		}
		else return nullptr;
	}
	PractRand::RNGs::vRNG *clcg_factory(std::vector<std::string> &params) {
		if (params.size() != 2) {params.push_back("wrong number of parameters to clcg - should be clcg(out_bits,total_bits)");return nullptr;}
		int out_bits = atoi(params[0].c_str());
		int total_bits = atoi(params[1].c_str());
		if (out_bits!=8 && out_bits != 16 && out_bits != 32) {params.push_back("clcg out_bits must be 8, 16, or 32 bits");return nullptr;}
		if (total_bits < out_bits || total_bits > 96) {params.push_back("clcg total_bits invalid: must be out_bits+32 <= total_bits <= 96");return nullptr;}
		if (out_bits == 8) return new PractRand::RNGs::Polymorphic::NotRecommended::clcg8of96_varqual(total_bits - out_bits - 32);
		else if (out_bits == 16) return new PractRand::RNGs::Polymorphic::NotRecommended::clcg16of96_varqual(total_bits - out_bits - 32);
		else if (out_bits == 32) return new PractRand::RNGs::Polymorphic::NotRecommended::clcg32of96_varqual(total_bits - out_bits - 32);
		else return nullptr;
	}
	PractRand::RNGs::vRNG *cxlcg_factory(std::vector<std::string> &params) {
		if (params.size() != 2) {params.push_back("wrong number of parameters to clcg - should be clcg(out_bits,total_bits)");return nullptr;}
		int out_bits = atoi(params[0].c_str());
		int total_bits = atoi(params[1].c_str());
		if (out_bits!=8 && out_bits != 16 && out_bits != 32) {params.push_back("cxlcg out_bits must be 8, 16, or 32 bits");return nullptr;}
		if (total_bits < out_bits || total_bits > 96) {params.push_back("cxlcg total_bits invalid: must be out_bits+32 <= total_bits <= 96");return nullptr;}
		if (out_bits == 8) return new PractRand::RNGs::Polymorphic::NotRecommended::cxlcg8of96_varqual(total_bits - out_bits - 32);
		else if (out_bits == 16) return new PractRand::RNGs::Polymorphic::NotRecommended::cxlcg16of96_varqual(total_bits - out_bits - 32);
		else if (out_bits == 32) return new PractRand::RNGs::Polymorphic::NotRecommended::cxlcg32of96_varqual(total_bits - out_bits - 32);
		else return nullptr;
	}
	PractRand::RNGs::vRNG *bigbadlcg_factory(std::vector<std::string> &params) {
		if (params.size() != 3) { params.push_back("wrong number of parameters to bigbadlcg - should be bigbadlcg(out_bits,total_bits,shift)"); return nullptr; }
		int out_bits = atoi(params[0].c_str());
		int total_bits = atoi(params[1].c_str());
		int shift = atoi(params[2].c_str());
		if (out_bits != 8 && out_bits != 16 && out_bits != 32 && out_bits != 64) { params.push_back("cxlcg out_bits must be 8, 16, 32, or 64 bits"); return nullptr; }
		if (total_bits < out_bits || total_bits > 1024) { params.push_back("bigbadlcg total_bits invalid: must be out_bits <= total_bits <= 1024"); return nullptr; }
		if (out_bits == 8) return new PractRand::RNGs::Polymorphic::NotRecommended::bigbadlcg8X(total_bits - out_bits, shift);
		else if (out_bits == 16) return new PractRand::RNGs::Polymorphic::NotRecommended::bigbadlcg16X(total_bits - out_bits, shift);
		else if (out_bits == 32) return new PractRand::RNGs::Polymorphic::NotRecommended::bigbadlcg32X(total_bits - out_bits, shift);
		else if (out_bits == 64) return new PractRand::RNGs::Polymorphic::NotRecommended::bigbadlcg64X(total_bits - out_bits, shift);
		else return nullptr;
	}
	PractRand::RNGs::vRNG *chacha_factory(std::vector<std::string> &params) {
		if (params.size() > 1) {params.push_back("too many parameters for chacha");return nullptr;}
		if (params.size() == 0) return new PractRand::RNGs::Polymorphic::chacha(PractRand::SEED_NONE);
		std::string &param = params[0];
		int value = 0;
		for (unsigned long i = 0; i < param.size(); i++) {
			if (param[i] < '0' || param[i] > '9') {params.push_back("parameter must be an integer between 1 and 255, inclusive");return nullptr;}
			value = value * 10 + (param[i] - '0');
			if (value > 255) {params.push_back("number of rounds too high for chacha");return nullptr;}
		}
		if (value < 1) {params.push_back("number of rounds too low for chacha");return nullptr;}
		auto *rng = new PractRand::RNGs::Polymorphic::chacha(PractRand::SEED_NONE);
		rng->set_rounds(value);
		return rng;
	}
	PractRand::RNGs::vRNG *salsa_factory(std::vector<std::string> &params) {
		if (params.size() > 1) {params.push_back("too many parameters for salsa");return nullptr;}
		if (params.size() == 0) return new PractRand::RNGs::Polymorphic::salsa(PractRand::SEED_NONE);
		std::string &param = params[0];
		int value = 0;
		for (unsigned long i = 0; i < param.size(); i++) {
			if (param[i] < '0' || param[i] > '9') {params.push_back("parameter must be an integer between 1 and 255, inclusive");return nullptr;}
			value = value * 10 + (param[i] - '0');
			if (value > 255) {params.push_back("number of rounds too high for salsa");return nullptr;}
		}
		if (value < 1) {params.push_back("number of rounds too low for salsa");return nullptr;}
		auto *rng = new PractRand::RNGs::Polymorphic::salsa(PractRand::SEED_NONE);
		rng->set_rounds(value);
		return rng;
	}
	PractRand::RNGs::vRNG *SelfShrink_factory(std::vector<std::string> &params) {
		//SShrink(BaseRNG)
		if (params.size() != 1) return nullptr;
		PractRand::RNGs::vRNG *rng = create_rng(params[0]);
		if (!rng) return nullptr;
		return PractRand::RNGs::Polymorphic::NotRecommended::apply_SelfShrinkTransform(rng);
	}
	PractRand::RNGs::vRNG *BDS_factory(std::vector<std::string> &params) {
		//BDS(BaseRNG,log2_of_table_size)
		if (params.size() < 2 || params.size() > 3) return nullptr;
		int L2 = std::atoi(params[1].c_str());
		if (L2 < 1 || L2 > 16) return nullptr;
		int shift = -1;
		if (params.size() == 3) shift = std::atoi(params[2].c_str());
		PractRand::RNGs::vRNG *rng = create_rng(params[0]);
		if (!rng) return nullptr;
		return PractRand::RNGs::Polymorphic::NotRecommended::apply_BaysDurhamShuffle(rng, L2, shift);
	}
	template<class RNG>
	PractRand::RNGs::vRNG *_generic_single_parameter_transform_RNG_factory(std::vector<std::string> &params) {
		if (params.size() != 1) return nullptr;
		PractRand::RNGs::vRNG *rng = create_rng(params[0]);
		if (!rng) return nullptr;
		return new RNG(rng);
	}
	template<class RNG>
	PractRand::RNGs::vRNG *_generic_variable_parameter_transform_RNG_factory(std::vector<std::string> &params) {
		if (params.size() < 1) return nullptr;
		std::vector<PractRand::RNGs::vRNG *> source_rngs;
		PractRand::RNGs::vRNG *rng;
		for (unsigned int i = 0; i < params.size(); i++) {
			rng = create_rng(params[i]);
			if (!rng) return nullptr;
			source_rngs.push_back(rng);
		}
		rng = new RNG(source_rngs);
		if (rng) return rng;
		for (unsigned int i = 0; i < params.size(); i++) delete source_rngs[i];
		return nullptr;
	}
	template<class RNG>
	PractRand::RNGs::vRNG *_generic_single_parameter_RNG_factory(std::vector<std::string> &params) {
		if (params.size() != 1) return nullptr;
		return new RNG(atoi(params[0].c_str()));
	}
	template<class RNG>
	PractRand::RNGs::vRNG *_generic_two_parameter_RNG_factory(std::vector<std::string> &params) {
		if (params.size() != 2) return nullptr;
		return new RNG(atoi(params[0].c_str()),atoi(params[1].c_str()));
	}
	void register_recommended_RNGs() {
#define REGISTER_RNG_0(RNG) RNG_Factories::RNG_factory_index[ #RNG ] = _generic_recommended_RNG_factory<PractRand::RNGs::Polymorphic:: RNG>;
//#define REGISTER_RNG_0(RNG,DESC) register_rng( #RNG , _generic_recommended_RNG_factory<PractRand::RNGs::Polymorphic:: RNG>, DESC);
	//{ name, full_name, factory, description, param_description, quality, speed, word, size, flags };
//		RNG_Description rngs[] = {
//			{
//		};
		RNG_Factories::RNG_factory_index["recommended_rng"] = rngset_lookup_recommended;
		REGISTER_RNG_0(xsm32)
		REGISTER_RNG_0( xsm64 )
		REGISTER_RNG_0( hc256 )
		REGISTER_RNG_0( isaac32x256 )
		REGISTER_RNG_0( isaac64x256 )
		REGISTER_RNG_0( efiix8x48 )
		REGISTER_RNG_0( efiix16x48 )
		REGISTER_RNG_0( efiix32x48 )
		REGISTER_RNG_0( efiix64x48 )
		REGISTER_RNG_0( mt19937 )
		REGISTER_RNG_0( jsf32 )
		REGISTER_RNG_0( jsf64 )
		REGISTER_RNG_0( sfc16 )
		REGISTER_RNG_0( sfc32 )
		REGISTER_RNG_0( sfc64 )
		REGISTER_RNG_0( rarns16 )
		REGISTER_RNG_0( rarns32 )
		REGISTER_RNG_0( rarns64 )
		REGISTER_RNG_0( arbee )
		REGISTER_RNG_0( trivium )
		REGISTER_RNG_0( sha2_based_pool )
		RNG_Factories::RNG_factory_index[ "chacha" ] = chacha_factory;
		RNG_Factories::RNG_factory_index[ "salsa" ] = salsa_factory;
#undef REGISTER_RNG_0
	}
	void register_input_RNGs() {
#define REGISTER_RNG_0(RNG) RNG_Factories::RNG_factory_index[ #RNG ] = _generic_notrecommended_RNG_factory<  Special_RNGs:: RNG_ ## RNG>;
		RNG_Factories::RNG_factory_index["stdin"] = _generic_notrecommended_RNG_factory<Special_RNGs::RNG_stdin>;
		REGISTER_RNG_0(stdin8)
		REGISTER_RNG_0(stdin16)
		REGISTER_RNG_0(stdin32)
		REGISTER_RNG_0(stdin64)
#undef REGISTER_RNG_0
	}
	void register_nonrecommended_RNGs() {
		RNG_Factories::RNG_factory_index["nonrecommended_rng_simple"] = rngset_lookup_nonrecommended_simple;
		RNG_Factories::RNG_factory_index["nonrecommended_rng_nonlcg"] = rngset_lookup_nonrecommended_nonlcg;
		RNG_Factories::RNG_factory_index["nonrecommended_rng_lcgish"] = rngset_lookup_nonrecommended_lcgish;
		RNG_Factories::RNG_factory_index["nonrecommended_rng_cbuf"] = rngset_lookup_nonrecommended_cbuf;
		RNG_Factories::RNG_factory_index["nonrecommended_rng_indirect"] = rngset_lookup_nonrecommended_indirect;
		RNG_Factories::RNG_factory_index["nonrecommended_rng"] = rngset_lookup_nonrecommended;
#define REGISTER_RNG_0(RNG) RNG_Factories::RNG_factory_index[ #RNG ] = _generic_notrecommended_RNG_factory<  PractRand::RNGs::Polymorphic::NotRecommended:: RNG>;
#define REGISTER_RNG_1(RNG) RNG_Factories::RNG_factory_index[ #RNG ] = _generic_single_parameter_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended:: RNG>;
#define REGISTER_RNG_2(RNG) RNG_Factories::RNG_factory_index[ #RNG ] = _generic_two_parameter_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended:: RNG>;
		// include/PractRand/RNGs/other/simple.h
		REGISTER_RNG_0(xsalta16x3)
		REGISTER_RNG_0(xsaltb16x3)
		REGISTER_RNG_0(xsaltc16x3)
		REGISTER_RNG_0(xorshift32)
		REGISTER_RNG_0(xorshift64)
		REGISTER_RNG_0(xorshift32of128)
		REGISTER_RNG_0(xorshift64of128)
		REGISTER_RNG_0(xorshift16of32)
		REGISTER_RNG_0(xorshift32of64)
		REGISTER_RNG_0(xorshift32x4)
		REGISTER_RNG_0(xorwow32plus32)
		REGISTER_RNG_0(xorwow32of96)
		REGISTER_RNG_0(xorwow32x6)
		REGISTER_RNG_0(xorshift128plus)
		REGISTER_RNG_0(xoroshiro128plus)
		REGISTER_RNG_0(xoroshiro128plus_2p64)
		REGISTER_RNG_0(tinyMT)
		REGISTER_RNG_0(sapparot)
		REGISTER_RNG_0(sap16of48)
		REGISTER_RNG_0(sap32of96)
		REGISTER_RNG_0(flea32x1)
		REGISTER_RNG_0(jsf16)
		REGISTER_RNG_0(tyche)
		REGISTER_RNG_0(tyche16)
		REGISTER_RNG_0(sfc_v1_16)
		REGISTER_RNG_0(sfc_v1_32)
		REGISTER_RNG_0(sfc_v2_16)
		REGISTER_RNG_0(sfc_v2_32)
		REGISTER_RNG_0(sfc_v3_16)
		REGISTER_RNG_0(sfc_v3_32)
		REGISTER_RNG_0(simpleA)
		REGISTER_RNG_0(simpleB)
		REGISTER_RNG_0(simpleC)
		REGISTER_RNG_0(simpleD)
		REGISTER_RNG_0(simpleE)
		REGISTER_RNG_0(simpleF)
		REGISTER_RNG_0(simpleG)
		REGISTER_RNG_0(trivium_weakenedA)
		REGISTER_RNG_0(trivium_weakenedB)
		REGISTER_RNG_0(mo_Lesr32)
		REGISTER_RNG_0(mo_ResrRers32)
		REGISTER_RNG_0(mo_Rers32of64)
		REGISTER_RNG_0(mo_Resr32of64)
		REGISTER_RNG_0(mo_Resdra32of64)
		REGISTER_RNG_0(murmlacish)
		REGISTER_RNG_0(gjishA)
		REGISTER_RNG_0(gjishB)
		REGISTER_RNG_0(gjishC)
		REGISTER_RNG_0(gjishD)
		REGISTER_RNG_0(ara16)
		REGISTER_RNG_0(ara32)
		REGISTER_RNG_0(arx16)
		REGISTER_RNG_0(arx32)
		REGISTER_RNG_0(hara16)
		REGISTER_RNG_0(harx16)
		REGISTER_RNG_0(learx16)
		REGISTER_RNG_0(hlearx16)
		REGISTER_RNG_0(alearx16)
		REGISTER_RNG_0(arac16)
		REGISTER_RNG_0(arxc16)

		// include/PractRand/RNGs/other/mult.h
		RNG_Factories::RNG_factory_index["lcg"] = lcg_factory;
		RNG_Factories::RNG_factory_index["xlcg"] = xlcg_factory;
		RNG_Factories::RNG_factory_index["clcg"] = clcg_factory;
		RNG_Factories::RNG_factory_index["cxlcg"] = cxlcg_factory;
		REGISTER_RNG_0(lcg16of32_extended)
		REGISTER_RNG_0(lcg32_extended)
		REGISTER_RNG_0(pcg32)
		REGISTER_RNG_0(pcg32_norot)
		REGISTER_RNG_0(cmrg32of192)
		REGISTER_RNG_0(xsh_lcg_bad)
		REGISTER_RNG_0(mmr16)
		REGISTER_RNG_0(mmr32)
		REGISTER_RNG_0(garthy16)
		REGISTER_RNG_0(garthy32)
		REGISTER_RNG_0(binarymult16)
		REGISTER_RNG_0(binarymult32)
		REGISTER_RNG_0(rxmult16)
		REGISTER_RNG_0(multish2x64)
		REGISTER_RNG_0(multish3x32)
		REGISTER_RNG_0(multish4x16)
		REGISTER_RNG_0(mwrca16)
		REGISTER_RNG_0(mwrca32)
		REGISTER_RNG_0(mwrcc16)
		REGISTER_RNG_0(mwrcc32)
		REGISTER_RNG_0(mwrcca16)
		REGISTER_RNG_0(mwrcca32)
		REGISTER_RNG_0(old_mwlac16)
		REGISTER_RNG_0(mwlac_varA)
		REGISTER_RNG_0(mwlac_varB)
		REGISTER_RNG_0(mwlac_varC)
		REGISTER_RNG_0(mwlac_varD)
		REGISTER_RNG_0(mwlac_varE)
		REGISTER_RNG_0(mwc64x)
		REGISTER_RNG_0(mo_Cmfr32)
		REGISTER_RNG_0(mo_Cmr32)
		REGISTER_RNG_0(mo_Cmr32of64)
		REGISTER_RNG_1(murmlac32)
		REGISTER_RNG_0(mulcr16)
		REGISTER_RNG_0(mulcr32)
		REGISTER_RNG_0(mulcr64)

		// include/PractRand/RNGs/other/fibonacci.h
		RNG_Factories::RNG_factory_index["bigbadlcg"] = bigbadlcg_factory;
		RNG_Factories::RNG_factory_index["bblcg"] = bigbadlcg_factory;
		REGISTER_RNG_0(lfsr_medium)
		REGISTER_RNG_0(mm32)
		REGISTER_RNG_0(mm16of32)
		REGISTER_RNG_0(mm32_awc)
		REGISTER_RNG_0(mm16of32_awc)
		REGISTER_RNG_0(mwc4691)
		REGISTER_RNG_0(cbuf_accum)
		REGISTER_RNG_0(cbuf_accum_big)
		REGISTER_RNG_0(cbuf_2accum_small)
		REGISTER_RNG_0(cbuf_2accum)
		REGISTER_RNG_0(dual_cbuf_small)
		REGISTER_RNG_0(dual_cbuf)
		REGISTER_RNG_0(dual_cbufa_small)
		REGISTER_RNG_0(dual_cbuf_accum)
		REGISTER_RNG_0(fibmul16of32)
		REGISTER_RNG_0(fibmul32of64)
		REGISTER_RNG_0(fibmulmix16)
		REGISTER_RNG_0(ranrot32small)
		REGISTER_RNG_0(ranrot32)
		REGISTER_RNG_0(ranrot32big)
		REGISTER_RNG_0(ranrot3tap32small)
		REGISTER_RNG_0(ranrot3tap32)
		REGISTER_RNG_0(ranrot3tap32big)
		REGISTER_RNG_0(ranrot32hetsmall)
		REGISTER_RNG_0(ranrot32het)
		REGISTER_RNG_0(ranrot32hetbig)
		REGISTER_RNG_0(mt19937_unhashed)

		// include/PractRand/RNGs/other/indirection.h
		REGISTER_RNG_0(rc4)
		REGISTER_RNG_0(rc4_weakenedA)
		REGISTER_RNG_0(rc4_weakenedB)
		REGISTER_RNG_0(rc4_weakenedC)
		REGISTER_RNG_0(rc4_weakenedD)
		REGISTER_RNG_1(ibaa8)
		REGISTER_RNG_1(ibaa16)
		REGISTER_RNG_1(ibaa32)
		REGISTER_RNG_1(isaac32_varqual)
		REGISTER_RNG_1(isaac16_varqual)
		REGISTER_RNG_2(efiix4_varqual)
		REGISTER_RNG_2(efiix8_varqual)
		REGISTER_RNG_1(genindA)
		REGISTER_RNG_1(genindB)
		REGISTER_RNG_1(genindC)
		REGISTER_RNG_1(genindD)
		REGISTER_RNG_1(genindE)
		REGISTER_RNG_1(genindF)

		// include/PractRand/RNGs/other/transform.h
		RNG_Factories::RNG_factory_index["xor"] = _generic_variable_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Xor>;
		RNG_Factories::RNG_factory_index["BDS"] = BDS_factory;
		RNG_Factories::RNG_factory_index["SShrink"] = SelfShrink_factory;
		RNG_Factories::RNG_factory_index["AsUnknown"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::ReinterpretAsUnknown>;
		RNG_Factories::RNG_factory_index["As8" ] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::ReinterpretAs8>;
		RNG_Factories::RNG_factory_index["As16"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::ReinterpretAs16>;
		RNG_Factories::RNG_factory_index["As32"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::ReinterpretAs32>;
		RNG_Factories::RNG_factory_index["As64"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::ReinterpretAs64>;
		RNG_Factories::RNG_factory_index["Discard16to8"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Discard16to8>;
		RNG_Factories::RNG_factory_index["Discard32to8"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Discard32to8>;
		RNG_Factories::RNG_factory_index["Discard64to8"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Discard64to8>;
		RNG_Factories::RNG_factory_index["Discard32to16"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Discard32to16>;
		RNG_Factories::RNG_factory_index["Discard64to16"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Discard64to16>;
		RNG_Factories::RNG_factory_index["Discard64to32"] = _generic_single_parameter_transform_RNG_factory<PractRand::RNGs::Polymorphic::NotRecommended::Discard64to32>;
#undef REGISTER_RNG_0
#undef REGISTER_RNG_1
#undef REGISTER_RNG_2
	}
}
