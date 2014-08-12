#ifndef ALPHABET_H
#define ALPHABET_H

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <string>


typedef uint64_t BaseCount;

//const std::string BWT_ALPHABET("$ACGT");

namespace BWT_ALPHABET {
		const uint8_t ALPHABET_SIZE = 5;
		const char RANK_ALPHABET[ALPHABET_SIZE] = {'$', 'A', 'C', 'G', 'T'};
    static const uint8_t s_bwtLexoRankLUT[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
        0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    inline static uint8_t getRank(char b) {return s_bwtLexoRankLUT[static_cast<uint8_t>(b)];}
    inline char getChar(size_t idx) {return RANK_ALPHABET[idx];}
};


inline std::string reverse(const std::string& str) {
	std::string rstr(str);
  std::reverse(rstr.begin(),rstr.end());
  return(rstr);
}

inline std::string complement(std::string str) {
    auto complement_bp = [](char bp) -> char {
    	switch(bp) {
        case 'A':return 'T';
        case 'C':return 'G';
        case 'G':return 'C';
        case 'T':return 'A';
        case 'N':return 'N';
        default:
            assert(false && "Unknown base!");
            return 'N';
			};
    };
    std::transform(str.begin(), str.end(), str.begin(), complement_bp);
    return(str);
}

#endif
