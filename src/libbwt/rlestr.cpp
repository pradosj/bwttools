
#include "rlestr.h"
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>


enum BWFlag {BWF_NOFMI = 0,BWF_HASFMI};

static void readHeader(std::istream& is,size_t& num_strings, size_t& num_symbols, size_t& num_runs, BWFlag& flag) {
    uint16_t magic_number;
    is.read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
    if(magic_number != 0xCACA) {
        std::cerr << "BWT file is not properly formatted, aborting\n";
        exit(EXIT_FAILURE);
    }
    is.read(reinterpret_cast<char*>(&num_strings), sizeof(num_strings));
    is.read(reinterpret_cast<char*>(&num_symbols), sizeof(num_symbols));
    is.read(reinterpret_cast<char*>(&num_runs), sizeof(num_runs));
    is.read(reinterpret_cast<char*>(&flag), sizeof(flag));
}


void RLEString::read(const std::string& filename) {
		std::ifstream is(filename, std::ios::binary);
		BWFlag flag;
    size_t numStrings, numRuns;
    readHeader(is, numStrings, m_numSymbols, numRuns, flag);
    assert(numRuns > 0);
    
		resize(numRuns);
		is.read(reinterpret_cast<char*>(&(*this)[0]), numRuns*sizeof(RLUnit));
		
#ifndef NDEBUG
		std::cerr << "An RLEString have been read:" << size() << std::endl;
		std::cerr << "Num run:" << size() << std::endl;
		std::cerr << "Num symbol:" << m_numSymbols << std::endl;
		std::cerr << "Sum of lengths:" << std::accumulate(begin(),end(),0,[](uint64_t s,RLUnit u) -> uint64_t{return s+u.length();}) << std::endl;
#endif
		assert(m_numSymbols == std::accumulate(begin(),end(),0,[](uint64_t s,RLUnit u){return s+u.length();}));
}

