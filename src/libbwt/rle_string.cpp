
#include "rle_string.h"
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <exception>


enum BWFlag {BWF_NOFMI = 0,BWF_HASFMI};

static void readHeader(std::istream& is,size_t& num_strings, size_t& num_symbols, size_t& num_runs, BWFlag& flag) {
    uint16_t magic_number;
    is.read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
    if (magic_number != 0xCACA) throw "BWT file is not properly formatted: the magic number provided in file header doesn't correspond to the expected one";
    is.read(reinterpret_cast<char*>(&num_strings), sizeof(num_strings));
    is.read(reinterpret_cast<char*>(&num_symbols), sizeof(num_symbols));
    is.read(reinterpret_cast<char*>(&num_runs), sizeof(num_runs));
    is.read(reinterpret_cast<char*>(&flag), sizeof(flag));
}


void rle_string::read(const std::string& filename) {
		std::ifstream is(filename, std::ios::binary);
		BWFlag flag;
    size_t numStrings, numSymbols, numRuns;
    readHeader(is, numStrings, numSymbols, numRuns, flag);
    assert(numRuns > 0);
    
		resize(numRuns);
		is.read(reinterpret_cast<char*>(&(*this)[0]), numRuns*sizeof(rle_unit));
		
		m_numSymbols = std::accumulate(begin(),end(),0,[](uint64_t s,rle_unit u){return s+u.length();});
		if (numSymbols!=m_numSymbols) throw "BWT file corrupted: the number of symbol provided in the header do not match with the effective number of symbol";
		
#ifndef NDEBUG
		std::cerr << "An RLEString have been read:" << size() << std::endl;
		std::cerr << "Num run:" << size() << std::endl;
		std::cerr << "Num symbol:" << m_numSymbols << std::endl;
		std::cerr << "Sum of lengths:" << std::accumulate(begin(),end(),0,[](uint64_t s,rle_unit u) -> uint64_t{return s+u.length();}) << std::endl;
#endif
		
}

