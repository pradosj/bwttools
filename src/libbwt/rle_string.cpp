
#include "rle_string.h"
#include <iostream>
#include <istream>
#include <string>
#include <numeric>
#include <exception>


namespace bwt {

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


std::istream& operator>>(std::istream& is,rle_string& str) {
		BWFlag flag;
    size_t numStrings, numSymbols, numRuns;
    readHeader(is, numStrings, numSymbols, numRuns, flag);
    assert(numRuns > 0);
    
    str.runs.clear();
		str.runs.resize(numRuns);
		is.read(reinterpret_cast<char*>(&str.runs[0]), numRuns*sizeof(rle_unit));
		
		str.update_size();
#ifndef NDEBUG
		std::cerr << "An rle_string have been read:" << std::endl;
		std::cerr << "  size:" << numSymbols << std::endl;
		std::cerr << "  #run:" << str.runs.size() << std::endl;
		std::cerr << "  avg run size:" << ((double) numSymbols) / str.runs.size() << std::endl;
#endif		
		if (numSymbols!=str.size()) throw "BWT file corrupted: the number of symbol provided in the header do not match with the effective number of symbol";
		return is;
}


std::ostream& rle_string::print(std::ostream& os) {
		for(auto r:runs) os << r.length() << 'x' << r.value() << ' ';
		return os;
}



};