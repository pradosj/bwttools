#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include <string>
#include "bwt.h"


// functions
namespace BWTAlgorithms {

// get the interval(s) in pBWT that corresponds to the string w using a backward search algorithm
// Find the interval in pBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
bwt::interval findInterval(const bwt::fm_index& fm, const std::string& w) {
		if (w.size()<1) return bwt::interval();
    bwt::interval range = fm.initInterval(w.back());
    for(auto i=w.rbegin()+1;i!=w.rend();i++) {
        fm.updateInterval(range,*i);
        if (range.empty()) return range;
    }
    return range;
}

// Return the string from the BWT at idx
std::string extractString(const bwt::fm_index& fm, size_t idx) {
    // The range [0,n) in the BWT contains all the terminal
    // symbols for the reads. Search backwards from one of them
    // until the '$' is found gives a full string.
    std::string out;
    bwt::interval range(idx, idx);
    while(1) {
        assert(!range.empty());
        uint8_t b = fm.symbol(range.lower);
        if (b == 0) break;
        out.push_back(b);
        fm.updateInterval(range, b);
    }
    std::reverse(out.begin(),out.end());
    return out;
}


};

#endif
