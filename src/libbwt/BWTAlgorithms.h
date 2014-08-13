#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include <string>
#include "bwt.h"
#include "BWTInterval.h"


// functions
namespace BWTAlgorithms {

// get the interval(s) in pBWT that corresponds to the string w using a backward search algorithm
// Find the interval in pBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
BWTInterval findInterval(const bwt& BWT, const std::string& w) {
		if (w.size()<1) return BWTInterval();
    BWTInterval interval = BWT.initInterval(w.back());
    for(auto i=w.rbegin()+1;i!=w.rend();i++) {
        BWT.updateInterval(interval,*i);
        if (!interval.isValid()) return interval;
    }
    return interval;
}

// Return the string from the BWT at idx
std::string extractString(const bwt& BWT, size_t idx) {
    // The range [0,n) in the BWT contains all the terminal
    // symbols for the reads. Search backwards from one of them
    // until the '$' is found gives a full string.
    std::string out;
    BWTInterval interval(idx, idx);
    while(1) {
        assert(interval.isValid());
        uint8_t b = BWT.symbol(interval.lower);
        if (b == '$') break;
        out.push_back(b);
        BWT.updateInterval(interval, b);
    }
    std::reverse(out.begin(),out.end());
    return out;
}


};

#endif
