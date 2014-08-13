//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "BWTAlgorithms.h"


BWTInterval BWTAlgorithms::findInterval(const bwt* pBWT, const std::string& w) {
    int len = w.size();
    int j = len - 1;
    char curr = w[j];
    BWTInterval interval = pBWT->initInterval(curr);
    for(--j;j >= 0; --j) {
        curr = w[j];
        pBWT->updateInterval(interval, curr);
        if(!interval.isValid()) return interval;
    }
    return interval;
}



// Return the string from the BWT at idx
std::string BWTAlgorithms::extractString(const bwt& BWT, size_t idx) {
    // The range [0,n) in the BWT contains all the terminal
    // symbols for the reads. Search backwards from one of them
    // until the '$' is found gives a full string.
    std::string out;
    BWTInterval interval(idx, idx);
    while(1) {
        assert(interval.isValid());
        uint8_t b = BWT.symbol(interval.lower);
        if(b == '$') break;
        out.push_back(b);
        BWT.updateInterval(interval, b);
    } 
    return reverse(out);
}



