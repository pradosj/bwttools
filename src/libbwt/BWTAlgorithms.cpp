//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "BWTAlgorithms.h"


BWTInterval BWTAlgorithms::findInterval(const bwt* pBWT, const std::string& w) {
    int len = w.size();
    int j = len - 1;
    char curr = w[j];
    BWTInterval interval;
    initInterval(interval, curr, pBWT);
    --j;

    for(;j >= 0; --j) {
        curr = w[j];
        updateInterval(interval, curr, pBWT);
        if(!interval.isValid()) return interval;
    }
    return interval;
}



// Return the string from the BWT at idx
std::string BWTAlgorithms::extractString(const bwt& BWT, size_t idx) {
    assert(idx < BWT.getNumStrings());

    // The range [0,n) in the BWT contains all the terminal
    // symbols for the reads. Search backwards from one of them
    // until the '$' is found gives a full string.
    std::string out;
    BWTInterval interval(idx, idx);
    while(1) {
        assert(interval.isValid());
        char b = BWT[interval.lower];
        if(b == '$')
            break;
        else
            out.push_back(b);
        updateInterval(interval, b, &BWT);
    } 
    return reverse(out);
}


// Recursive traversal to extract all the strings needed for the above function
void _extractRankedPrefixes(const bwt* pBWT, BWTInterval interval, const std::string& curr, RankedPrefixVector* pOutput) {
    AlphaCount64 extensions = BWTAlgorithms::getExtCount(interval, pBWT);

    for(size_t i = 0; i < 4; ++i) {
        char b = "ACGT"[i];

        if(extensions[b] > 0) {
            BWTInterval ni = interval;
            BWTAlgorithms::updateInterval(ni, b, pBWT);
            _extractRankedPrefixes(pBWT, ni, curr + b, pOutput);
        }

    }

    // If we have extended the prefix as far as possible, stop
    BWTAlgorithms::updateInterval(interval, '$', pBWT);
    for(int64_t i = interval.lower; i <= interval.upper; ++i) {
        // backwards search gives a reversed prefix, fix it
        RankedPrefix rp = { (size_t)i, reverse(curr) };
        pOutput->push_back(rp);
    }
}

// Extract all strings found from a backwards search starting at the given interval
RankedPrefixVector BWTAlgorithms::extractRankedPrefixes(const bwt& BWT, BWTInterval interval) {
    std::string curr;
    RankedPrefixVector output;
    output.reserve(interval.size());
    _extractRankedPrefixes(&BWT, interval, curr, &output);
    return output;
}

std::string BWTAlgorithms::extractUntilInterval(const bwt& BWT, int64_t start, const BWTInterval& check) {
    std::string out;
    BWTInterval interval(start, start);
    while(interval.lower < check.lower || interval.lower > check.upper) {
        assert(interval.isValid());
        char b = BWT[interval.lower];
        if(b == '$') return "";
        out.push_back(b);
        updateInterval(interval, b, &BWT);
    } 
    return reverse(out);        
}
