//
// BWTAlgorithms.h - Algorithms for aligning to a bwt structure
//
#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include "bwt.h"
#include "BWTInterval.h"
#include <queue>
#include <list>


// structures

// A (partial) prefix of a string contained in the BWT
// and its lexicographic rank
struct RankedPrefix {
    size_t rank;
    std::string prefix;
};
typedef std::vector<RankedPrefix> RankedPrefixVector;

// functions
namespace BWTAlgorithms {

// Update the given interval using backwards search
// If the interval corresponds to string S, it will be updated for string bS
inline void updateInterval(BWTInterval& interval, char b, const bwt* pBWT) {
    size_t pb = pBWT->getPC(b);
    interval.lower = pb + pBWT->getFullOcc(interval.lower - 1)[b];
    interval.upper = pb + pBWT->getFullOcc(interval.upper)[b] - 1;
}

// get the interval(s) in pBWT that corresponds to the string w using a backward search algorithm
// Find the interval in pBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
BWTInterval findInterval(const bwt* pBWT, const std::string& w);


// Initialize the interval of index idx to be the range containining all the b suffixes
inline void initInterval(BWTInterval& interval, char b, const bwt* pB) {
    interval.lower = pB->getPC(b);
    interval.upper = interval.lower + pB->getFullOcc(pB->getBWLen() - 1)[b] - 1;
}

// Return the counts of the bases between the lower and upper interval in pBWT
inline AlphaCount64 getExtCount(const BWTInterval& interval, const bwt* pBWT) {
    return pBWT->getOccDiff(interval.lower - 1, interval.upper);
}

// Extract the complete string starting at idx in the BWT
std::string extractString(const bwt& BWT, size_t idx);

// Extract the substring from start, start+length of the sequence starting at position idx
std::string extractSubstring(const bwt& BWT, uint64_t idx, size_t start, size_t length = std::string::npos);

// Extract all prefixes of the suffixes for the given interval, along
// with their lexicographic rank.
RankedPrefixVector extractRankedPrefixes(const bwt& BWT, BWTInterval interval);

// Extract symbols from the starting index of the BWT until the index lies
// within the given interval. If the extraction hits the start of a string
// without finding a prefix, the empty string is returned.
std::string extractUntilInterval(const bwt& BWT, int64_t start, const BWTInterval& interval);

};

#endif
