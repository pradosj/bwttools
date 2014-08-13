#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include "bwt.h"
#include "BWTInterval.h"


// functions
namespace BWTAlgorithms {

// get the interval(s) in pBWT that corresponds to the string w using a backward search algorithm
// Find the interval in pBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
BWTInterval findInterval(const bwt* pBWT, const std::string& w);

// Extract the complete string starting at idx in the BWT
std::string extractString(const bwt& BWT, size_t idx);


};

#endif
