//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "BWTAlgorithms.h"

// Find the interval in pBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
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
        if(!interval.isValid())
            return interval;
    }
    return interval;
}



// Find the intervals in pBWT/pRevBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
BWTIntervalPair BWTAlgorithms::findIntervalPair(const bwt* pBWT, const bwt* pRevBWT, const std::string& w) {
    BWTIntervalPair intervals;    
    int len = w.size();
    int j = len - 1;
    char curr = w[j];
    initIntervalPair(intervals, curr, pBWT, pRevBWT);
    --j;

    for(;j >= 0; --j)
    {
        curr = w[j];
        updateBothL(intervals, curr, pBWT);
        if(!intervals.isValid())
            return intervals;
    }
    return intervals;
}



// Count the number of occurrences of string w
size_t BWTAlgorithms::countSequenceOccurrences(const std::string& w, const bwt* pBWT) {
    BWTInterval fwd_interval = findInterval(pBWT, w);
    return (fwd_interval.isValid())?fwd_interval.size():0;
}



// Return the count of all the possible one base extensions of the string w.
// This returns the number of times the suffix w[i, l]A, w[i, l]C, etc 
// appears in the FM-index for all i s.t. length(w[i, l]) == overlapLen.
AlphaCount64 BWTAlgorithms::calculateExactExtensions(const unsigned int overlapLen, const std::string& w, const bwt* pBWT, const bwt* pRevBWT)
{
    // The algorithm is as follows:
    // We perform a backward search on the FM-index of w.
    // For each signficant suffix (length w[i,l] >= minOverlap)
    // we determine the proper prefixes that match w[i,l]. For each proper prefix matching, 
    // we compute the number of extensions of A,C,G,T for those prefix.
    AlphaCount64 ext_counts;
    BWTIntervalPair ranges;
    size_t l = w.length();
    int start = l - 1;
    BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);

    for(int i = start - 1; i >= 0; --i)
    {
        // Compute the range of the suffix w[i, l]
        BWTAlgorithms::updateBothL(ranges, w[i], pBWT);

        // Break if the suffix is no longer found
        if(!(ranges.interval[0].isValid() && ranges.interval[1].isValid())) 
            break;

        if((l - i) == overlapLen)
        {
            if(ranges.interval[1].isValid())
            {
                assert(ranges.interval[1].lower > 0);
                // The count for each extension is the difference between rank(B, upper) and rank(B, lower - 1)
                AlphaCount64 ac = pRevBWT->getOccDiff(ranges.interval[1].lower - 1, ranges.interval[1].upper);
                ext_counts += ac;
            }
        }
    }
    return ext_counts;
}



// Return the string from the BWT at idx
std::string BWTAlgorithms::extractString(const bwt* pBWT, size_t idx) {
    assert(idx < pBWT->getNumStrings());

    // The range [0,n) in the BWT contains all the terminal
    // symbols for the reads. Search backwards from one of them
    // until the '$' is found gives a full string.
    std::string out;
    BWTInterval interval(idx, idx);
    while(1)
    {
        assert(interval.isValid());
        char b = pBWT->getChar(interval.lower);
        if(b == '$')
            break;
        else
            out.push_back(b);
        updateInterval(interval, b, pBWT);
    } 
    return reverse(out);
}

// Extract the substring from start, start+length of the sequence starting at position idx
std::string BWTAlgorithms::extractSubstring(const bwt* pBWT, uint64_t idx, size_t start, size_t length) {
    std::string s = extractString(pBWT, idx);
    return s.substr(start, length);
}

// Return the next len bases of the string starting at index idx of the BWT
std::string BWTAlgorithms::extractString(const bwt* pBWT, size_t idx, size_t len) {
    std::string out;
    BWTInterval interval(idx, idx);
    while(out.length() < len)
    {
        assert(interval.isValid());
        char b = pBWT->getChar(interval.lower);
        if(b == '$')
            break;
        else
            out.push_back(b);
        updateInterval(interval, b, pBWT);
    } 
    return reverse(out);
}


// Recursive traversal to extract all the strings needed for the above function
void _extractRankedPrefixes(const bwt* pBWT, BWTInterval interval, const std::string& curr, RankedPrefixVector* pOutput)
{
    AlphaCount64 extensions = BWTAlgorithms::getExtCount(interval, pBWT);

    for(size_t i = 0; i < 4; ++i)
    {
        char b = "ACGT"[i];

        if(extensions.get(b) > 0)
        {
            BWTInterval ni = interval;
            BWTAlgorithms::updateInterval(ni, b, pBWT);
            _extractRankedPrefixes(pBWT, ni, curr + b, pOutput);
        }

    }

    // If we have extended the prefix as far as possible, stop
    BWTAlgorithms::updateInterval(interval, '$', pBWT);
    for(int64_t i = interval.lower; i <= interval.upper; ++i)
    {
        // backwards search gives a reversed prefix, fix it
        RankedPrefix rp = { (size_t)i, reverse(curr) };
        pOutput->push_back(rp);
    }
}

// Extract all strings found from a backwards search starting at the given interval
RankedPrefixVector BWTAlgorithms::extractRankedPrefixes(const bwt* pBWT, BWTInterval interval) {
    std::string curr;
    RankedPrefixVector output;
    output.reserve(interval.size());
    _extractRankedPrefixes(pBWT, interval, curr, &output);
    return output;
}

std::string BWTAlgorithms::extractUntilInterval(const bwt* pBWT, int64_t start, const BWTInterval& check) {
    std::string out;
    BWTInterval interval(start, start);
    while(interval.lower < check.lower || interval.lower > check.upper)
    {
        assert(interval.isValid());
        char b = pBWT->getChar(interval.lower);
        if(b == '$')
            return "";
        else
            out.push_back(b);
        updateInterval(interval, b, pBWT);
    } 
    return reverse(out);        
}
