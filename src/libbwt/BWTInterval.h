#ifndef BWTINTERVAL_H
#define BWTINTERVAL_H

#include <inttypes.h>

/*! \class BWTInterval
 *  \brief Data structures for holding and manipulating the coordinates in a BWT/FM-index
 * 
 *  A BWTInterval holds a pair of integers which delineate an alignment of some string to a BWT/Suffix Array
 */
struct BWTInterval {
    // Functions
    BWTInterval() : lower(0), upper(0) {}
    BWTInterval(uint64_t l, uint64_t u) : lower(l), upper(u) {}

    inline bool isValid() const {return lower <= upper;}
    inline uint64_t size() const {return isValid()?upper-lower+1:0;}

    // Data
    uint64_t lower;
    uint64_t upper;
};


#endif

