#ifndef INTERVAL_H
#define INTERVAL_H

#include <inttypes.h>


namespace bwt {

/*! \class interval
 *  \brief Data structures for holding and manipulating the coordinates in a BWT/FM-index
 * 
 *  An interval holds a pair of integers which delineate an alignment of some string to a BWT/Suffix Array
 */
struct interval {
    uint64_t lower;
    uint64_t upper;

    interval() : lower(0), upper(0) {}
    interval(uint64_t l, uint64_t u) : lower(l), upper(u) {}

    inline bool empty() const {return lower > upper;}
    inline uint64_t size() const {return empty()?0:upper-lower+1;}	
};

};	


#endif