#ifndef BWTALGO_H
#define BWTALGO_H

#include "fm_index.h"
#include <inttypes.h>	


/*! \class interval
 *  \brief An interval holds a range in a suffix array
 */
struct interval {
    uint64_t first;
    uint64_t last;
    interval() : first(0), last(0) {}
    interval(uint64_t l, uint64_t u) : first(l), last(u) {}
    inline bool empty() const {return last <= first;}
    inline uint64_t size() const {return empty()?0:last-first;}
    friend inline bool operator==(const interval& a, const interval&b) {return a.first==b.first && a.last==b.last;}
};



//! return the full suffix array interval
template <size_t AlphabetSize>
inline interval sa_interval(const fm_index<AlphabetSize>& fm) {
		return interval(0,fm.size());
}

//! return the suffix array interval for character b
template <size_t AlphabetSize>
inline interval sa_interval(const fm_index<AlphabetSize>& fm, const uint8_t b) {
		return interval(fm.C(b),fm.C(b>=fm.alphabet_size()?fm.size():b+1));
}

//! \brief update a suffix array interval using backwards search
//! if the given interval corresponds to string S, it will be updated for string bS
template <size_t AlphabetSize>
inline void update_sa_interval(interval& interval, const fm_index<AlphabetSize>& fm, const uint8_t b) {
		//NOTE: computation of occ(.,last) might be faster using occ(.,first)
    interval.first = fm.C(b) + (interval.first>0?fm.occ(b,interval.first-1):0);
    interval.last = fm.C(b) + fm.occ(b,interval.last-1);
}



#endif