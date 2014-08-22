#ifndef BWTALGO_H
#define BWTALGO_H

#include "fm_index.h"
#include <cinttypes>	


namespace bwt {


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



/*! \brief update a suffix array interval using backwards search
 *  if the interval is empty, it will be initialized for string "b"
 *  if the interval corresponds to string "S", it will be updated for string "bS"
*/
template <size_t AlphabetSize>
inline void update_interval(interval& interval, const fm_index<AlphabetSize>& fm, const uint8_t b) {
		if (interval.empty()) {
				interval.first = fm.C(b);
				interval.last = fm.C(b+1);
		} else {
			//NOTE: computation of occ(.,last) might be faster using occ(.,first)
			interval.first = fm.C(b) + (interval.first>0?fm.occ(b,interval.first-1):0);
	    interval.last = fm.C(b) + fm.occ(b,interval.last-1);
		}
}

};

#endif