#ifndef FMINDEX_H
#define FMINDEX_H


#include "interval.h"
#include <algorithm>
#include <vector>
#include <cassert>
#include <string>
#include <array>
	
	

/*! \class fm_index
 *  \brief Run-length encoded Burrows Wheeler transform
 */
template <typename String,size_t AlphabetSize>
class fm_index {	
	public:
 		typedef typename String::allocator_type::value_type uint8_t;

    // Reference to the Burrow-Wheeler encoded string
    const String& bwt;

    //! \brief constructor
    fm_index(const String& bwt): bwt(bwt) {
    		// compute the # of occurence of a character c in C[c]
    		// and put in _occ[i][c] the # of occurence of c in bwt[0..i]
    		std::fill(_C.begin(),_C.end(),0);
    		_lmark.reserve(bwt.size()>>16);
    		_smark.reserve(bwt.size()>>9);
    		
				uint64_t k=0;
				alpha_count n;
				std::fill(n.begin(),n.end(),0);
    		for(auto c:bwt) {
    				++_C[c];
    				++n[c];
    				if (k & 0xFFFF) _lmark.push_back(std::make_pair(0,n));
    				if (k & 0x01FF) {
    						alpha_count16 m;
    						std::transform(n.begin(),n.end(),_lmark.back().second.begin(),m.begin(),std::minus<uint64_t>());
    						_smark.push_back(std::make_pair(0,m));
    				}
    				++k;
    		}
    		// update C[c] to contain the # of occurence of all lexicography smaller characters
				uint64_t s = 0;
    		for(auto& i:_C) {
    				auto v = i;
    				i = s;
    				s += v;
    		}
    		assert(_C.back()==bwt.size());
    }

    //! \brief last to first mapping
    //! \return suffix array rank for character bwt[i]
    inline uint64_t lf_rank(uint64_t i) const {
    		auto c = bwt[i];
        return C(c) + occ(c,i);
    }

    //! \return number of occurence of symbol c in bwt[0..i]
    inline uint64_t C(const uint8_t c) const {return _C[c];}

    //! \return number of occurence of symbol c in bwt[0..i]
    inline uint64_t occ(const uint8_t c, uint64_t i) const {
    		return 0;
    }

		//! return the suffix array interval for character b
		inline interval sa_interval(const uint8_t b) const {return interval(C(b),C(b+1));}
	
		//! \brief update a suffix array interval using backwards search
		//! if the given interval corresponds to string S, it will be updated for string bS
		inline void update_sa_interval(interval& interval, const uint8_t b) const {
				assert(interval.lower>0);
		    interval.lower = C(b) + occ(b,interval.lower-1);
		    interval.upper = C(b) + occ(b,interval.upper-1);
		}
		
		
		
		
		
	private:
			typedef std::array<uint64_t,AlphabetSize> alpha_count;
			typedef std::array<uint16_t,AlphabetSize> alpha_count16;
	    std::array<uint64_t,AlphabetSize+1> _C;

		
			// set one largeMark every 65536 indices, and one smallMark every 512
			// the large mark _lmark[k] stores occ(.,k*65536), and run_index(bwt[k*65536])
			// the small mark _smark[k] stores the difference occ(.,k*512) - occ(.,k*512 - k*512 % 65536), and the difference run_index(bwt[k*512]) - run_index(bwt[k*512 % 65536])
			std::vector< std::pair<uint64_t,alpha_count> > _lmark;
			std::vector< std::pair<uint16_t,alpha_count16> > _smark;
			
			//! \return the interpolated mark preceding i
			inline std::pair<uint64_t,alpha_count> previous_mark(uint64_t i) const {
					auto lm = _lmark[i>>16];
					auto sm = _smark[i>>9];
					lm.first += sm.first;
					std::transform(lm.second.begin(),lm.second.end(),sm.second.begin,lm.begin(),std::plus<uint64_t>());
					return lm;
			}

};


#endif
