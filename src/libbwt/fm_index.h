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
 		typedef typename String::allocator_type::size_type size_type;
 		typedef typename String::allocator_type::value_type value_type;
		typedef std::array<size_type,AlphabetSize> alpha_count;
		std::vector<alpha_count> _occ;
    
    // The C(a) array
    std::array<size_type,AlphabetSize+1> C;
    
    // Reference to the Burrow-Wheeler encoded string
    const String& bwt;

    //! \brief constructor
    fm_index(const String& bwt): bwt(bwt),_occ(bwt.size()) {
    		// compute the # of occurence of a character c in C[c]
    		// and put in _occ[i][c] the # of occurence of c in bwt[0..i]
    		size_type k=0;
    		std::fill(C.begin(),C.end(),0);
				auto i = _occ.begin();					
    		std::fill(i->begin(),i->end(),0);
    		for(auto c:bwt) {
    				++C[c];
    				++(i->operator[](c));
    				auto j = i;
    				++i;
    				*i = *j;
    		}        		
    		// update C[c] to contain the # of occurence of all lexicography smaller characters
				size_type s = 0;
    		for(auto& i:C) {
    				auto v = i;
    				i = s;
    				s += v;
    		}
    		assert(C.back()==bwt.size());
    }

    //! \brief last to first mapping
    //! \return suffix array rank for character bwt[i]
    inline size_type lf_rank(size_type i) const {
    		auto c = bwt[i];
        return C[c] + occ(i)[c];
    }

    //! \return number of occurence of symbol c in bwt[0..i]
    inline size_type occ(const value_type c, size_type i) const {return _occ[i][c];}

		//! return the suffix array interval for character b
		inline interval sa_interval(const value_type b) const {return interval(C[b],C[b+1]);}
	
		//! \brief update a suffix array interval using backwards search
		//! if the given interval corresponds to string S, it will be updated for string bS
		inline void update_sa_interval(interval& interval, const value_type b) const {
				assert(interval.lower>0);
		    interval.lower = C[b] + occ(b,interval.lower-1);
		    interval.upper = C[b] + occ(b,interval.upper-1);
		}
};


#endif
