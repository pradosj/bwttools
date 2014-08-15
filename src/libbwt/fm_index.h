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
		typedef std::array<uint64_t,AlphabetSize> alpha_count64;
		std::vector<alpha_count64> _occ;
				
    public:
	      // The C(a) array
	      std::array<uint64_t,AlphabetSize+1> C;
	      
	      // Reference to the Burrow-Wheeler encoded string
	      const String& bwt;
    
        // Constructors
        fm_index(const String& bwt): bwt(bwt),_occ(bwt.size()) {
        		// compute the # of occurence of a character c in C[c]
        		// and put in _occ[i][c] the # of occurence of c in bwt[0..i]
        		uint64_t k=0;
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
						uint64_t s = 0;
        		for(auto& i:C) {
        				auto v = i;
        				i = s;
        				s += v;
        		}
        		assert(C.back()==bwt.size());
        }

        // Return the first letter of the suffix starting at idx
        inline uint64_t LF(uint64_t i) const {
        		auto c = bwt[i];
            return C[c] + occ(i)[c];
        }

        // Return the number of times each symbol in the alphabet appears in bwt[0, idx]
        inline alpha_count64 occ(uint64_t i) const {return _occ[i];}

				// Initialize the interval of index idx to be the range containining all the b suffixes
				inline interval init_interval(const uint8_t b) const {return interval(C[b],C[b+1]);}
			
				// Update the given interval using backwards search
				// If the interval corresponds to string S, it will be updated for string bS
				inline void update_interval(interval& interval, uint8_t b) const {
						assert(interval.lower>0);
				    interval.lower = C[b] + occ(interval.lower-1)[b];
				    interval.upper = C[b] + occ(interval.upper-1)[b];
				}	
};


#endif
