#ifndef ALGO_H
#define ALGO_H

#include "fm_index.h"

namespace bwt {
		
		/*! \brief initialize the intervals [low[c],high[c]) with the 1-character string "c"
		 *         for each character of the alphabet
     */
    template<size_t sz>
    inline void alpha_range(const fm_index<sz>& fm, typename fm_index<sz>::alpha_count64& low, typename fm_index<sz>::alpha_count64& high) {				
      low = fm.C();
      std::copy(low.begin()+1,low.end(),high.begin());
      high.back() = fm.bwt().size();
    }
	
    /*! \brief 1-character prefix extension of the interval [first,last) corresponding to string "S"
     *         Extension is done with all characters of the alphabet.
     *         Bounds of the extended intervals are set in [low[c],high[c])
     */
    template<size_t sz>
    inline void extend_lhs(const fm_index<sz>& fm, typename fm_index<sz>::alpha_count64& low, typename fm_index<sz>::alpha_count64& high, uint64_t first, uint64_t last) {
      if (first>=last) {
      	std::fill(low.begin(),low.end(),first);
      	std::fill(high.begin(),high.end(),last);
      } else {
				//NOTE: computation of occ(.,last) might be faster using occ(.,first)
				low = high = fm.C();
				if (first>0) std::transform(low.begin(),low.end(),fm.occ(first-1).begin(),low.begin(),std::plus<uint64_t>());
				std::transform(high.begin(),high.end(),fm.occ(last-1).begin(),high.begin(),std::plus<uint64_t>());
      }
    }
    
    /*! \brief 1-character prefix extension of the interval [low[b],high[b]) corresponding to string "S"
     *         Extension is done with all characters of the alphabet.
     *         Bounds of the extended intervals are set in [low[c],high[c])
     */
    template<size_t sz>
    inline void extend_lhs(const fm_index<sz>& fm, typename fm_index<sz>::alpha_count64& low,typename fm_index<sz>::alpha_count64& high, uint8_t b) {
      extend_lhs(fm,low,high,low[b],high[b]);
    }
        
};



#endif
