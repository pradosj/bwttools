#ifndef ALGO_H
#define ALGO_H

#include "fm_index.h"

namespace bwt {
		
		/*! \brief initialize the intervals [low[c],high[c]) with the 1-character string "c"
		 *         for each character of the alphabet
     */
    template<size_t sz>
    inline void init_char_range(const fm_index<sz>& fm, typename fm_index<sz>::alpha_count64& low, typename fm_index<sz>::alpha_count64& high) {				
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
				init_char_range(fm,low,high);
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
    
    /*! \brief 1-character suffix extension of the interval [first,last) corresponding to string "S"
     *         Extension is done with all characters of the alphabet.
     *         Bounds of the extended intervals are set in [low[c],high[c])
     */
    template<size_t sz>
    inline void extend_rhs(const fm_index<sz>& fm, typename fm_index<sz>::alpha_count64& low, typename fm_index<sz>::alpha_count64& high, uint64_t first, uint64_t last) {
      if (first>=last) {
				init_char_range(fm,low,high);
      } else {
				// compute diff = fm.occ(last-1) - fm.occ(first-1)
				high = fm.occ(last-1);
				if (first>0) std::transform(high.begin(),high.end(),fm.occ(first-1).begin(),high.begin(),std::minus<uint64_t>());
				
				// update low and high using diff
				uint64_t s = 0;
				auto i = low.begin();
				for(auto& v:high) {
					*i += s;
					s += v;
					v += *i;
					++i;
				}
			}
    }
    
};



#endif
