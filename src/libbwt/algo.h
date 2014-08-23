#include "fm_index.h"

namespace bwt {

		/*! \brief initialize the intervals [low[c],high[c]) with a the 1 character string "c"
     */
    template<size_t sz>
    inline void init_char_range(const fm_index<sz>& fm, typename fm_index<sz>::alpha_count64& low, typename fm_index<sz>::alpha_count64& high) {				
      low = fm.C;
      std::copy(fm.C.begin()+1,fm.C.end(),high.begin());
      high.back() = fm.size();
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
				low = high = fm.C;
				if (first>0) low += fm.occ(first-1);
				high += fm.occ(last-1);
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