#ifndef FMINDEX_H
#define FMINDEX_H


#include "interval.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include <cassert>
#include <string>
#include <array>
	




/*! \class fm_index
 *  \brief Run-length encoded Burrows Wheeler transform
 *  set one largeMark every 65536 indices, and one smallMark every 512
 *  _marks64[i] stores occ(.,k) and run_index(bwt[k]), for k=i*65536
 *  _marks16[i] stores occ(.,k) and run_index(bwt[k]), for k=i*512 express relatively to the corresponding _marks64
 */
template <size_t AlphabetSize>
class fm_index {	
	public:
		typedef std::array<uint64_t,AlphabetSize> alpha_count64;

		fm_index() {
				//read_runs();
				init_fm();
		}



		//! \return total number of symbol in the bwt string
		inline uint64_t size() const {return _bwt_size;}
    //! \return number of occurence of symbols [0..c) in bwt
    inline uint64_t C(uint8_t c) const {return _C[c];}

    //! \return number of occurence of symbol c in bwt[0..i]
    inline uint64_t occ(uint8_t c, uint64_t i) const {return occ(i)[c];}
    inline alpha_count64 occ(uint64_t i) const {
    		auto mark = previous_mark(i);
    		auto run = _runs.begin() + mark.run_index;
    		auto run_pos = std::accumulate(mark.counts.begin(),mark.counts.end(),0) - 1;
    		while(i >= run_pos + rle_length(*run)) {
    				run_pos += rle_length(*run);
    				mark.counts[rle_value(*run)] += rle_length(*run);
    				//++mark.run_index;
    				++run;
    		}
    		mark.counts[rle_value(*run)] += (i-run_pos);
    		return mark.counts;
    }

		//! return the suffix array interval for character b
		inline interval sa_interval(const uint8_t b) const {return interval(C(b),C(b>=AlphabetSize?size():b+1));}
	
		//! \brief update a suffix array interval using backwards search
		//! if the given interval corresponds to string S, it will be updated for string bS
		inline void update_sa_interval(interval& interval, const uint8_t b) const {
				assert(interval.lower>0);
		    interval.lower = C(b) + occ(b,interval.lower-1);
		    //NOTE: computation of occ(.,upper) might be faster using using occ(.,lower)
		    interval.upper = C(b) + occ(b,interval.upper-1);
		}
		
		
		
		
		
	private:
			//
			// types definitions
			//
			typedef std::array<uint16_t,AlphabetSize> alpha_count16;
			template <typename T1,typename T2> 
			struct mark_t {
					T1 run_index;
					T2 counts;
					mark_t(T1 i,T2 c):run_index(i),counts(c) {}
					mark_t(T1 i):run_index(i) {}
			};
			typedef mark_t<uint64_t,alpha_count64> mark64_t;
			typedef mark_t<uint16_t,alpha_count16> mark16_t;


			//
			// attributes
			//
			uint64_t _bwt_size=0;
			alpha_count64 _C;
			std::vector<uint8_t> _runs; // rle encoded bwt string
			std::vector<mark64_t> _marks64;
			std::vector<mark16_t> _marks16;
			
			//
			// RLE access methods
			//
			static inline uint8_t rle_length(uint8_t r) {return r & 0x1F;}
			static inline uint8_t rle_value(uint8_t r) {return (r & 0xE0)>>5;}
			

			/*! \return the run index containg the mark preceding i and 
			 *          occ(.,k) where k is the index of the first symbol of the run
			*/
			inline mark64_t previous_mark(uint64_t i) const {
					auto m64 = _marks64[i>>16];
					const auto& m16 = _marks16[i>>9];
					m64.run_index += m16.run_index;
					std::transform(m64.counts.begin(),m64.counts.end(),m16.counts.begin(),m64.counts.begin(),std::plus<uint64_t>());
					return m64;
			}
			
			//! \brief iterate over _runs to initialize _C, bwt_size, _marks16 and marks64
	    void init_fm() {
	    		//_marks64.reserve(bwt_size>>16);
	    		//_marks16.reserve(bwt_size>>9);
	    		std::fill(_C.begin(),_C.end(),0);
	    		uint64_t run_index=0;
	    		uint64_t run_pos=0;
	    		for(auto run:_runs) {
	    				if (run_pos >= _marks64.size()<<16) _marks64.push_back(mark64_t(run_index,_C));
	    				if (run_pos >= _marks16.size()<<9) {
	    						_marks16.push_back(mark16_t(run_index - _marks64.back().run_index));
	    						std::transform(_C.begin(),_C.end(),_marks64.back().counts.begin(),_marks16.back().counts.begin(),std::minus<uint64_t>());
	    				}
	    				_C[rle_value(run)] += rle_length(run);
	    				run_pos += rle_length(run);
	    				++run_index;
	    		}
	
					// C[c] is the count symbol c in bwt string
					// transform it into the count of lexicography smaller symbols [0..c)
					uint64_t s = 0;
	    		for(auto& i:_C) {
	    				auto v = i;
	    				i = s;
	    				s += v;
	    		}
	    		_bwt_size = s;
	    }
};


#endif
