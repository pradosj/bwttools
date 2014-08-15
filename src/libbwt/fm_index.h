#ifndef FMINDEX_H
#define FMINDEX_H


#include "interval.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include <cassert>
#include <string>
#include <array>
	




/*! \class run
 *  \brief A run-length encoded unit of the FM-index
 * 
 *  A unit of the rle_string is a pair (symbol,count). 
 *  The high 3 bits encodes the symbol to store. 
 *  The low 5 bits encodes the length of the run.
*/
struct run {
		// constant masks
		static const uint8_t value_shift = 5;
		static const uint8_t value_mask = 0xE0;  //11100000
		static const uint8_t length_mask = 0x1F; //00011111
		static const uint8_t max_length = 31;
		
    // variable where the data is stored
    uint8_t data;
		
    run() : data(0) {}
    run(uint8_t b) : data(1) {set_value(b);}

    //! \return true if the count cannot be incremented
    inline bool full() const {return ((data ^ length_mask) & length_mask) == 0;}
    inline bool empty() const {return (data & length_mask) == 0;}

    // 
    inline run& operator++() {
        assert(!full());
        ++data;
        return *this;
    }

    // 
    inline run& operator--() {
        assert(!empty());
        --data;
        return *this;
    }    

    inline uint8_t length() const {return data & length_mask;}

    //! \brief Set the symbol
    inline void set_value(uint8_t symbol) {
        data &= length_mask; // Clear the current symbol
        uint8_t code = symbol;
        code <<= value_shift;
        data |= code;
    }

    //! \brief Get the symbol
    inline uint8_t value() const {
        uint8_t code = data & value_mask;
        code >>= value_shift;
        return code;
    }
};






/*! \class fm_index
 *  \brief Run-length encoded Burrows Wheeler transform
 */
template <typename String,size_t AlphabetSize>
class fm_index {	
	public:
		typedef std::array<uint64_t,AlphabetSize> alpha_count64;
 		

    // Reference to the Burrow-Wheeler encoded string
    const String& bwt;

    //! \brief constructor
    fm_index(const String& bwt): bwt(bwt) {
    		// compute the # of occurence of a character c in C[c]
    		// and put in _occ[i][c] the # of occurence of c in bwt[0..i]
				alpha_count64 n;
    		std::fill(_C.begin(),_C.end(),0);
    		_marks64.reserve(bwt.size()>>16);
    		_marks16.reserve(bwt.size()>>9);
    		
				uint64_t k=0;
				std::fill(n.begin(),n.end(),0);
    		for(auto c:bwt) {
    				++_C[c];
    				++n[c];
    				if (k & 0xFFFF) _marks64.push_back(mark64_t());
    				if (k & 0x01FF) {
    						alpha_count16 m;
    						std::transform(n.begin(),n.end(),_marks64.back().counts.begin(),m.begin(),std::minus<uint64_t>());
    						_marks16.push_back(mark16_t());
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
    inline uint64_t occ(const uint8_t c, uint64_t i) const {return occ(i)[c];}
    inline alpha_count64 occ(uint64_t i) const {
    		auto m = previous_mark(i);
    		auto m_pos = std::accumulate(m.counts.begin(),m.counts.end(),0);
    		auto r = _runs.begin() + m.run_index;
    		while(i >= m_pos + r->length()) {
    				m_pos += r->length();
    				//++m.run_index;
    				m.counts[r->value()] += r->length();
    				++r;
    		}
    		m.counts[r->value()] += (i-m_pos);
    		return m.counts;
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
			//
			// types definitions
			//
			typedef std::array<uint16_t,AlphabetSize> alpha_count16;
			template <typename T1,typename T2> 
			struct mark_t {
					T1 run_index;
					T2 counts;
			};
			typedef mark_t<uint64_t,alpha_count64> mark64_t;
			typedef mark_t<uint16_t,alpha_count16> mark16_t;


			//
			// attributes
			//
			std::array<uint64_t,AlphabetSize+1> _C;
			std::vector<run> _runs; // rle encoded bwt string

			// set one largeMark every 65536 indices, and one smallMark every 512
			// _mark64[i] stores occ(.,k) and run_index(bwt[k]), for k=i*65536
			// _mark16[i] stores occ(.,k) and run_index(bwt[k]), for k=i*512 express relatively to the corresponding _mark64
			std::vector<mark64_t> _marks64;
			std::vector<mark16_t> _marks16;

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

};


#endif
