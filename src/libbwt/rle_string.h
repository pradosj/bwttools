#ifndef RLESTRING_H
#define RLESTRING_H


#include <cassert>
#include <cstdint>
#include <vector>
#include <numeric>
#include <functional>
#include <iostream>

namespace bwt {


/*! \class rle_unit
 *  \brief A run-length encoded unit of the FM-index
 * 
 *  A unit of the rle_string is a pair of a symbol and its count. The high 3 bits encodes the symbol to store. The low 5 bits encodes the length of the run.
*/
struct rle_unit {
		// constant masks
		static const uint8_t value_shift = 5;
		static const uint8_t value_mask = 0xE0;  //11100000
		static const uint8_t length_mask = 0x1F; //00011111
		static const uint8_t max_length = 31;
		
    // variable where the data is stored
    uint8_t data;
		
    rle_unit() : data(0) {}
    rle_unit(uint8_t b) : data(1) {set_value(b);}

    //! \return true if the count cannot be incremented
    inline bool full() const {return ((data ^ length_mask) & length_mask) == 0;}
    inline bool empty() const {return (data & length_mask) == 0;}

    // 
    inline rle_unit& operator++() {
        assert(!full());
        ++data;
        return *this;
    }

    // 
    inline rle_unit& operator--() {
        assert(!empty());
        --data;
        return *this;
    }    

    inline uint8_t length() const {
        assert((data & length_mask) != 0);
        return data & length_mask;
    }

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




class rle_string {
    size_t m_size = 0;
    
    // marks used for random access
    const size_t marks_sampling = 1024;
    std::vector<size_t> m_marks; // m_marks[i] is the index of the run containing symbol str[i*marks_sampling]
    std::vector<uint8_t> m_marks_offset; // m_marks_offset[i] is where in the run is symbol str[i*marks_sampling]
    
		//! \brief append a symbol to the end of rle_string
		void push_back(uint8_t b) {
		    if (runs.empty()) {
		   		 	runs.push_back(b);
		    } else {
		        rle_unit& lastUnit = runs.back();
		        if (lastUnit.value() == b && !lastUnit.full()) {
		            ++lastUnit;
		        } else {
		        	runs.push_back(b);
		        }
		    }
		    ++m_size;
		}    
    
    inline void update_size(){
    		m_size = std::accumulate(runs.begin(),runs.end(),0,[](uint64_t s,rle_unit u){return s+u.length();});
    }
    
    inline void init_marks() {
    		m_marks.resize(size()>0?(size()-1)/marks_sampling+1:0);
    		m_marks_offset.resize(m_marks.size());
    		
    		// iterate over the runs to initialize the marks
    		uint64_t ub,lb = 0;
    		size_t mark_idx = 0;
    		uint64_t mark_pos = 0;
    		size_t run_idx = 0;
    		for(auto run:runs) {
    				ub = lb + run.length();
    				if (ub > mark_pos) {
    						m_marks[mark_idx] = run_idx;
    						m_marks_offset[mark_idx] = mark_pos - lb;
    						mark_pos += marks_sampling;
    						++mark_idx;
    				}
    				lb = ub;
    				++run_idx;
    		}
    }
  public:
			std::vector<rle_unit> runs;
			
			// create an empty string
			rle_string() {}
			template<typename T> rle_string(std::initializer_list<T> l) : rle_string(l.begin(),l.end()) {}
			template<typename ForwardIterator> rle_string(ForwardIterator first,ForwardIterator last) {
					for(;first!=last;++first) push_back(*first);
					init_marks();					
			}
			
			
			
			
			
			//! \return total length of the string
			inline size_t size() const { return m_size; }
			
			//! \brief read a rle_string from a binary stream
			friend std::istream& operator>>(std::istream& is,rle_string& str);
			
			std::ostream& print(std::ostream& os);
			
      inline uint8_t operator[](uint64_t i) const {
      		size_t mark_idx = i/marks_sampling;
      		auto run_idx = m_marks[mark_idx];

      		// [lb..ub) is the index range of run_idx
      		uint64_t ub, lb = (mark_idx * marks_sampling) - m_marks_offset[mark_idx];
      		while(1) {
      			ub = lb + runs[run_idx].length();
      			if (ub>i) break;
      			lb = ub;
      			++run_idx;
      		}
					return runs[run_idx].value();
      }
      //const_iterator begin() const {return 0;}
      //const_iterator end() const {return 0;}
};


};




#endif
