#ifndef RLESTR_H
#define RLESTR_H


#include <cassert>
#include <cstdint>
#include <vector>
#include <numeric>
#include <functional>

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




struct rle_string {
    // The total length of the bw string
    size_t m_numSymbols;
		std::vector<rle_unit> runs;
		
		// Append a symbol to the bw string
		void read(const std::string& filename);

		// Append a symbol to the bw string
		void append(uint8_t b) {
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
		    ++m_numSymbols;
		}
		
};







#endif
