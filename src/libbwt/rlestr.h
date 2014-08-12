#ifndef RLESTR_H
#define RLESTR_H


#include <cassert>
#include <cstdint>
#include <vector>
#include <numeric>
#include <functional>

/*! \class RLUnit
 *  \brief A run-length encoded unit of the FM-index
 * 
 *  A unit of the RLBWT is a pair of a symbol and its count. The high 3 bits encodes the symbol to store. The low 5 bits encodes the length of the run.
*/
struct RLUnit {
		// constant masks
		static const uint8_t RL_COUNT_MASK = 0x1F;  //00011111
		static const uint8_t RL_SYMBOL_MASK = 0xE0; //11100000
		static const uint8_t RL_FULL_COUNT = 31;
		static const uint8_t RL_SYMBOL_SHIFT = 5;
		
    // variable where the data is stored
    uint8_t data;
		
    RLUnit() : data(0) {}
    RLUnit(uint8_t b) : data(1) {setChar(b);}

    //! \return true if the count cannot be incremented
    inline bool isFull() const {return (data & RL_COUNT_MASK) == RL_FULL_COUNT;}
    inline bool isEmpty() const {return (data & RL_COUNT_MASK) == 0;}

    // 
    inline RLUnit& operator++() {
        assert(!isFull());
        ++data;
        return *this;
    }

    // 
    inline RLUnit& operator--() {
        assert(!isEmpty());
        --data;
        return *this;
    }    

    inline uint8_t length() const {
        assert((data & RL_COUNT_MASK) != 0);
        return data & RL_COUNT_MASK;
    }

    //! \brief Set the symbol
    inline void setChar(uint8_t symbol) {
        data &= RL_COUNT_MASK; // Clear the current symbol
        uint8_t code = symbol;
        code <<= RL_SYMBOL_SHIFT;
        data |= code;
    }

    //! \brief Get the symbol
    inline uint8_t value() const {
        uint8_t code = data & RL_SYMBOL_MASK;
        code >>= RL_SYMBOL_SHIFT;
        return code;
    }
};




struct RLEString:std::vector<RLUnit>{
    // The total length of the bw string
    size_t m_numSymbols;

		// Append a symbol to the bw string
		void append(char b);
};







#endif
