#ifndef RLUNIT_H
#define RLUNIT_H


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
    RLUnit(char b) : data(1) {setChar(b);}

    //! \return true if the count cannot be incremented
    inline bool isFull() const {return (data & RL_COUNT_MASK) == RL_FULL_COUNT;}
    inline bool isEmpty() const {return (data & RL_COUNT_MASK) == 0;}
    inline bool isInitialized() const {return data > 0;}

    // Add this run to the AlphaCount
    // Only add up to maxCount symbols. Returns the number
    // of symbols added
    inline size_t addAlphaCount(AlphaCount64& ac, size_t max) const {
        size_t count = std::min<decltype(max)>(getCount(),max);
        ac.add(getChar(), count);
        return count;
    }

    // Add this run to the count of base b if it matches
    // Only add up to maxCount symbols. Returns the number
    // of symbols in the current run, up to max
    inline size_t addCount(char b, size_t& base_count, size_t max) const {
        size_t run_len = std::min<decltype(max)>(getCount(),max);
        if (getChar() == b) base_count += run_len;
        return run_len;
    }    

    // Subtract this run from the AlphaCount
    // Only subtract up to maxCount symbols. Returns the number
    // of symbols added
    inline size_t subtractAlphaCount(AlphaCount64& ac, size_t max) const {
        size_t count = std::min<decltype(max)>(getCount(),max);
        ac.subtract(getChar(), count);
        return count;
    }

		// Subtract this run from the count of base b if it matches
    // Only add up to maxCount symbols. Returns the number
    // of symbols in the current run, up to max
    inline size_t subtractCount(char b, size_t& base_count, size_t max) const {
        size_t run_len = std::min<decltype(max)>(getCount(),max);
        if (getChar() == b) base_count -= run_len;
        return run_len;
    }
        
    // 
    inline void incrementCount() {
        assert(!isFull());
        ++data;
    }

    // 
    inline void decrementCount() {
        assert(!isEmpty());
        --data;
    }    

    inline uint8_t getCount() const {
        assert((data & RL_COUNT_MASK) != 0);
        return data & RL_COUNT_MASK;
    }

    //! \brief Set the symbol
    inline void setChar(char symbol) {
        // Clear the current symbol
        data &= RL_COUNT_MASK;        
        uint8_t code = BWT_ALPHABET::getRank(symbol);
        code <<= RL_SYMBOL_SHIFT;
        data |= code;
    }

    //! \brief Get the symbol
    inline char getChar() const {
        uint8_t code = data & RL_SYMBOL_MASK;
        code >>= RL_SYMBOL_SHIFT;
        return BWT_ALPHABET::getChar(code);
    }
};

#endif
