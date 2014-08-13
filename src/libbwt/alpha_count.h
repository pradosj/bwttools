#ifndef ALPHACOUNT_H
#define ALPHACOUNT_H

#include <cstdint>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <array>


namespace bwt {
	
/*! \class alpha_count
 *  \brief A simple class holding the count for each character of an integer alphabet.
 *         The size of the alphabet is given by a the template parameter AlphabetSize
 */
template<typename Storage,uint8_t AlphabetSize=5>
struct alpha_count: public std::array<Storage,AlphabetSize> {
        inline alpha_count() {clear();}
        inline void clear() {std::fill(this->begin(),this->end(),0);}

				//! \return sum of all elements
				inline Storage sum() const {return std::accumulate(this->begin(),this->end(),0);}        

        // Operators
        inline friend bool operator==(const alpha_count& left, const alpha_count& right) {return std::mismatch(left.begin(),left.end(),right.begin()).first == left.end();}
        inline friend bool operator!=(const alpha_count& left, const alpha_count& right) {return !(left == right);}
        inline friend alpha_count operator+(const alpha_count<Storage>& left, const alpha_count<Storage>& right) {        		
            alpha_count out;
            std::transform(left.begin(),left.end(),right.begin(),out.begin(),std::plus<Storage>());
            return out;
        }
        // As the counts are unsigned integers, each value in left
        // must be larger or equal to value in right. The calling function
        // must guarentee this.
        friend alpha_count operator-(const alpha_count& left, const alpha_count& right) {
            alpha_count out;
            std::transform(left.begin(),left.end(),right.begin(),out.begin(),std::minus<Storage>());
            return out;
        }

        template<typename Storage2> inline alpha_count& operator+=(const alpha_count<Storage2>& other) {
            std::transform(this->begin(),this->end(),other.begin(),this->begin(),std::plus<Storage>());
            return *this;
        }
};

// Typedef commonly used AlphaCounts
typedef alpha_count<uint64_t> alpha_count64;
typedef alpha_count<uint16_t> alpha_count16;

};

#endif
