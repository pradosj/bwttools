#ifndef ALPHACOUNT_H
#define ALPHACOUNT_H

#include <cstdint>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <array>



/*! \class AlphaCount
 *  \brief A simple class holding the count for each character of an integer alphabet.
 *         The size of the alphabet is given by a the template parameter AlphabetSize
 */
template<typename Storage,uint8_t AlphabetSize=5>
class AlphaCount: public std::array<Storage,AlphabetSize> {
    public:
        inline AlphaCount() {clear();}
        inline void clear() {std::fill(this->begin(),this->end(),0);}

				//! \return sum of all elements
				inline Storage sum() const {return std::accumulate(this->begin(),this->end(),0);}
				//! \return the sum of elements lower than b
				inline Storage sumLT(uint8_t b) const {return std::accumulate(this->begin(),this->begin()+b,0);}        
        

        // Operators
        inline friend bool operator==(const AlphaCount& left, const AlphaCount& right) {return std::mismatch(left.begin(),left.end(),right.begin()).first == left.end();}
        inline friend bool operator!=(const AlphaCount& left, const AlphaCount& right) {return !(left == right);}
        inline friend AlphaCount operator+(const AlphaCount<Storage>& left, const AlphaCount<Storage>& right) {        		
            AlphaCount out;
            std::transform(left.begin(),left.end(),right.begin(),out.begin(),std::plus<Storage>());
            return out;
        }
        // As the counts are unsigned integers, each value in left
        // must be larger or equal to value in right. The calling function
        // must guarentee this.
        friend AlphaCount operator-(const AlphaCount& left, const AlphaCount& right) {
            AlphaCount out;
            std::transform(left.begin(),left.end(),right.begin(),out.begin(),std::minus<Storage>());
            return out;
        }

        template<typename Storage2> inline AlphaCount& operator+=(const AlphaCount<Storage2>& other) {
            std::transform(this->begin(),this->end(),other.begin(),this->begin(),std::plus<Storage>());
            return *this;
        }

        
        // I/O
        friend std::ostream& operator<<(std::ostream& out, const AlphaCount<Storage>& ac) {
            std::copy(ac.begin(), ac.end(), std::ostream_iterator<Storage>(out, " "));
            return out;
        }
        friend std::istream& operator>>(std::istream& in, AlphaCount<Storage>& ac) {
        		for(auto &e:ac) in >> e;
            return in;
        }
};

// Typedef commonly used AlphaCounts
typedef AlphaCount<uint64_t> AlphaCount64;
typedef AlphaCount<uint16_t> AlphaCount16;



#endif
