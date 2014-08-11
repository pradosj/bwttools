#ifndef ALPHABET_H
#define ALPHABET_H

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <string>

typedef uint64_t BaseCount;

namespace BWT_ALPHABET {
		const uint8_t ALPHABET_SIZE = 5;
		const char RANK_ALPHABET[ALPHABET_SIZE] = {'$', 'A', 'C', 'G', 'T'};
    static const uint8_t s_bwtLexoRankLUT[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
        0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    inline static uint8_t getRank(char b) {return s_bwtLexoRankLUT[static_cast<uint8_t>(b)];}
    inline char getChar(size_t idx) {return RANK_ALPHABET[idx];}
};


//
// A simple class holding the count for each base of a DNA string (plus the terminator)
// Note that this uses RANK_ALPHABET
template<typename Storage>
class AlphaCount {
    public:
        //
        inline AlphaCount() {clear();}

        inline void clear() {memset(m_counts, 0, BWT_ALPHABET::ALPHABET_SIZE * sizeof(Storage));}

        //
        inline void set(char b, Storage v) {m_counts[BWT_ALPHABET::getRank(b)] = v;}

        //
        inline void setByIdx(size_t i, Storage v) {m_counts[i] = v;}


        // 
        inline void increment(char b) {
            int br = BWT_ALPHABET::getRank(b);
            assert(m_counts[br] != maxValue);
            m_counts[br]++;
        }

        //
        inline void add(char b, Storage v) {
            int br = BWT_ALPHABET::getRank(b);
            assert(m_counts[br] + v <= maxValue);
            m_counts[br] += v;
        }

        //
        inline void subtract(char b, Storage v) {
            int br = BWT_ALPHABET::getRank(b);
            assert(m_counts[br] > v);
            m_counts[br] -= v;
        }

        // 
        inline Storage get(char b) const {return m_counts[BWT_ALPHABET::getRank(b)];}

        //
        inline Storage getByIdx(const int i) const {return m_counts[i];}

        // Return the base for index i
        static char getBase(size_t i) {return BWT_ALPHABET::RANK_ALPHABET[i];}
        
        // Return the maximum possible count for a symbol
        static size_t getMaxValue() {return maxValue;}

        // Swap the (A,T) and (C,G) entries, which turns the AlphaCount
        // into the AlphaCount for the complemented sequence
        inline void complement() {
            Storage tmp;

            // A,T
            tmp = m_counts[4];
            m_counts[4] = m_counts[1];
            m_counts[1] = tmp;

            // C,G
            tmp = m_counts[3];
            m_counts[3] = m_counts[2];
            m_counts[2] = tmp;
        }

        // Return the sum of the basecounts for characters lexo. lower than b
        inline size_t getLessThan(char b) const {
            size_t out = 0;
            int stop = BWT_ALPHABET::getRank(b);
            for(int i = 0; i < stop; ++i)
                out += m_counts[i];
            return out;
        }

        // Returns the number of non-zero counts
        uint8_t getNumNonZero() const {
            uint8_t count = 0;
            for(int i = 0; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
                if(m_counts[i] > 0)
                    count += 1;
            }
            return count;
        }

        // Returns true if only one of the DNA characters
        // has a non-zero count
        inline bool hasUniqueDNAChar()
        {
            // index 0 is the '$' character, which we skip
            bool nonzero = false;
            for(int i = 1; i < BWT_ALPHABET::ALPHABET_SIZE; ++i)
            {
                if(m_counts[i] > 0)
                {
                    // if nonzero is set, there is some other nonzero character, return false
                    if(nonzero)
                        return false;
                    else
                        nonzero = true;
                }
            }
            return nonzero;
        }

        // Returns true if any DNA char 
        // has a non-zero count
        inline bool hasDNAChar() {
            // index 0 is the '$' character, which we skip
            for(int i = 1; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
                if(m_counts[i] > 0)
                    return true;
            }
            return false;
        }

        //
        inline char getMaxBase() const {
            char base;
            Storage val;
            getMax(base, val);
            return base;
        }

        // 
        inline char getMaxDNABase() const {
            Storage max = 0;
            int maxIdx = 0;
            for(int i = 1; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
                if (m_counts[i] > max) {
                    max = m_counts[i];
                    maxIdx = i;
                }
            }

            assert(max > 0);
            return BWT_ALPHABET::RANK_ALPHABET[maxIdx];
        }
        //
        inline Storage getMaxCount() const {
            char base;
            Storage val;
            getMax(base, val);
            return val;
        }

        //
        inline void getMax(char& base, Storage& val) const {
            Storage max = 0;
            int maxIdx = 0;
            for(int i = 0; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
                if(m_counts[i] > max) {
                    max = m_counts[i];
                    maxIdx = i;
                }
            }
            base = BWT_ALPHABET::RANK_ALPHABET[maxIdx];
            val = max;
        }

        //
        inline size_t getSum() const {
            size_t sum = m_counts[0];
            sum += m_counts[1];
            sum += m_counts[2];
            sum += m_counts[3];
            sum += m_counts[4];
            return sum;
        }

        // Sort the DNA bases into count order and
        // write them to pOut. This does not include the '$' symbol
        inline void getSorted(char* pOut, size_t len) const {
            assert(len >= BWT_ALPHABET::ALPHABET_SIZE);
            (void)len;
            for(size_t i = 0; i < BWT_ALPHABET::ALPHABET_SIZE; ++i)
                pOut[i] = BWT_ALPHABET::RANK_ALPHABET[i];
            std::sort(pOut, pOut+BWT_ALPHABET::ALPHABET_SIZE, AlphaCountCompareDesc(this));
        }

        // Return the unique DNA character described by the alphacount
        // Returns '$' if no such character exists
        // Asserts if more than one character is described by the alphacount
        inline char getUniqueDNAChar() {
            char r = '$';
            for(int i = 1; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
                if(m_counts[i] > 0) {
                    assert(r == '$');
                    // lookup the character in the ranked alphabet, since the elements
                    // are stored in lexographic order
                    r = BWT_ALPHABET::RANK_ALPHABET[i];
                }
            }
            return r;
        }

        // Operators
        friend std::ostream& operator<<(std::ostream& out, const AlphaCount<Storage>& ac) {
            std::copy(ac.m_counts, ac.m_counts+BWT_ALPHABET::ALPHABET_SIZE, std::ostream_iterator<Storage>(out, " "));
            return out;
        }

        friend std::istream& operator>>(std::istream& in, AlphaCount<Storage>& ac) {
            for(size_t i = 0; i < BWT_ALPHABET::ALPHABET_SIZE; ++i)
                in >> ac.m_counts[i];
            return in;
        }
        
        inline friend AlphaCount operator+(const AlphaCount<Storage>& left, const AlphaCount<Storage>& right) {
            AlphaCount out;
            out.m_counts[0] = left.m_counts[0] + right.m_counts[0];
            out.m_counts[1] = left.m_counts[1] + right.m_counts[1];
            out.m_counts[2] = left.m_counts[2] + right.m_counts[2];
            out.m_counts[3] = left.m_counts[3] + right.m_counts[3];
            out.m_counts[4] = left.m_counts[4] + right.m_counts[4];
            return out;
        }

        inline AlphaCount& operator+=(const AlphaCount& other) {
            m_counts[0] += other.m_counts[0];
            m_counts[1] += other.m_counts[1];
            m_counts[2] += other.m_counts[2];
            m_counts[3] += other.m_counts[3];
            m_counts[4] += other.m_counts[4];
            return *this;
        }

        // Specialization for add/subtracing a small AlphaCount to/from a large alphacount
        inline friend void alphacount_add(AlphaCount<uint64_t>& lhs, const AlphaCount<uint8_t>& rhs);
        inline friend void alphacount_subtract(AlphaCount<uint64_t>& lhs, const AlphaCount<uint8_t>& rhs);
        inline friend void alphacount_add16(AlphaCount<uint64_t>& lhs, const AlphaCount<uint16_t>& rhs);
        inline friend void alphacount_subtract16(AlphaCount<uint64_t>& lhs, const AlphaCount<uint16_t>& rhs);

        // As the counts are unsigned integers, each value in left
        // must be larger or equal to value in right. The calling function
        // must guarentee this.
        friend AlphaCount operator-(const AlphaCount& left, const AlphaCount& right)
        {
            AlphaCount out;
            out.m_counts[0] = left.m_counts[0] - right.m_counts[0];
            out.m_counts[1] = left.m_counts[1] - right.m_counts[1];
            out.m_counts[2] = left.m_counts[2] - right.m_counts[2];
            out.m_counts[3] = left.m_counts[3] - right.m_counts[3];
            out.m_counts[4] = left.m_counts[4] - right.m_counts[4];
            return out;
        }

        inline friend bool operator==(const AlphaCount& left, const AlphaCount& right)
        {
            return left.m_counts[0] == right.m_counts[0] && left.m_counts[1] == right.m_counts[1] &&
                   left.m_counts[2] == right.m_counts[2] && left.m_counts[3] == right.m_counts[3] &&
                   left.m_counts[4] == right.m_counts[4];   
        }
        
        inline friend bool operator!=(const AlphaCount& left, const AlphaCount& right)
        {
            return !(left == right);
        }
                

    private:
        Storage m_counts[BWT_ALPHABET::ALPHABET_SIZE];
        const static size_t maxValue;

        struct AlphaCountCompareDesc
        {
            AlphaCountCompareDesc(const AlphaCount* p) : pAC(p) {}
            bool operator()(char a, char b)
            {
                return pAC->get(a) > pAC->get(b);
            }
            const AlphaCount* pAC;
        };

};

// Typedef commonly used AlphaCounts
typedef AlphaCount<uint64_t> AlphaCount64;
typedef AlphaCount<uint16_t> AlphaCount16;

//
inline void alphacount_add16(AlphaCount64& lhs, const AlphaCount16& rhs) {
    lhs.m_counts[0] += rhs.m_counts[0];
    lhs.m_counts[1] += rhs.m_counts[1];
    lhs.m_counts[2] += rhs.m_counts[2];
    lhs.m_counts[3] += rhs.m_counts[3];
    lhs.m_counts[4] += rhs.m_counts[4];
}

inline void alphacount_subtract16(AlphaCount64& lhs, const AlphaCount16& rhs) {
    // This function should only be used when lhs is larger the rhs
    // the calling function must guarentee this
    lhs.m_counts[0] -= rhs.m_counts[0];
    lhs.m_counts[1] -= rhs.m_counts[1];
    lhs.m_counts[2] -= rhs.m_counts[2];
    lhs.m_counts[3] -= rhs.m_counts[3];
    lhs.m_counts[4] -= rhs.m_counts[4];
}


inline std::string reverse(const std::string& str) {
	std::string rstr(str);
  std::reverse(rstr.begin(),rstr.end());
  return(rstr);
}

inline char complement_bp(char bp) {
    switch(bp) {
        case 'A':return 'T';
        case 'C':return 'G';
        case 'G':return 'C';
        case 'T':return 'A';
        case 'N':return 'N';
        default:
            assert(false && "Unknown base!");
            return 'N';
    }
}

inline std::string complement(const std::string& str) {
    std::string s(str);
    std::transform(s.begin(), s.end(), s.begin(), complement_bp);
    return(s);
}



#endif
