//
// Occurrence.h - Data structure holding the number of times
// the letter b appears in the string S before position i
//
#ifndef OCCURRENCE_H
#define OCCURRENCE_H

#include "AlphaCount.h"
#include "Alphabet.h"
#include <vector>

// Power of 2 macros

// return true if x is a power of 2
#define IS_POWER_OF_2(x) ((x) & ((x) - 1)) == 0

// return the x % y given that y is a power of 2
#define MOD_POWER_2(x, y) (x) & ((y) - 1)

class Occurrence {
    public:
        
        // Constructors
        Occurrence() : m_sampleRate(1) {}

        void set(char a, size_t i, BaseCount s);
        void print() const;
        size_t getByteSize() const;
        size_t size() const { return m_values.size(); }
        size_t getSampleRate() const { return m_sampleRate; }

        // Calculate the amount a value should be shifted to perform a division
        // by divisor. The divisor must be a power of 2
        static int calculateShiftValue(int divisor);

        friend std::ostream& operator<<(std::ostream& out, const Occurrence& o);
        friend std::istream& operator>>(std::istream& in, Occurrence& o);

    private:

        int m_shift;
        int m_sampleRate;
        std::vector<AlphaCount64> m_values;
};


#endif
