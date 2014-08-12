//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Occurrence.cpp - Data structure holding the number of times
// the letter b appears in the string S from S[0..i] (inclusive)
//
#include "Occurrence.h"

// 
int Occurrence::calculateShiftValue(int divisor) {
    assert(divisor > 0);
    assert(IS_POWER_OF_2(divisor));

    // m_sampleRate is a power of 2, count what bit is set
    unsigned int v = divisor;
    unsigned int c = 0; // c accumulates the total bits set in v

    while(v != 1) {
        v >>= 1;
        ++c;
    }
    assert(1 << c == divisor);
    return c;
}

//
void Occurrence::set(char a, size_t i, BaseCount s) {
    m_values[i][a] = s;
}

//
size_t Occurrence::getByteSize() const {
    return m_values.size() * sizeof(AlphaCount64);
}

std::ostream& operator<<(std::ostream& out, const Occurrence& o) {
    out << o.m_sampleRate << "\n";
    out << o.m_values.size() << "\n";
    for(size_t i = 0; i < o.m_values.size(); ++i)
        out << o.m_values[i] << "\n";
    return out;
}

std::istream& operator>>(std::istream& in, Occurrence& o) {
    in >> o.m_sampleRate;
    size_t n;
    in >> n;
    o.m_values.resize(n);
    for(size_t i = 0; i < n; ++i)
        in >> o.m_values[i];
    o.m_shift = Occurrence::calculateShiftValue(o.m_sampleRate);
    return in;
}


//
void Occurrence::print() const {
    for(size_t i = 0; i < m_values.size(); i++) {
        std::cout << m_values[i];
    }
}
