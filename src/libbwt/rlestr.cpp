
#include "rlestr.h"




void RLEString::append(char b) {
    bool increment = false;
    if(!m_rlString.empty()) {
        RLUnit& lastUnit = m_rlString.back();
        if(lastUnit.value() == b && !lastUnit.isFull()) {
            ++lastUnit;
            increment = true;
        }
    }

    if(!increment) {
        // Must add a new unit to the string
        m_rlString.push_back(RLUnit(b));
    }
    ++m_numSymbols;
}



