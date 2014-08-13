//
// FMMarkers - Marker classes used in the FM-index implementation
//
#ifndef FMMARKERS_H
#define FMMARKERS_H

#include "AlphaCount.h"


/*! \class LargeMarker 
 * \brief To allow random access to the BWT symbols and implement the occurrence array we keep a vector of symbol counts every D1 symbols.
 *        These counts are the absolute number of times each symbol has been seen up to that point.
 */
struct LargeMarker {
    LargeMarker() : unitIndex(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getActualPosition() const {return counts.sum();}

    // Returns true if the data in the markers is identical
    bool operator==(const LargeMarker& rhs) {return counts==rhs.counts && unitIndex == rhs.unitIndex;}

    // The number of times each symbol has been seen up to this marker
    AlphaCount64 counts; 

    // The index in the RLVector of the run that starts after
    // this marker. That is, if C = getActualPosition(), then
    // the run containing the B[C] is at unitIndex. This is not necessary
    // a valid index if there is a marker after the last symbol in the BWT
    size_t unitIndex;
};

/*! \class SmallMarker
 * \brief Small markers contain the counts within an individual block of the BWT. In other words the small marker contains the count for the last D2 symbols
 */
struct SmallMarker {
    SmallMarker() : unitCount(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getCountSum() const {return counts.sum();}

    // The number of times each symbol has been seen up to this marker
    AlphaCount16 counts; 

    // The number of RL units in this block
    uint16_t unitCount;
};


#endif
