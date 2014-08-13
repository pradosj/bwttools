#ifndef BWT_H
#define BWT_H


#include "rlestr.h"
#include "Alphabet.h"
#include "FMMarkers.h"
#include <algorithm>
#include <vector>


// return true if x is a power of 2
#define IS_POWER_OF_2(x) ((x) & ((x) - 1)) == 0

// return the x % y given that y is a power of 2
#define MOD_POWER_2(x, y) (x) & ((y) - 1)






/*! \class bwt
 *  \brief Run-length encoded Burrows Wheeler transform
 */
class bwt {
    public:
    
        // Constructors
        bwt(const std::string& filename, int sampleRate = 128);

        //
        void initializeFMIndex();


        inline char operator[](size_t idx) const {
            // Calculate the Marker who's position is not less than idx
            const LargeMarker& upper = getUpperMarker(idx);
            size_t current_position = upper.getActualPosition();
            assert(current_position >= idx);

            size_t symbol_index = upper.unitIndex; 

            // Search backwards (towards 0) until idx is found
            while(current_position > idx) {
                assert(symbol_index != 0);
                symbol_index -= 1;
                current_position -= m_rlString[symbol_index].length();
            }

            // symbol_index is now the index of the run containing the idx symbol
            const RLUnit& unit = m_rlString[symbol_index];
            assert(current_position <= idx && current_position + unit.length() >= idx);
            return unit.value();
        }


        inline AlphaCount64::value_type getPC(char b) const { return m_predCount[b]; }









        // Return the number of times each symbol in the alphabet appears in bwt[0, idx]
        inline AlphaCount64 getFullOcc(size_t idx) const { 
            // The counts in the marker are not inclusive (unlike the Occurrence class)
            // so we increment the index by 1.
            ++idx;

            const LargeMarker& marker = getNearestMarker(idx);
            size_t current_position = marker.getActualPosition();
            bool forwards = current_position < idx;

            AlphaCount64 running_count = marker.counts;
            size_t symbol_index = marker.unitIndex; 

            if(forwards)
                accumulateForwards(running_count, symbol_index, current_position, idx);
            else
                accumulateBackwards(running_count, symbol_index, current_position, idx);
            return running_count;
        }

        // Adds to the count of symbol b in the range [targetPosition, currentPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateBackwards(AlphaCount64& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const {
            // Search backwards (towards 0) until idx is found
            while(currentPosition != targetPosition) {
                size_t diff = currentPosition - targetPosition;
                assert(currentUnitIndex >= 0);
                --currentUnitIndex;
                size_t count = std::min<size_t>(m_rlString[currentUnitIndex].length(),diff);                
    						running_count[m_rlString[currentUnitIndex].value()] -= count;
                currentPosition -= count;
            }
        }

        // Adds to the count of symbol b in the range [currentPosition, targetPosition)
        // Precondition: currentPosition <= targetPosition
        inline void accumulateForwards(AlphaCount64& running_count, size_t currentUnitIndex, size_t currentPosition, const size_t targetPosition) const {
            // Search backwards (towards 0) until idx is found
            while(currentPosition != targetPosition) {
                size_t diff = targetPosition - currentPosition;
                assert(currentUnitIndex < m_rlString.size());
                size_t count = std::min<size_t>(m_rlString[currentUnitIndex].length(),diff);
                running_count[m_rlString[currentUnitIndex].value()] += count;
                currentPosition += count;
                ++currentUnitIndex;
            }
        }
        
        

        // Return the number of times each symbol in the alphabet appears ins bwt[idx0, idx1]
        inline AlphaCount64 getOccDiff(size_t idx0, size_t idx1) const {return getFullOcc(idx1) - getFullOcc(idx0);}
        inline size_t getBWLen() const { return m_rlString.m_numSymbols; }
        inline size_t getNumRuns() const { return m_rlString.size(); }

        // Return the first letter of the suffix starting at idx
        inline char getF(size_t idx) const {
            size_t ci = 0;
            while(ci < BWT_ALPHABET::ALPHABET_SIZE && m_predCount[ci] <= idx)
                ci++;
            assert(ci != 0);
            return BWT_ALPHABET::RANK_ALPHABET[ci - 1];
        }

        // Print the size of the BWT
        void printInfo() const;
        void print() const;
        void printRunLengths() const;

        // IO
        friend class BWTReaderBinary;

    private:
        // Default constructor is not allowed
        bwt() {}
        
        // Calculate the number of markers to place
        size_t getNumRequiredMarkers(size_t n, size_t d) const;

        // Get the index of the marker nearest to position in the bwt
        inline size_t getNearestMarkerIdx(size_t position, size_t sampleRate, size_t shiftValue) const {
            size_t offset = MOD_POWER_2(position, sampleRate); // equivalent to position % sampleRate
            size_t baseIdx = position >> shiftValue;
            return (offset < (sampleRate >> 1))?baseIdx:baseIdx+1;
        }

        // Return a LargeMarker with values that are interpolated by adding
        // the relative count nearest to the requested position to the last
        // LargeMarker
        inline LargeMarker getInterpolatedMarker(size_t target_small_idx) const {
            // Calculate the position of the LargeMarker closest to the target SmallMarker
            size_t target_position = target_small_idx << m_smallShiftValue;
            size_t curr_large_idx = target_position >> m_largeShiftValue;

            LargeMarker absoluteMarker = m_largeMarkers[curr_large_idx];
            const SmallMarker& relative = m_smallMarkers[target_small_idx];
            absoluteMarker.counts += relative.counts;
            absoluteMarker.unitIndex += relative.unitCount;
            return absoluteMarker;
        }

        // Get the interpolated marker with position closest to position
        inline LargeMarker getNearestMarker(size_t position) const {
            size_t nearest_small_idx = getNearestMarkerIdx(position, m_smallSampleRate, m_smallShiftValue);
            return getInterpolatedMarker(nearest_small_idx);
        }

        // Get the greatest interpolated marker whose position is less than or equal to position
        inline LargeMarker getLowerMarker(size_t position) const {
            size_t target_small_idx = position >> m_smallShiftValue;
            return getInterpolatedMarker(target_small_idx);
        }

        // Get the lowest interpolated marker whose position is strictly greater than position
        inline LargeMarker getUpperMarker(size_t position) const {
            size_t target_small_idx = (position >> m_smallShiftValue) + 1;
            return getInterpolatedMarker(target_small_idx);
        }


        // The C(a) array
        AlphaCount64 m_predCount;
        
        // The run-length encoded string
        RLEString m_rlString;

        // The marker vector
        std::vector<LargeMarker> m_largeMarkers;
        std::vector<SmallMarker> m_smallMarkers;

        // The sample rate used for the markers
        size_t m_largeSampleRate;
        size_t m_smallSampleRate;

        // The amount to shift values by to divide by m_sampleRate
        int m_smallShiftValue;
        int m_largeShiftValue;

};
#endif
