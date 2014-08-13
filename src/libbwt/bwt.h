#ifndef BWT_H
#define BWT_H


#include "rlestr.h"
#include "Alphabet.h"
#include "FMMarkers.h"
#include "BWTInterval.h"
#include <algorithm>
#include <vector>


// return the x % y given that y is a power of 2
#define MOD_POWER_2(x, y) (x) & ((y) - 1)






/*! \class bwt
 *  \brief Run-length encoded Burrows Wheeler transform
 */
class bwt {
    public:
    
        // Constructors
        bwt(const std::string& filename, int smallShift = 7);

				inline size_t getBWLen() const { return m_rlString.m_numSymbols; }
				inline size_t getNumRuns() const { return m_rlString.size(); }

        inline uint8_t symbol(size_t idx) const {
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
        
        

				// Initialize the interval of index idx to be the range containining all the b suffixes
				inline BWTInterval initInterval(const uint8_t b) const {
						return BWTInterval(getPC(b),getPC(b) + getFullOcc(getBWLen() - 1)[b] - 1);
				}
				
				// Update the given interval using backwards search
				// If the interval corresponds to string S, it will be updated for string bS
				inline void updateInterval(BWTInterval& interval, uint8_t b) const {
				    size_t pb = getPC(b);
				    interval.lower = pb + getFullOcc(interval.lower - 1)[b];
				    interval.upper = pb + getFullOcc(interval.upper)[b] - 1;
				}

        inline AlphaCount64::value_type getPC(uint8_t b) const { return m_predCount[b]; }






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

        

        // Return the first letter of the suffix starting at idx
        inline uint8_t getF(uint64_t idx) const {
        		auto i = std::lower_bound(m_predCount.begin(),m_predCount.end(),idx,std::less_equal<uint64_t>());
            return std::distance(m_predCount.begin(),i);
        }

        // IO
        friend class BWTReaderBinary;

    private:
        // Default constructor is not allowed
        bwt() {}
        
        //
        void initializeFMIndex();

				// get the number of markers required to cover the n symbols at sample rate of d
		    // we place a marker at the beginning (with no accumulated counts), every m_sampleRate
		    // bases and one at the very end (with the total counts)
				static inline size_t getNumRequiredMarkers(size_t n, size_t d) {return (n % d == 0) ? (n / d) + 1 : (n / d) + 2;}
        

        // Get the index of the marker nearest to position in the bwt
        inline size_t getNearestMarkerIdx(size_t position) const {
        		size_t baseIdx = position >> m_smallShift;
        		size_t sampleRate = 1<<m_smallShift;
            size_t offset = position & (sampleRate-1);
            return (offset < (sampleRate >> 1))?baseIdx:baseIdx+1;
        }

        // Return a LargeMarker with values that are interpolated by adding
        // the relative count nearest to the requested position to the last
        // LargeMarker
        inline LargeMarker getInterpolatedMarker(size_t target_small_idx) const {
            // Calculate the position of the LargeMarker closest to the target SmallMarker
            size_t target_position = target_small_idx << m_smallShift;
            size_t curr_large_idx = target_position >> m_largeShift;

            LargeMarker absoluteMarker = m_largeMarkers[curr_large_idx];
            const SmallMarker& relative = m_smallMarkers[target_small_idx];
            absoluteMarker.counts += relative.counts;
            absoluteMarker.unitIndex += relative.unitCount;
            return absoluteMarker;
        }

        // Get the interpolated marker with position closest to position
        inline LargeMarker getNearestMarker(size_t position) const {return getInterpolatedMarker(getNearestMarkerIdx(position));}

        // Get the greatest interpolated marker whose position is less than or equal to position
        inline LargeMarker getLowerMarker(size_t position) const {return getInterpolatedMarker(position >> m_smallShift);}

        // Get the lowest interpolated marker whose position is strictly greater than position
        inline LargeMarker getUpperMarker(size_t position) const {return getInterpolatedMarker((position >> m_smallShift) + 1);}




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
 


        // The C(a) array
        AlphaCount64 m_predCount;
        
        // The run-length encoded string
        RLEString m_rlString;

        // The marker vector
        std::vector<LargeMarker> m_largeMarkers;
        std::vector<SmallMarker> m_smallMarkers;

        // The amount to shift values by to divide by m_sampleRate
        int m_smallShift;
        int m_largeShift;

};
#endif
