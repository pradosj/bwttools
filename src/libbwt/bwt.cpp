

#include "bwt.h"
#include <istream>
#include <fstream>


namespace bwt {
	
		// Parse a BWT from a file
		fm_index::fm_index(const std::string& filename, int smallShift): m_largeShift(13), m_smallShift(smallShift) {
				std::fstream(filename,std::ios::binary) >> bwt;
		    init();
		}
		

		
		// Fill in the FM-index data structures
		void fm_index::init() {
				const size_t m_smallSampleRate = 1<<m_smallShift;
				const size_t m_largeSampleRate = 1<<m_largeShift;
			
		    // initialize the marker vectors,
		    // LargeMarkers are placed every 2048 bases (by default) containing the absolute count
		    // of symbols seen up to that point. SmallMarkers are placed every 128 bases with the
		    // count over the last 128 symbols. From these relative counts the absolute count
		    // every 128 symbols can be interpolated.
		    size_t num_large_markers = getNumRequiredMarkers(bwt.size(), m_largeSampleRate);
		    size_t num_small_markers = getNumRequiredMarkers(bwt.size(), m_smallSampleRate);
		    m_largeMarkers.resize(num_large_markers);
		    m_smallMarkers.resize(num_small_markers);
		
		    // Fill in the marker values
		    // We wish to place markers every sampleRate symbols however since a run may
		    // not end exactly on sampleRate boundaries, we place the markers AFTER
		    // the run crossing the boundary ends
		
		    // Place a blank markers at the start of the data
		    m_largeMarkers[0].unitIndex = 0;
		    m_smallMarkers[0].unitCount = 0;
		
		    // State variables for the number of markers placed,
		    // the next marker to place, etc
		    size_t curr_large_marker_index = 1;
		    size_t curr_small_marker_index = 1;
		
		    size_t next_small_marker = m_smallSampleRate;
		    size_t next_large_marker = m_largeSampleRate;
		
		    size_t prev_small_marker_unit_index = 0;
		    size_t running_total = 0;
		    alpha_count64 running_ac;
		
		    for(size_t i = 0; i < bwt.runs.size(); ++i) {
		        // Update the count and advance the running total
		        rle_unit& unit = bwt.runs[i];
		
		        char symbol = unit.value();
		        uint8_t run_len = unit.length();
		        running_ac[symbol] += run_len;
		        running_total += run_len;
		
		        size_t curr_unit_index = i + 1;
		        bool last_symbol = i == bwt.runs.size() - 1;
		
		        // Check whether to place a new large marker
		        bool place_last_large_marker = last_symbol && curr_large_marker_index < num_large_markers;
		        while(running_total >= next_large_marker || place_last_large_marker) {
		            size_t expected_marker_pos = curr_large_marker_index * m_largeSampleRate;
		
		            // Sanity checks
		            // The marker position should always be less than the running total unless 
		            // the number of symbols is smaller than the sample rate
		            assert(expected_marker_pos <= running_total || place_last_large_marker);
		            assert((running_total - expected_marker_pos) <= rle_unit::max_length || place_last_large_marker);
		            assert(curr_large_marker_index < num_large_markers);
		            assert(running_ac.sum() == running_total);
		
		            LargeMarker& marker = m_largeMarkers[curr_large_marker_index];
		            marker.unitIndex = i + 1;
		            marker.counts = running_ac;
		
		            next_large_marker += m_largeSampleRate;
		            curr_large_marker_index += 1;
		            place_last_large_marker = last_symbol && curr_large_marker_index < num_large_markers;
		        }    
		
		        // Check whether to place a new small marker
		        bool place_last_small_marker = last_symbol && curr_small_marker_index < num_small_markers;
		        while(running_total >= next_small_marker || place_last_small_marker) {
		            // Place markers
		            size_t expected_marker_pos = curr_small_marker_index * m_smallSampleRate;
		
		            // Sanity checks
		            // The marker position should always be less than the running total unless 
		            // the number of symbols is smaller than the sample rate
		            assert(expected_marker_pos <= running_total || place_last_small_marker);
		            assert((running_total - expected_marker_pos) <= rle_unit::max_length || place_last_small_marker);
		            assert(curr_small_marker_index < num_small_markers);
		            assert(running_ac.sum() == running_total);
		    
		            // Calculate the number of rl units that are contained in this block
		            if(curr_unit_index - prev_small_marker_unit_index > std::numeric_limits<uint16_t>::max()) {
		                std::cerr << "Error: Number of units in occurrence array block " << curr_small_marker_index 
		                          << " exceeds the maximum value.\n";
		                exit(EXIT_FAILURE);
		
		            }
		
		            // Calculate the large marker to set the relative count from
		            // This is generally the most previously placed large block except it might 
		            // be the second-previous in the case that we placed the last large marker.
		            size_t large_marker_index = expected_marker_pos >> m_largeShift;
		            assert(large_marker_index < curr_large_marker_index); // ensure the last has ben placed
		            LargeMarker& prev_large_marker = m_largeMarkers[large_marker_index];
		
		            // Set the 8bit AlphaCounts as the sum since the last large (superblock) marker
		            alpha_count16 smallAC;
		            for(size_t j = 0; j < smallAC.size(); ++j) {
		                size_t v = running_ac[j] - prev_large_marker.counts[j];
		                if(v > std::numeric_limits<alpha_count16::value_type>::max()) {
		                    std::cerr << "Error: Number of symbols in occurrence array block " << curr_small_marker_index 
		                              << " exceeds the maximum value (" << v << " > " << std::numeric_limits<alpha_count16::value_type>::max() << ")\n";
		                    exit(EXIT_FAILURE);
		                }
		                smallAC[j] = v;
		            }
		            
		            // Set the small marker
		            SmallMarker& small_marker = m_smallMarkers[curr_small_marker_index];
		            small_marker.unitCount = curr_unit_index - prev_large_marker.unitIndex;
		            small_marker.counts = smallAC;
		
		            // Update state variables
		            next_small_marker += m_smallSampleRate;
		            curr_small_marker_index += 1;
		            prev_small_marker_unit_index = curr_unit_index;
		            place_last_small_marker = last_symbol && curr_small_marker_index < num_small_markers;
		        }    
		    }
		
		    assert(curr_small_marker_index == num_small_markers);
		    assert(curr_large_marker_index == num_large_markers);
		
		    // Initialize C(a)
		    m_predCount[0] = 0;
		    for(uint8_t i=1;i<m_predCount.size();i++) {
		    		m_predCount[i] = m_predCount[i-1] + running_ac[i-1];
		    }
		    assert(running_ac.sum()==bwt.size());
		}

};


