#ifndef FMINDEX_H
#define FMINDEX_H


#include <algorithm>
#include <numeric>
#include <vector>
#include <cassert>
#include <array>
#include <istream>
#include <iostream>
#include <cinttypes>

#include "rle_string.h"

namespace bwt {


  /*! \class fm_index
   *  \brief FM index with an internal run-length encoded Burrows Wheeler Transform string
   *         To speed up random access, the class store one largeMark every 65536 indices, and one smallMark every 512
   */
  template <size_t AlphabetSize>
  class fm_index {
  public:
    //
    // public types definitions
    //  
    //! \brief define an array of numbers for each alphabet character
    typedef std::array<uint64_t,AlphabetSize> alpha_count64;
        
    //
    // constructors
    //
    fm_index(const rle_string& bwt);

    //
    // methods
    //
    //! \return size of the alphabet
    size_t alphabet_size() const {return AlphabetSize;}
		
    //! \return number of occurence of symbol c in bwt[0..i]
    inline alpha_count64 occ(const uint64_t i) const {return mark_at(i).counts;}
    
    //! \brief the rle_string indexed by the object and storing the BWT
    const rle_string& bwt() const {return _bwt;}
    
    //! \brief number of occurence of symbols [0..c) in bwt string.
    const alpha_count64& C() const {return _C;}
    
    //! \return bwt[i], the ith character of bwt string
    inline uint8_t operator[](const uint64_t i) const {return _bwt.runs()[mark_at(i).run_index].value();}
    
    //! \brief output debugging informations to the given stream
    void print_debug_info(std::ostream& os) const;
    
  private:
    //
    // internal constants
    //
    static const uint8_t shift64 = 16;
    static const uint8_t shift16 = 7;
    
    //
    // internal types
    //
    typedef std::array<uint16_t,AlphabetSize> alpha_count16;
    template <typename T1,typename T2> 
    struct mark_t {
			T1 run_index;
			T2 counts;
      mark_t(T1 i,T2 c):run_index(i),counts(c) {}
      mark_t(T1 i):run_index(i) {}
    };
    typedef mark_t<uint64_t,alpha_count64> mark64_t;
    typedef mark_t<uint16_t,alpha_count16> mark16_t;    
    
    //
    // internal attributes
    //
    std::vector<mark64_t> _marks64; // _marks64[i] stores the index I=run_index(bwt[k]) for k=i*65536, and occ(.,I)
    std::vector<mark16_t> _marks16; // _marks16[i] stores the index I=run_index(bwt[k]) for k=i*512, and occ(.,I) expressed relatively to the preceeding _marks64
    const rle_string& _bwt;
    alpha_count64 _C;
    
    //
    // internal methods
    //
    //! \return the run index containing symbol i and occ(.,i)
    inline mark64_t mark_at(const uint64_t i) const {
    	// retreive the preeceding mark
      auto m64 = _marks64[i>>shift64];
      const auto& m16 = _marks16[i>>shift16];
      m64.run_index += m16.run_index;
      std::transform(m64.counts.begin(),m64.counts.end(),m16.counts.begin(),m64.counts.begin(),std::plus<uint64_t>());
    	
    	// interpolate the mark to the requested position
      auto run = _bwt.runs().begin() + m64.run_index;
      auto run_first = std::accumulate(m64.counts.begin(),m64.counts.end(),0);
      while(true) {
				auto run_len = run->length();
				if (i < run_first + run_len) break;
				run_first += run_len;
				m64.counts[run->value()] += run_len;
				++m64.run_index;
				++run;
      }
      m64.counts[run->value()] += (i+1-run_first);
      return m64;
    }
  };
  
  




  
  ////////////////////////////////////////////////
  //
  // fm_index class implementation
  //
  ////////////////////////////////////////////////

  template <size_t AlphabetSize>
  fm_index<AlphabetSize>::fm_index(const rle_string& bwt): _bwt(bwt) {
    _marks64.reserve((bwt.size()>>shift64) + 1);
    _marks16.reserve((bwt.size()>>shift16) + 1);
    
    std::fill(_C.begin(),_C.end(),0);
    uint64_t run_index=0;
    uint64_t run_pos=0;
    for(auto run:bwt.runs()) {
      if (run_pos >= _marks64.size()<<shift64) _marks64.push_back(mark64_t(run_index,C()));
      if (run_pos >= _marks16.size()<<shift16) {
				_marks16.push_back(mark16_t(run_index - _marks64.back().run_index));
				std::transform(_C.begin(),_C.end(),_marks64.back().counts.begin(),_marks16.back().counts.begin(),std::minus<uint64_t>());
      }
      _C[run.value()] += run.length();
      run_pos += run.length();
      ++run_index;
    }
    
    // C[c] is the count symbol c in bwt string
    // transform it into the count of lexicography smaller symbols [0..c)
    uint64_t s = 0;
    for(auto& i:_C) {
      auto v = i;
      i = s;
      s += v;
    }
    assert(bwt.size()==s);
  }
  

  template <size_t AlphabetSize>
  void fm_index<AlphabetSize>::print_debug_info(std::ostream& os) const {
		_bwt.print_debug_info(os);
    os << "#marks64:" << _marks64.size() << " (" << (double) _marks64.size() * sizeof(mark64_t)/1024/1024 << "Mo)" << std::endl;
    os << "#marks16:" << _marks16.size() << " (" << (double) _marks16.size() * sizeof(mark16_t)/1024/1024 << "Mo)" << std::endl;
  }
  
};

#endif
