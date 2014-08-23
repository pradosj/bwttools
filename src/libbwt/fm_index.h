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


namespace bwt {


  /*! \class fm_index
   *  \brief FM index with an internal run-length encoded Burrows Wheeler Transform string
   *         To speed up random access, the class store one largeMark every 65536 indices, and one smallMark every 512
   */
  template <size_t AlphabetSize>
  class fm_index {
  public:
    //! \brief define an array of numbers for each alphabet character
    typedef std::array<uint64_t,AlphabetSize> alpha_count64;
    
    //! \brief construct an object reading the bwt from the input range
    template<typename InputIterator>
      fm_index(InputIterator first,InputIterator last);
    
    //! \brief construct an object by reading the rle-encoded bwt string from a binary stream
    fm_index(std::istream& is) {
      _bwt_size = read_runs(is);
      init_fm_from_runs();
    }
    
    //! \return size of the alphabet
    size_t alphabet_size() const {return AlphabetSize;}
    
    //! \return total number of symbol in the bwt string
    inline uint64_t size() const {return _bwt_size;}
    
    /*! \brief number of occurence of symbols [0..c) in bwt string. 
     *         This is a public class attribute but should not be modified externaly
     */
    alpha_count64 C; 
    
    //! \return number of occurence of symbol c in bwt[0..i]
    inline alpha_count64 occ(const uint64_t i) const {return mark_at(i).counts;}
    
    //! \return bwt[i], the ith character of bwt string
    inline uint8_t operator[](const uint64_t i) const {return _runs[mark_at(i).run_index].value();}
    
    
    //! \brief output debugging informations to the given stream
    void print_debug_info(std::ostream& os);
    
  private:
    //
    // internal constants
    //
    static const uint8_t shift64 = 16;
    static const uint8_t shift16 = 7;
    
    //
    // internal types definitions
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
    struct run_t {
      uint8_t _data;
      run_t():_data(0) {}
      run_t(uint8_t val):_data((val<<5) | 1) {}
      inline uint8_t length() const {return _data & 0x1F;}
      inline uint8_t value() const {return (_data & 0xE0)>>5;}
      inline bool full() const {return length()>=31;}
      inline run_t& operator++() {++_data;return *this;}
    };
    
    
    //
    // attributes
    //
    uint64_t _bwt_size=0;
    std::vector<run_t> _runs; // rle encoded bwt string
    std::vector<mark64_t> _marks64; // _marks64[i] stores the index I=run_index(bwt[k]) for k=i*65536, and occ(.,I)
    std::vector<mark16_t> _marks16; // _marks16[i] stores the index I=run_index(bwt[k]) for k=i*512, and occ(.,I) expressed relatively to the preceeding _marks64
    
    
        
    //! \return the run index containing symbol i and occ(.,i)
    inline mark64_t mark_at(const uint64_t i) const {
    	// retreive the preeceding mark
      auto m64 = _marks64[i>>shift64];
      const auto& m16 = _marks16[i>>shift16];
      m64.run_index += m16.run_index;
      std::transform(m64.counts.begin(),m64.counts.end(),m16.counts.begin(),m64.counts.begin(),std::plus<uint64_t>());
    	
    	// interpolate the mark to the requested position
      auto run = _runs.begin() + m64.run_index;
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
    
    
    //! \brief iterate over _runs to initialize _C, bwt_size, _marks16 and marks64
    void init_fm_from_runs();
    
    //! \brief read rle encoded runs from an input stream
    size_t read_runs(std::istream& is);	    
  };
  
  




  
  ////////////////////////////////////////////////
  //
  // fm_index class implementation
  //
  ////////////////////////////////////////////////
  
  
  template <size_t AlphabetSize>
  template<typename InputIterator>
  fm_index<AlphabetSize>::fm_index(InputIterator first,InputIterator last) {
    for(;first != last;first++) {
      if (_runs.empty()) {
				_runs.push_back(*first);
      } else if (!_runs.back().full() && _runs.back().value()==*first) {
				++_runs.back();
      } else {
				_runs.push_back(*first);
      }
      _bwt_size++;
    }
    init_fm_from_runs();
  }
  

  template <size_t AlphabetSize>
  void fm_index<AlphabetSize>::init_fm_from_runs() {
    _marks64.reserve((size()>>shift64) + 1);
    _marks16.reserve((size()>>shift16) + 1);
    
    std::fill(C.begin(),C.end(),0);
    uint64_t run_index=0;
    uint64_t run_pos=0;
    for(auto run:_runs) {
      if (run_pos >= _marks64.size()<<shift64) _marks64.push_back(mark64_t(run_index,C));
      if (run_pos >= _marks16.size()<<shift16) {
				_marks16.push_back(mark16_t(run_index - _marks64.back().run_index));
				std::transform(C.begin(),C.end(),_marks64.back().counts.begin(),_marks16.back().counts.begin(),std::minus<uint64_t>());
      }
      C[run.value()] += run.length();
      run_pos += run.length();
      ++run_index;
    }
    
    // C[c] is the count symbol c in bwt string
    // transform it into the count of lexicography smaller symbols [0..c)
    uint64_t s = 0;
    for(auto& i:C) {
      auto v = i;
      i = s;
      s += v;
    }
    _bwt_size = s;
  }
  

  template <size_t AlphabetSize>
    size_t fm_index<AlphabetSize>::read_runs(std::istream& is) {
    enum {BWF_NOFMI = 0,BWF_HASFMI} flag;
    uint16_t magic_number;
    size_t num_strings, num_symbols, num_runs;
    is.read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
    if (magic_number != 0xCACA) throw std::runtime_error("BWT file is not properly formatted: the magic number provided in file header doesn't correspond to the expected one");
    is.read(reinterpret_cast<char*>(&num_strings), sizeof(num_strings));
    is.read(reinterpret_cast<char*>(&num_symbols), sizeof(num_symbols));
    is.read(reinterpret_cast<char*>(&num_runs), sizeof(num_runs));
    is.read(reinterpret_cast<char*>(&flag), sizeof(flag));
    std::cerr << "#symbols:" << num_symbols << std::endl;
    std::cerr << "#strings:" << num_strings << std::endl;
    std::cerr << "#run:" << num_runs << std::endl;
    _runs.clear();
    _runs.resize(num_runs);
    is.read(reinterpret_cast<char*>(&_runs[0]), num_runs*sizeof(_runs[0]));
    return num_symbols;
  }
  

  template <size_t AlphabetSize>
    void fm_index<AlphabetSize>::print_debug_info(std::ostream& os) {
    os << "size:" << size() << std::endl;
    os << "#run:" << _runs.size() << std::endl;
    os << "avg run size:" << (double)size() / _runs.size() << std::endl;
    os << "#marks64:" << _marks64.size() << " (" << (double) _marks64.size() * sizeof(mark64_t)/1024/1024 << "Mo)" << std::endl;
    os << "#marks16:" << _marks16.size() << " (" << (double) _marks16.size() * sizeof(mark16_t)/1024/1024 << "Mo)" << std::endl;
  }
  
};

#endif
