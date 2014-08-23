#ifndef RLESTRING_H
#define RLESTRING_H

#include <cinttypes>

namespace bwt {
	
	
	/*! \class run_t
	 *  \brief a unit of the rle_string 
	 */
	struct run_t {
	  uint8_t _data;
	  run_t():_data(0) {}
	  run_t(uint8_t val):_data((val<<5) | 1) {}
	  inline uint8_t length() const {return _data & 0x1F;}
	  inline uint8_t value() const {return (_data & 0xE0)>>5;}
	  inline bool full() const {return length()>=31;}
	  inline run_t& operator++() {++_data;return *this;}
	};



	/*! \class rle_string
	 *  \brief run length encoded string
	 */
	struct rle_string {
    std::vector<run_t> runs;
        
    //! \brief construct an object reading the bwt from the input range
    template<typename InputIterator>
    rle_string(InputIterator first,InputIterator last);
    
    //! \brief construct an object by reading the rle-encoded bwt string from a binary stream
    rle_string(std::istream& is) {_size = read_runs(is);}
    
    //! \return total number of symbol in the bwt string
    inline uint64_t size() const {return _size;}

    //! \brief empty the string
    void clear() {runs.clear();_size = 0;}
    
    //! \brief append a character at the end of the string
    void push_back(uint8_t v) {
      if (runs.empty()) {
				runs.push_back(v);
      } else if (!runs.back().full() && runs.back().value()==v) {
				++runs.back();
      } else {
				runs.push_back(v);
      }
      _size++;
    }

	  void print_debug_info(std::ostream& os) const {
	    os << "size:" << size() << std::endl;
	    os << "#run:" << runs.size() << std::endl;
	    os << "avg run size:" << (double) size() / runs.size() << std::endl;
	  }

  private:
  	uint64_t _size = 0;
  	
    //! \brief read rle encoded runs from an input stream
    size_t read_runs(std::istream& is);
	};
	
	
  ////////////////////////////////////////////////
  //
  // rle_string class implementation
  //
  ////////////////////////////////////////////////
  
  template<typename InputIterator>
  rle_string::rle_string(InputIterator first,InputIterator last) {
    for(;first != last;first++) push_back(*first);
  }
	
  size_t rle_string::read_runs(std::istream& is) {
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
    runs.clear();
    runs.resize(num_runs);
    is.read(reinterpret_cast<char*>(&runs[0]), num_runs*sizeof(runs[0]));
    return num_symbols;
  }
  
};



#endif
