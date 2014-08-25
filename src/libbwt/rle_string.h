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
		//! \brief construct an empty string
		rle_string() {}
		
    //! \brief construct an object reading the string from the input range
    template<typename InputIterator>
    rle_string(InputIterator first,InputIterator last);
        
    //! \return total number of symbol in the bwt string
    inline uint64_t size() const {return _size;}

		//! return the collection of rle runs
		inline const std::vector<run_t>& runs() const {return _runs;}

    //! \brief empty the string
    void clear() {_runs.clear();_size = 0;}
    
    //! \brief append a character at the end of the string
    void push_back(uint8_t v) {
      if (_runs.empty()) {
				_runs.push_back(v);
      } else if (!_runs.back().full() && _runs.back().value()==v) {
				++_runs.back();
      } else {
				_runs.push_back(v);
      }
      _size++;
    }

	  void print_debug_info(std::ostream& os) const {
	    os << "size:" << size() << std::endl;
	    os << "#run:" << _runs.size() << std::endl;
	    os << "#full run:" << std::count_if(_runs.begin(),_runs.end(),[](const run_t& r){return r.full();}) << std::endl;
	    os << "avg run size:" << (double) size() / _runs.size() << std::endl;
	  }

  private:
  	uint64_t _size = 0;
  	std::vector<run_t> _runs;
  	
    //! \brief read rle encoded runs from an input stream
    friend rle_string read_rle_bwt(const std::string& filename);
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
	
  rle_string read_rle_bwt(const std::string& filename) {
  	std::ifstream is(filename,std::ios::binary);
  	
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
    
    rle_string bwt;
    bwt._runs.resize(num_runs);
    is.read(reinterpret_cast<char*>(&bwt._runs[0]), num_runs*sizeof(bwt._runs[0]));
    bwt._size = num_symbols;
    
    return bwt;
  }
  
};



#endif
