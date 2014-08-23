
#include <iostream>
#include <fstream>
#include <algorithm>
#include <fm_index.h>
#include <algo.h>

int test_fm() {
		const std::string alphabet("$abcdr");
		auto encode = [&](char c){return alphabet.find(c);};
		auto decode = [&](char c){return alphabet[c];};

		std::string bwt("ard$rcaaaabb");
		std::transform(bwt.begin(),bwt.end(),bwt.begin(),encode);
		bwt::fm_index<6> fm(bwt.begin(),bwt.end());
		fm.print_debug_info(std::cout);
		
		std::cout << "C[]:" << std::endl;
		for(auto c:fm.C) std::cout << c << ' ';
		std::cout << std::endl;

		std::cout << "occ[c,i]:" << std::endl;
		for(auto c:{0,1,2,3,4,5}) {
				for(auto i=0;i<bwt.size();i++) std::cout << fm.occ(i)[c];
				std::cout << std::endl;
		}
		
		
		
		uint8_t c;
		decltype(fm)::alpha_count64 low;
		decltype(fm)::alpha_count64 high;
		
		bwt::init_char_range(fm,low,high);
		
		c = encode('a');
		assert(low[c]==1 && high[c]==6);
		extend_lhs(fm,low,high,c);
		
		c = encode('r');
		assert(low[c]==10 && high[c]==12);
		
		
		extend_lhs(fm,low,high,c);
		c = encode('b');
		assert(low[c]==6 && high[c]==8);
		
		return 0;
}




int main(int argc, char* argv[]) {
		test_fm();
		return 0;
}
