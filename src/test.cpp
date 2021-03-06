
#include <iostream>
#include <fstream>
#include <algorithm>
#include <fm_index.h>
#include <algo.h>



void test_fm() {
	// define alphabet
	const std::string alphabet("$abcdr");
	auto encode = [&](char c){return alphabet.find(c);};
	auto decode = [&](char c){return alphabet[c];};
	
	// load bwt string
	std::string str("ard$rcaaaabb");
	std::transform(str.begin(),str.end(),str.begin(),encode);
	bwt::rle_string bwt(str.begin(),str.end());
	bwt.print_debug_info(std::cout);
	
	// create fm index
	bwt::fm_index<6> fm(bwt);
	fm.print_debug_info(std::cout);
	
	
	{// display C[] and occ() tables
		std::cout << "C[]:" << std::endl;
		for(auto c:fm.C()) std::cout << c << ' ';
		std::cout << std::endl;

		std::cout << "occ[c,i]:" << std::endl;
		for(auto c:{0,1,2,3,4,5}) {
				for(auto i=0;i<bwt.size();i++) std::cout << fm.occ(i)[c];
				std::cout << std::endl;
		}
	}

	{// test extend_lhs
		decltype(fm)::alpha_count64 low;
		decltype(fm)::alpha_count64 high;
		bwt::alpha_range(fm,low,high);
		
		uint8_t c = encode('a');
		assert(low[c]==1 && high[c]==6);
		extend_lhs(fm,low,high,c);
		
		c = encode('r');
		assert(low[c]==10 && high[c]==12);
		extend_lhs(fm,low,high,c);
		
		c = encode('b');
		assert(low[c]==6 && high[c]==8);
	}
}




int main(int argc, char* argv[]) {
	test_fm();
	return 0;
}
