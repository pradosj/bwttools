
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
		for(auto c:{0,1,2,3,4,5}) std::cout << fm.C(c) << ' ';
		std::cout << std::endl;

		std::cout << "occ[c,i]:" << std::endl;
		for(auto c:{0,1,2,3,4,5}) {
				for(auto i=0;i<bwt.size();i++) std::cout << fm.occ(c,i);
				std::cout << std::endl;
		}
		
		bwt::interval rg1;
		bwt::update_interval(rg1, fm, encode('a'));
		assert(rg1.first==1 && rg1.last==6);
		bwt::update_interval(rg1, fm, encode('r'));
		assert(rg1.first==10 && rg1.last==12);
		bwt::update_interval(rg1, fm, encode('b'));
		assert(rg1.first==6 && rg1.last==8);
		
		return 0;
}




int main(int argc, char* argv[]) {
		test_fm();
		return 0;
}
