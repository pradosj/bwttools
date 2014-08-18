
#include <iostream>
#include <fstream>
#include <algorithm>
#include <fm_index.h>
#include <interval.h>



int main(int argc, char* argv[]) {
		const std::string alphabet("$abcdr");
		std::string bwt("ard$rcaaaabb");
		std::transform(bwt.begin(),bwt.end(),bwt.begin(),[&](char c){return alphabet.find(c);});
		fm_index<6> fm(bwt.begin(),bwt.end());
		
/*
		std::ifstream f(argv[1],std::ios::binary);
		fm_index<5> fm(f);
*/
		fm.print_info(std::cout);		
		
		std::cout << "C[]:" << std::endl;
		for(auto c:{0,1,2,3,4,5}) std::cout << fm.C(c) << ' ';
		std::cout << std::endl;

		std::cout << "occ[c,i]:" << std::endl;
		for(auto c:{0,1,2,3,4,5}) {
				for(auto i=0;i<bwt.size();i++) std::cout << fm.occ(c,i);
				std::cout << std::endl;
		}
		
		interval range = fm.sa_interval(1);
		std::cout << range.lower << ':' << range.upper << std::endl;
		fm.update_sa_interval(range,5);
		std::cout << range.lower << ':' << range.upper << std::endl;
		fm.update_sa_interval(range,2);
		std::cout << range.lower << ':' << range.upper << std::endl;

		return 0;
}


