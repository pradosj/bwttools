
#include <iostream>
#include <fm_index.h>



int main(int argc, char* argv[]) {
		std::string bwt("ard$rcaaaabb");
		fm_index<std::string,256> fm(bwt);
		
		std::cout << "C[]:" << std::endl;
		for(auto c:fm.C) std::cout << c << ' ';
		std::cout << std::endl;

		std::cout << "occ[c,i]:" << std::endl;
		for(auto c:{'$','a','b','c','d','r'}) {
				for(auto i=0;i<bwt.size();i++) std::cout << fm.occ(i)[c];
				std::cout << std::endl;
		}
		
		interval range = fm.init_interval('a');
		std::cout << range.lower << ':' << range.upper << std::endl;
		fm.update_interval(range,'r');
		std::cout << range.lower << ':' << range.upper << std::endl;
		fm.update_interval(range,'b');
		std::cout << range.lower << ':' << range.upper << std::endl;

		return 0;
}


