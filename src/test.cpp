
#include <iostream>
#include <rle_string.h>



int main(int argc, char* argv[]) {
	bwt::rle_string s = {1,2,3,3,3,2,3,0,0,0,0,0,0,1,1,2,2};
	std::cout << s.size() << std::endl;
	s.print(std::cout) << std::endl;
	
	for(auto i=0;i<s.size();i++) {std::cout << (int) s[i] << ' ';}
	std::cout << std::endl;
	
	return 0;
}


