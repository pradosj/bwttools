
#include <iostream>
#include <rle_string.h>



int main(int argc, char* argv[]) {
	bwt::rle_string s = {1,2,3,3,2,3};
	s.print(std::cout) << std::endl;
	return 0;
}


