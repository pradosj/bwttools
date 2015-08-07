#include <Rcpp.h>
#include <array>
#include <algorithm>
#include <cassert>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

typedef struct {
	uint64_t u,v;
	//v: line number in the original text with the character pointed by the structure
	//u: number of occurence of the character in bwt before its current position
} pair64_t;



void bcr(const uint8_t* text_begin, const uint8_t* text_end, uint8_t* bwt_begin) {
	if (text_begin == text_end) return;
	
	// find EOL, and initialize an array of iterators on the last character
	std::vector<const uint8_t*> eol;
	for(auto i=text_begin;i!=text_end;++i) {
		i = std::find(i,text_end,0);
		if (i==text_end) break;
		eol.push_back(i-1);
	}
	
	// initialize
	std::vector<pair64_t> a(eol.size()),aa(eol.size());
	uint64_t k=0;for(auto &x:a) x.u = x.v = k++;
	uint8_t *bwt_end = bwt_begin + std::distance(text_begin,text_end);
	uint8_t *bwti_begin(bwt_end),*bwtj_begin(bwt_end);
	std::array<uint64_t,256> mc,mc2,b,ac;

	// core loop
	while (!a.empty()) {
		// initialize loop variables
		uint64_t pre = 0;
		bwtj_begin -= a.size();
		const uint8_t *p = bwti_begin;
		uint8_t *q = bwtj_begin;
		mc.fill(0);
		mc2.fill(0);
		
		// iterate over last characters of the sorted lines
		auto n = a.begin();
		for (auto &u:a) {
			for (uint64_t l = 0; l != u.u - pre; ++l) {
				++mc[*p];
				*q++ = *p++;
			}
			
			auto c = (eol[u.v]>=text_begin?*eol[u.v]:0);
			*q++ = c;
			pre = u.u + 1;
			u.u = mc[c]++;
			if (c) {
				*n++ = u;
				++mc2[c];
			}
		}
		a.resize(std::distance(a.begin(),n));
		
		assert(p==q);
		while(p < bwt_end) ++mc[*(p++)];
		
		// compute pos for next round
		ac[0] = 0;for(int c = 1; c != ac.size(); ++c) ac[c] = ac[c-1] + mc[c-1];
		for(auto &x:a) {
			auto c = (eol[x.v]>=text_begin?*eol[x.v]:0);
			x.u += ac[c] + a.size();
		}

		// stable counting sort
		aa.resize(a.size());
		b[0] = 0;for (int c = 1; c != b.size(); ++c) b[c] = b[c-1] + mc2[c-1];
		for(auto x:a) {
			auto c = (eol[x.v]>=text_begin?*eol[x.v]:0);
			aa[b[c]++] = x;
		}
		std::swap(a,aa);

		// move EOL iterators to previous character
		for(auto &l:eol) if (l>=text_begin && *l) --l;

		bwti_begin = bwtj_begin;
	}
}

// [[Rcpp::export]]
CharacterVector testBcr() {
	const int len = 26;
	const uint8_t* str = (const uint8_t*) "BANA\0BANANA\0ANANAS\0RANANA\0";
	uint8_t bwt[len];
	bcr(str,str+len,bwt);
	std::string BWT((const char*) bwt,len);
	std::replace(BWT.begin(),BWT.end(),'\0','$');
	return Rcpp::CharacterVector::create(BWT);
}



/*** R
testBcr()
*/
