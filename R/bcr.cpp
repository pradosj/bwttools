#include <Rcpp.h>
#include <array>
#include <algorithm>
#include <cassert>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

typedef struct {
	uint64_t u,v;
	uint8_t w;
} pair64_t;



void bcr(const uint8_t* text_begin, const uint8_t* text_end, uint8_t* bwt_begin) {
	
	// split input text into lines
	if (text_begin == text_end) return;
	std::vector<const uint8_t*> lines;
	for(auto i=text_begin;i!=text_end;++i) {
		auto line_end = std::find(i,text_end,0);
		// if the text is not ending with eol, the last line is skipped
		if (line_end==text_end) break;
		lines.push_back(i);
		i = line_end;
	}
	lines.push_back(text_end);

	
	// initialize
	std::vector<pair64_t> a(lines.size()-1),aa;
	uint64_t k=0;for(auto &x:a) x.u = x.v = k++;
	uint8_t *bwt_end = bwt_begin + std::distance(text_begin,text_end);
	uint8_t *bwt0_begin = bwt_begin = bwt_end;

		
	// core loop
	for (long i = 0; !a.empty(); ++i) {
		// initialize loop variables
		uint64_t pre = 0;
		bwt_begin -= a.size();
		const uint8_t *p = bwt0_begin;
		uint8_t *q = bwt_begin;
		std::array<uint64_t,256> mc;mc.fill(0);
		
		// iterate over last characters of the lines ordered according to a[].v
		aa.clear();
		for (auto &u:a) {
			u.w = lines[u.v+1]-2-i >= lines[u.v]? *(lines[u.v+1]-2-i) : 0;
			for (uint64_t l = 0; l != u.u - pre; ++l) {
				++mc[*p];
				*q++ = *p++;
			}
			*q++ = u.w;
			pre = u.u + 1;
			u.u = mc[u.w]++;
			if (u.w) aa.push_back(u);
		}
		a.resize(aa.size());
		
		assert(p==q);//std::copy(p,end,q);
		while(p < bwt_end) ++mc[*(p++)];
		
		std::array<uint64_t,256> ac;
		ac[0] = 0;for(int c = 1; c != ac.size(); ++c) ac[c] = ac[c-1] + mc[c-1];
		for(auto &x:aa) x.u += ac[x.w] + aa.size(); // compute pos for next round

		// stable counting sort ($a[k].w);
		std::array<uint64_t,256> mc2,b;
		mc2.fill(0);for(auto x:aa) ++mc2[x.w];
		b[0] = 0;for (int c = 1; c != b.size(); ++c) b[c] = b[c-1] + mc2[c-1];
		for(auto x:aa) a[b[x.w]++] = x; // this works because $a is already partially sorted
		
		bwt0_begin = bwt_begin;
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
