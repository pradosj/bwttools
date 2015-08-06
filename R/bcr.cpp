#include <Rcpp.h>
#include <array>
#include <algorithm>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

typedef struct {
	uint64_t u,v;
	uint8_t w;
} pair64_t;



void bcr(const uint8_t* text_begin, const uint8_t* text_end, uint8_t* bwt_begin) {
	auto Tlen = std::distance(text_begin,text_end);
	
	// split $T into short strings at sentinels
	if (text_begin == text_end) return;
	std::vector<const uint8_t*> lines;
	while(text_begin != text_end) {
		auto line_end = std::find(text_begin,text_end,0);
		if (line_end==text_end) break; // the text is not ending with eol, the last line is skipped
		lines.push_back(text_begin);
		text_begin = line_end + 1;
	}
	lines.push_back(text_begin);

	// initialize
	std::vector<pair64_t> a(lines.size()-1),aa(lines.size() - 1);
	uint64_t k=0;for(auto &x:a) x.u = x.v = k++;
	uint8_t *bwt0_begin = bwt_begin = bwt_begin + Tlen;
	
	// core loop
	long Blen = 0;	
	for (long i = 0; !a.empty(); ++i) {
		long pre = 0;
		const uint8_t *end = bwt0_begin + Blen;
		Blen += a.size();
		bwt_begin -= a.size();
		const uint8_t *p = bwt0_begin;
		uint8_t *q = bwt_begin;
		
		std::array<uint64_t,256> ac,mc,mc2,b;
		mc.fill(0);mc2.fill(0);
		
		auto n = a.begin();
		for (auto k = 0; k < a.size(); ++k) {
			pair64_t& u = a[k];
			u.w = lines[u.v+1]-2-i >= lines[u.v]? *(lines[u.v+1]-2-i) : 0; // symbol to insert
			for (long l = 0; l != u.u - pre; ++l) // copy ($u->u - $pre - 1) symbols from bwt0 to bwt
				++mc[*p], *q++ = *p++; // $mc: marginal counts of all processed symbols
			*q++ = u.w;
			pre = u.u + 1; u.u = mc[u.w]++;
			if (u.w) *(n++) = u, ++mc2[u.w]; // $mc2: marginal counts of the current column
		}
		a.resize(std::distance(a.begin(),n));
		aa.resize(a.size());
		
		std::copy(p,end,q); // copy the rest of $bwt0 to $bwt
		while(p < end) ++mc[*(p++)];
		ac[0] = 0;for(int c = 1; c != ac.size(); ++c) ac[c] = ac[c-1] + mc[c-1]; // accumulative count
		for(auto &x:a) x.u += ac[x.w] + a.size(); // compute positions for the next round

		// stable counting sort ($a[k].v&0xff); also possible with an in-place non-stable radix sort, which is slower
		b[0] = 0;for (int c = 1; c != b.size(); ++c) b[c] = b[c-1] + mc2[c-1];
		for(auto x:a) aa[b[x.w]++] = x; // this works because $a is already partially sorted
		
		std::swap(a,aa); // $aa now becomes $a
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
