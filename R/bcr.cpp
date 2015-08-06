#include <Rcpp.h>
#include <array>
#include <algorithm>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

typedef struct {
	uint64_t u,v;
	uint8_t w;
} pair64_t;



void bcr(const uint8_t* text_begin, const uint8_t* text_end,uint8_t* bwt_begin) {
	long Tlen = std::distance(text_begin,text_end);
	
	// split $T into short strings at sentinels
	if (text_begin == text_end) return;
	std::vector<const uint8_t*> P;
	while(text_begin != text_end) {
		auto line_end = std::find(text_begin,text_end,0);
		if (line_end==text_end) break; // the text is not ending with eol, the last line is skipped
		P.push_back(text_begin);
		text_begin = line_end + 1;
	}
	P.push_back(text_begin);
	long n = P.size() - 1;
	
	// initialize
	std::vector<pair64_t> a(n),aa(n);
	for (long k = 0; k < n; ++k) a[k].u = k, a[k].v = k;
	uint8_t *B0 = bwt_begin = bwt_begin + Tlen;
	
	// core loop
	long n0 = n;
	long Blen = 0;	
	for (long i = 0; n0; ++i) {
		long pre = 0;
		const uint8_t *end = B0 + Blen;
		Blen += n0;
		bwt_begin -= n0;
		const uint8_t *p = B0;
		uint8_t *q = bwt_begin;
		
		std::array<long,256> ac,mc,mc2;
		mc.fill(0);mc2.fill(0);
		
		for (long k = n = 0; k < n0; ++k) {
			pair64_t *u = &a[k];
			int c = P[u->v+1]-2-i >= P[u->v]? *(P[u->v+1]-2-i) : 0; // symbol to insert
			u->w = c;
			for (long l = 0; l != u->u - pre; ++l) // copy ($u->u - $pre - 1) symbols from B0 to B
				++mc[*p], *q++ = *p++; // $mc: marginal counts of all processed symbols
			*q++ = c;
			pre = u->u + 1; u->u = mc[c]++;
			if (c) a[n++] = a[k], ++mc2[c]; // $mc2: marginal counts of the current column
		}
		std::copy(p,end,q); // copy the rest of $B0 to $B
		while(p < end) ++mc[*(p++)];
		ac[0] = 0;for (int c = 1; c != 256; ++c) ac[c] = ac[c-1] + mc[c-1]; // accumulative count
		for(long k = 0; k < n; ++k) a[k].u += ac[a[k].w] + n; // compute positions for the next round

		// stable counting sort ($a[k].v&0xff); also possible with an in-place non-stable radix sort, which is slower
		std::array<uint64_t,256> b;
		b[0] = 0;for (int c = 1; c != 256; ++c) b[c] = b[c-1] + mc2[c-1];
		for (long k = 0; k < n; ++k) aa[b[a[k].w]++] = a[k]; // this works because $a is already partially sorted
		
		std::swap(a,aa); // $aa now becomes $a
		B0 = bwt_begin; n0 = n;
	}
}

// [[Rcpp::export]]
CharacterVector testBcr() {
  const uint8_t* str = (const uint8_t*) "BANANA\0BANANA\0";
  uint8_t bwt[14];
  bcr(str,str+14,bwt);
  return Rcpp::CharacterVector::create((const char*) bwt);
}



/*** R
testBcr()
*/
