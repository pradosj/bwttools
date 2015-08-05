#include <Rcpp.h>
#include <array>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


struct bcr_iterator_t {
	const uint8_t* last;
	uint64_t n;
};



void bcr(const uint8_t* text_begin, const uint8_t* text_end,uint8_t* bwt_begin) {
	
  // search for lines in the text, and fill last_char[] with a pointer onto the last character of each line
  std::vector<bcr_iterator_t> last_char;
  while(text_begin != text_end) {
  	auto line_end = std::find(text_begin,text_end,0);
  	if (line_end==text_end) break; // the text is not ending with eol, the last line is skipped
  	last_char.push_back({line_end-1,last_char.size()});
  	text_begin = line_end + 1;
  }
  
  // core loop
  auto bwt_end = bwt_begin + std::distance(text_begin,text_end);
  while(!last_char.empty()) {
  	uint64_t pre = 0;
  	bwt_begin = bwt_end - last_char.size();

  	std::array<uint64_t,256> mc;mc.fill(0);
    for (auto &x:last_char) {
    	Rprintf("%d\n",x.n - pre);
    	auto q = std::copy_n(bwt_end,x.n - pre,bwt_begin);
    	std::for_each(bwt_begin,q,[&mc](uint8_t p) {++mc[p];});
    	*q++ = *x.last;
    	pre = x.n + 1;
    	x.n = mc[*x.last]++;
    	if (*x.last) {
    		
    		--x.last;
    	}
    }
		break;
/*    
    long l, pre, ac[256], mc[256], mc2[256];
  	pair64_t *b[256], *aa;
  	for (c = 0; c != 256; ++c) mc[c] = mc2[c] = 0;
  	end = B0 + Blen; Blen += n0; B -= n0;
  	for (n = k = 0, p = B0, q = B, pre = 0; k < n0; ++k) {
  		pair64_t *u = &a[k];
  		c = P[(u->v>>8) + 1] - 2 - i >= P[u->v>>8]? *(P[(u->v>>8) + 1] - 2 - i) : 0; // symbol to insert
  		u->v = (u->v&~0xffULL) | c;
  		for (l = 0; l != u->u - pre; ++l) // copy ($u->u - $pre - 1) symbols from B0 to B
  			++mc[*p], *q++ = *p++; // $mc: marginal counts of all processed symbols
  		*q++ = c;
  		pre = u->u + 1; u->u = mc[c]++;
  		if (c) a[n++] = a[k], ++mc2[c]; // $mc2: marginal counts of the current column
  	}
  	
  	for(auto x:std::vector<pair64_t>(a,a+n)) Rprintf("(%d,%d),",x.u,x.v);Rprintf("\n");
  	
  	while (p < end) ++mc[*p], *q++ = *p++; // copy the rest of $B0 to $B
  	for (c = 1, ac[0] = 0; c != 256; ++c) ac[c] = ac[c-1] + mc[c-1]; // accumulative count
  	for (k = 0; k < n; ++k) a[k].u += ac[a[k].v&0xff] + n; // compute positions for the next round
  	// stable counting sort ($a[k].v&0xff); also possible with an in-place non-stable radix sort, which is slower
  	aa = (pair64_t *) malloc(sizeof(pair64_t) * n);
  	for (c = 1, b[0] = aa; c != 256; ++c) b[c] = b[c-1] + mc2[c-1];
  	for (k = 0; k < n; ++k) *b[a[k].v&0xff]++ = a[k]; // this works because $a is already partially sorted
  	free(a); a = aa; // $aa now becomes $a
  	B0 = B; n0 = n;
*/
	}
  
  
}

// [[Rcpp::export]]
CharacterVector testBcr() {
  const uint8_t* str = (const uint8_t*) "BANANA\0BANANA\0";
  uint8_t bwt[14];
  bcr(str,str+14,bwt);
  return Rcpp::CharacterVector::create((const char*) bwt);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
testBcr()
*/
