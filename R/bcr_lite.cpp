#include <Rcpp.h>
#include <array>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

typedef struct {
	uint64_t u,v;
} pair64_t;


/**
* Append $T to existing BWT $B.
*
* @param Blen    length of the existing BWT
* @param B       existing BWT; set to NULL if non-existing
* @param Tlen    length of input string
* @param T       input string; '\0' represents a sentinel
*
* @return  the new BWT string
*/
uint8_t *bcr_lite(long Blen, uint8_t *B, long Tlen, const uint8_t *T)
{
	long i, k, n, max, n0;
	uint8_t *p, *q, *B0;
	const uint8_t *end, **P = 0;
	pair64_t *a;
	int c;
	// split $T into short strings at sentinels
	if (T == 0 || Tlen == 0) return B;
	for (p = q = (uint8_t*)T, end = T + Tlen, n = max = 0; p != end; ++p) {
		if (*p) continue;
		if (n == max) {
			max = max? max<<1 : 256;
			P = (const uint8_t**) realloc(P, max * sizeof(void*));
		}
		P[n++] = q, q = p + 1;
	}
	P = (const uint8_t**) realloc(P, (n + 1) * sizeof(void*));
	P[n] = q;
	
	// initialize
	for (p = B, end = B + Blen, i = 0; p < end; ++p) i += (*p == 0); // count # of sentinels
	a = (pair64_t*) malloc(sizeof(pair64_t) * n);
	for (k = 0; k < n; ++k) a[k].u = k + i, a[k].v = k<<8;
	B = (uint8_t*) realloc(B, Blen + Tlen);
	memmove(B + Tlen, B, Blen); // finished BWT is always placed at the end of $B
	B = B0 = B + Tlen;
	
	// core loop
	for (i = 0, n0 = n; n0; ++i) {
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
	}
	free(P); free(a);
	return B;
}


// [[Rcpp::export]]
CharacterVector testBcr() {
	const int len = 12;
	const uint8_t* str = (const uint8_t*) "BANA\0BANANA\0ANANAS\0RANANA\0";
	const char* bwt = (const char*) bcr_lite(0,(uint8_t*) 0,len,str);
	std::string BWT((const char*) bwt,len);
	std::replace(BWT.begin(),BWT.end(),'\0','$');
	return Rcpp::CharacterVector::create(BWT);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
testBcr()
*/
