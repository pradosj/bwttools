
kmer-count:kmer-count.cpp libbwt/libbwt.a
	$(CXX) $(CXXFLAGS) -o $@ -Ilibbwt $^

libbwt/libbwt.a:
	$(MAKE) -C $(@D)

clean:
	rm -f kmer-count


