
export CXXFLAGS += -std=c++11 -O3
#export CXXFLAGS += -DNDEBUG

all:test kmer-count

test:test.cpp
	$(CXX) $(CXXFLAGS) -o $@ -Ilibbwt $^

kmer-count:kmer-count.cpp
	$(CXX) $(CXXFLAGS) -o $@ -Ilibbwt $^

libbwt/libbwt.a:
	$(MAKE) -C $(@D)

clean:
	$(MAKE) -C libbwt clean
	rm -f test kmer-count


