

OBJECTS := Occurrence.o BWTReader.o bwt.o BWTAlgorithms.o


all:kmer-count


kmer-count:kmer-count.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^


clean:
	rm -f kmer-count *.o


