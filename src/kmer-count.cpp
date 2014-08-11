
#include <iostream>
#include <sstream>
#include <stack>
#include <memory>
#include <getopt.h>
#include "bwt.h"
#include "BWTInterval.h"
#include "BWTAlgorithms.h"


template<> const size_t AlphaCount16::maxValue = std::numeric_limits<uint16_t>::max();
template<> const size_t AlphaCount64::maxValue = std::numeric_limits<uint64_t>::max();

//
// Getopt
//
static const char *KMERCOUNT_USAGE_MESSAGE =
"Usage: kmer-count [OPTION] src.bwt [test1.bwt] [test2.bwt]\n"
"Generate a table of the k-mers in src.bwt, and optionaly count the number of time they appears in testX.bwt.\n"
"Output on stdout the canonical kmers and their counts on forward and reverse strand\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display program version\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n";


namespace opt {
    static std::vector<std::string> bwtFiles;
    static int sampleRate = 128;
    static int kmerLength = 27;
}

static const char* shortopts = "d:k:x:";
enum { OPT_HELP = 1 };
static const struct option longopts[] = {
    { "sample-rate",           required_argument, NULL, 'd' },
    { "kmer-size",             required_argument, NULL, 'k' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};


void parseKmerCountOptions(int argc, char* argv[]) {

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case OPT_HELP:
                std::cout << KMERCOUNT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(opt::kmerLength <= 0 || opt::kmerLength % 2 == 0) {
        std::cerr << "kmer-count: invalid kmer length: " << opt::kmerLength << ", must be greater than zero and odd\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    for(;optind<argc;++optind) {
        opt::bwtFiles.push_back(argv[optind]);
    }

    if (opt::bwtFiles.size() < 1) {
        std::cerr << "kmer-count: missing arguments\n";
        std::cout << "\n" << KMERCOUNT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}



//
// BWT Traversal algorithm
//

// Stack structure used in the depth first search of kmers
struct stack_elt_t {
    size_t str_sz;
    char bp;
    BWTInterval range;
    stack_elt_t(size_t strsz,char bp):str_sz(strsz),bp(bp){}
    stack_elt_t(size_t strsz,char bp,const BWTInterval& range):str_sz(strsz),bp(bp),range(range){}
};



// extract all canonical kmers of a bwt by performing a backward depth-first-search
void traverse_kmer(size_t n, const bwt** pBWTs, unsigned int k) {
    std::stack< stack_elt_t > stack;
    std::string str; // string storing the current path

    // Intitialize the search with root elements
    for(size_t i = 1; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
        stack_elt_t e(str.size(),BWT_ALPHABET::RANK_ALPHABET[i]);
        BWTAlgorithms::initInterval(e.range,e.bp,pBWTs[0]);
        if (e.range.isValid()) stack.push(e);
    }

    // Perform the kmer search
    while(!stack.empty()) {
        // pop an element from the stack, and update the string path accordingly
        stack_elt_t top = stack.top();
        stack.pop();
        str.resize(top.str_sz);
        str.push_back(top.bp);
        if (str.length()>=k) {
            // we found a kmer, retreive the number of occurence
            std::string seq(reverse(str));
            int64_t seq_count = top.range.size();

            // look for the count of its reverse complement
            std::string seq_rc(complement(str));
            BWTInterval range_rc;// = BWTAlgorithms::findInterval(pBWTs[0],seq_rc);
            int64_t seq_rc_count = range_rc.isValid()?range_rc.size():0;

            // print the current kmer if canonical
            if (seq<seq_rc) {
                std::cout << seq << '\t' << seq_count << '\t' << seq_rc_count;
                for(size_t i=1;i<n;i++) {
                    std::cout << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(pBWTs[i],seq).size(),0) << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(pBWTs[i],seq_rc).size(),0);
                }
                std::cout << std::endl;
            } else if (seq_rc_count<=0) {
                // the current kmer is not canonical, but the reverse complement doesn't exists
                // so print it now as it will never be traversed by the searching algorithm
                std::cout << seq_rc << '\t' << seq_rc_count << '\t' << seq_count;
                for(size_t i=1;i<n;i++) {
                    std::cout << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(pBWTs[i],seq_rc).size(),0) << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(pBWTs[i],seq).size(),0);
                }
                std::cout << std::endl;
            }
        } else {
            // not yet a suffix of size k, push next candidates
            for(size_t i = 1; i < BWT_ALPHABET::ALPHABET_SIZE; ++i) {
                stack_elt_t e(str.size(),BWT_ALPHABET::RANK_ALPHABET[i],top.range);
                BWTAlgorithms::updateInterval(e.range,e.bp,pBWTs[0]);
                if (e.range.isValid()) stack.push(e);
            }
        }
    }
}




//
// Main
//
int main(int argc, char* argv[]) {
    // parse command line arguments
    parseKmerCountOptions(argc,argv);

    // allocate BWT objects
    const bwt** pBWTs = new const bwt*[opt::bwtFiles.size()];
    for(size_t i=0;i<opt::bwtFiles.size();++i) {
        pBWTs[i] = new bwt(opt::bwtFiles[i],opt::sampleRate);
    }

    // run kmer search
    traverse_kmer(opt::bwtFiles.size(),pBWTs,opt::kmerLength);

    // clean memory
    for(size_t i=0;i<opt::bwtFiles.size();++i) {
        delete pBWTs[i];
    }
    delete[] pBWTs;
    return 0;
}



