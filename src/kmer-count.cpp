
#include <iostream>
#include <sstream>
#include <stack>
#include <memory>
#include <getopt.h>

#include <bwt.h>
#include <BWTAlgorithms.h>



const std::string alphabet("$ACGT");

inline std::string reverse(const std::string& str) {
	std::string rstr(str);
  std::reverse(rstr.begin(),rstr.end());
  return(rstr);
}

inline std::string complement(std::string str) {
    auto complement_bp = [](char bp) -> char {
    	switch(bp) {
        case 'A':return 'T';
        case 'C':return 'G';
        case 'G':return 'C';
        case 'T':return 'A';
        case 'N':return 'N';
        default:
            assert(false && "Unknown base!");
            return 'N';
			};
    };
    std::transform(str.begin(), str.end(), str.begin(), complement_bp);
    return(str);
}



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
    uint8_t bp;
    bwt::interval range;
    stack_elt_t(size_t strsz,char bp,const bwt::interval& range):str_sz(strsz),bp(bp),range(range){}
};



// extract all canonical kmers of a bwt by performing a backward depth-first-search
void traverse_kmer(size_t n, const bwt::fm_index** pBWTs, unsigned int k) {
    std::stack< stack_elt_t > stack;
    std::string str; // string storing the current path

    // Intitialize the search with root elements
    for(uint8_t i = 1; i < alphabet.size(); ++i) {
        stack_elt_t e(str.size(),i,pBWTs[0]->initInterval(i));
        if (!e.range.empty()) stack.push(e);
    }

    // Perform the kmer search
    while(!stack.empty()) {
        // pop an element from the stack, and update the string path accordingly
        stack_elt_t top = stack.top();
        stack.pop();
        str.resize(top.str_sz);
        str.push_back(alphabet[top.bp]);
        if (str.length()>=k) {
            // we found a kmer, retreive the number of occurence
            std::string seq(reverse(str));
            int64_t seq_count = top.range.size();

            // look for the count of its reverse complement
            std::string seq_rc(complement(str));
            bwt::interval range_rc = BWTAlgorithms::findInterval(*pBWTs[0],seq_rc);
            int64_t seq_rc_count = range_rc.empty()?0:range_rc.size();

            // print the current kmer if canonical
            if (seq<seq_rc) {
                std::cout << seq << '\t' << seq_count << '\t' << seq_rc_count;
                for(size_t i=1;i<n;i++) {
                    std::cout << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(*pBWTs[i],seq).size(),0) << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(*pBWTs[i],seq_rc).size(),0);
                }
                std::cout << std::endl;
            } else if (seq_rc_count<=0) {
                // the current kmer is not canonical, but the reverse complement doesn't exists
                // so print it now as it will never be traversed by the searching algorithm
                std::cout << seq_rc << '\t' << seq_rc_count << '\t' << seq_count;
                for(size_t i=1;i<n;i++) {
                    std::cout << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(*pBWTs[i],seq_rc).size(),0) << '\t' << std::max<int64_t>(BWTAlgorithms::findInterval(*pBWTs[i],seq).size(),0);
                }
                std::cout << std::endl;
            }
        } else {
            // not yet a suffix of size k, push next candidates
            for(size_t i = 1; i < alphabet.size(); ++i) {
                stack_elt_t e(str.size(),i,top.range);
                pBWTs[0]->updateInterval(e.range,e.bp);
                if (!e.range.empty()) stack.push(e);
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
    const bwt::fm_index** pBWTs = new const bwt::fm_index*[opt::bwtFiles.size()];
    for(size_t i=0;i<opt::bwtFiles.size();++i) {
        pBWTs[i] = new bwt::fm_index(opt::bwtFiles[i],opt::sampleRate);
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



