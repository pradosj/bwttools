
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <memory>
#include <getopt.h>
#include <cinttypes>

#include <fm_index.h>
#include <algo.h>



typedef std::string dna_string;
const dna_string alphabet("$ACGT");
auto encode = [&](char c)->uint8_t{return alphabet.find(c);};
auto decode = [&](uint8_t c)->char{return alphabet[c];};

typedef bwt::fm_index<5> dna_index;
typedef std::vector< std::unique_ptr<dna_index> > dna_indices;



//
// Getopt
//
struct args_t {
    std::vector<std::string> bwtFiles;
    int sampleRate = 128;
    int kmerLength = 27;
};

args_t parseKmerCountOptions(int argc, char* argv[]) {
		static const char* usage_message =
		"Usage: kmer-count [OPTION] src.bwt [test1.bwt] [test2.bwt]\n"
		"Generate a table of the k-mers in src.bwt, and optionaly count the number of time they appears in testX.bwt.\n"
		"Output on stdout the canonical kmers and their counts on forward and reverse strand\n"
		"\n"
		"      --help                           display this help and exit\n"
		"      --version                        display program version\n"
		"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
		"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
		"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n";

		enum { OPT_HELP = 1 };
		static const struct option longopts[] = {
		    { "sample-rate",           required_argument, NULL, 'd' },
		    { "kmer-size",             required_argument, NULL, 'k' },
		    { "help",                  no_argument,       NULL, OPT_HELP },
		    { NULL, 0, NULL, 0 }
		};
		args_t args;
		

    for (char c; (c = getopt_long(argc, argv, "d:k:x:", longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'd': arg >> args.sampleRate; break;
            case 'k': arg >> args.kmerLength; break;
            case OPT_HELP:
                std::cout << usage_message;
                exit(EXIT_SUCCESS);
        }
    }

    if(args.kmerLength <= 0 || args.kmerLength % 2 == 0) {
        std::cerr << "kmer-count: invalid kmer length: " << args.kmerLength << ", must be greater than zero and odd\n";
        std::cout << "\n" << usage_message;
        exit(EXIT_FAILURE);
    }

    for(;optind<argc;++optind) {
        args.bwtFiles.push_back(argv[optind]);
    }

    if (args.bwtFiles.size() < 1) {
        std::cerr << "kmer-count: missing arguments\n";
        std::cout << "\n" << usage_message;
        exit(EXIT_FAILURE);
    }
    
    return args;
}



//
// BWT Traversal algorithm
//

// Stack structure used in the depth first search of kmers
struct stack_elt_t {
    dna_string path;
    dna_index::alpha_count64 lb;
    dna_index::alpha_count64 ub;
};



// extract all canonical kmers of a bwt by performing a backward depth-first-search
void traverse_kmer(dna_indices& bwts, unsigned int k) {
    std::stack< stack_elt_t > stack;
		stack.push(stack_elt_t());
		bwt::init_char_range(*bwts[0],stack.top().lb,stack.top().ub);

    // Perform the kmer search
    while(!stack.empty()) {
        // pop an element from the stack, and update the string path accordingly
        stack_elt_t top = stack.top();
        stack.pop();
        for(size_t i = 1; i < alphabet.size(); ++i) {
        		if (top.lb[i]<top.ub[i]) {
                stack_elt_t e = top;
                e.path.push_back(i);
				        if (e.path.length()>=k) {
				            std::reverse(e.path.begin(),e.path.end());
				            std::transform(e.path.begin(),e.path.end(),e.path.begin(),decode);
				            std::cout << e.path << '\t' << (e.ub[i]-e.lb[i]) << std::endl;
				        } else {
		                bwt::extend_lhs(*bwts[0],e.lb,e.ub,i);
		                stack.push(e);
				        }
        		}
        }
    }
}




//
// Main
//

int main(int argc, char* argv[]) {
	try {
    // parse command line arguments
    args_t args = parseKmerCountOptions(argc,argv);
		
		// load bwt from files
    dna_indices bwts;
    for(auto filename:args.bwtFiles) {
    		std::ifstream f(filename.c_str(),std::ios::binary);
        bwts.push_back(std::unique_ptr<dna_index>(new dna_index(f)));
    }
    
    // traverse the kmers
    traverse_kmer(bwts,args.kmerLength);
	} catch (std::exception e) {
			std::cerr << e.what() << std::endl;
	};
  return 0;
}



