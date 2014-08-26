
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <memory>
#include <getopt.h>
#include <cinttypes>
#include <cinttypes>

#include <thread>
#include <mutex>
#include <condition_variable>

#include <fm_index.h>
#include <algo.h>



//
// Define alphabet
//
typedef std::string dna_string;
const dna_string alphabet("$ACGT");
inline uint8_t encode(char c) {return alphabet.find(c);}
inline char decode(uint8_t c) {return alphabet[c];}
inline uint8_t complement(uint8_t c) {return 5-c;}

typedef bwt::fm_index<5> dna_index;
typedef std::vector<dna_index> dna_indices;



//
// Getopt
//
struct args_t {
  std::vector<std::string> bwtFiles;
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
	"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n";

	enum { OPT_HELP = 1 };
	static const struct option longopts[] = {
    { "kmer-size",             required_argument, NULL, 'k' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
	};
	args_t args;
	

  for (char c; (c = getopt_long(argc, argv, "d:k:x:", longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
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



dna_indices bwts;
std::stack< stack_elt_t > stack;
std::mutex mtx,io_mtx;
std::condition_variable cv;
unsigned int num_working_thread = 0;


// extract all canonical kmers of a bwt by performing a backward depth-first-search
void traverse_kmer(unsigned int k) {
	{
		std::unique_lock<std::mutex> lck(mtx);
		++num_working_thread;
	}
	
	while(true) {
		stack_elt_t top;
		{// pop one element from the stack
	    std::unique_lock<std::mutex> lck(mtx);
	    --num_working_thread;
	    while (stack.empty() && num_working_thread>0) cv.wait(lck);
	    if (stack.empty()) break;
	    ++num_working_thread;
	    top = stack.top();
	    stack.pop();
		}
    
    for(size_t i = 1; i < alphabet.size(); ++i) {
  		if (top.lb[i]<top.ub[i]) {
        stack_elt_t e = top;
        e.path.push_back(i);
        if (e.path.length()>=k) {
        	// extract forward an reverse sequence from the path
					std::string fwd(e.path);
					std::reverse(fwd.begin(),fwd.end());
					std::string rev(e.path);
					std::transform(rev.begin(),rev.end(),rev.begin(),complement);
					
					// count number of occurence of the reverse complement
					dna_index::alpha_count64 lb,ub;
					alpha_range(bwts[0],lb,ub);
					for(size_t i=rev.size()-1;i>1;--i) bwt::extend_lhs(bwts[0],lb,ub,rev[i]);
					uint64_t rev_count = ub[rev.front()]>lb[rev.front()]?ub[rev.front()]-lb[rev.front()]:0;
					
					// 
					uint64_t fwd_count = e.ub[i]>e.lb[i]?e.ub[i]-e.lb[i]:0;
					
					// output the counts 
          if (fwd<=rev) {
          	std::transform(fwd.begin(),fwd.end(),fwd.begin(),decode);
						std::unique_lock<std::mutex> lck(io_mtx);
						std::cout << fwd << '\t' << fwd_count << '\t' << rev_count << std::endl;
          } else {
          	std::transform(rev.begin(),rev.end(),rev.begin(),decode);
						std::unique_lock<std::mutex> lck(io_mtx);
						std::cout << rev << '\t' << rev_count << '\t' << fwd_count << std::endl;          	
          }
        } else {
          bwt::extend_lhs(bwts[0],e.lb,e.ub,i);
          std::unique_lock<std::mutex> lck(mtx);
          stack.push(e);
          cv.notify_one();
        }
  		}
    }
	}
	cv.notify_all();
}



//
// Main
//

int main(int argc, char* argv[]) {
	try {
    // parse command line arguments
    args_t args = parseKmerCountOptions(argc,argv);
		
		// load bwt from files
    for(const auto& filename:args.bwtFiles) {
      bwts.push_back(dna_index(bwt::read_rle_bwt(filename)));
    }
    
    // intialize kmer traversal
		stack.push(stack_elt_t());
		bwt::alpha_range(bwts[0],stack.top().lb,stack.top().ub);
		
		// launch the threads and wait for the end
    std::vector<std::thread> threads;
    for(auto i:{0,1,2,3}) threads.push_back(std::thread(traverse_kmer,args.kmerLength));
    for(auto& t:threads) t.join();
    
	} catch (std::exception e) {
			std::cerr << e.what() << std::endl;
	};
  return 0;
}



