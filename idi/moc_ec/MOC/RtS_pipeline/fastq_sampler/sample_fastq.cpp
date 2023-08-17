#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <random>
#include <cmath>
#include <boost/program_options.hpp>

#include "fastq_reader.hpp"
#include "fastq_writer.hpp"

namespace po = boost::program_options;

class sample_fastqc {

    public:
        void print_help();
        bool parse_args(int argc, char* argv[]); 
        void initialize();
        void main_func();
        void free_vars();
        unsigned long get_frag_count(std::string& infile_str);

    private:
        po::options_description desc;
        std::string infile_str;
        std::string outfile_str;
        unsigned long long top_seed;
        unsigned long interval;
        unsigned long sample_count;
        unsigned long infile_frag_count;
        long read_limit;
        double sample_p;

};

void sample_fastqc::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: sample_fastq --infile fastq --outfile fastq"
        " --top_seed <top_seed_number>  --sample_count <sample count>"
        "\n\n";
}

bool sample_fastqc::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input fastq.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output fastq.")
        ("top_seed,t", po::value(&top_seed),
            "Top seed for sampling.")
        ("sample_count,s", po::value(&sample_count),
            "Number of fragments to sample.")
        ("interval", po::value(&interval)->default_value(1000000),
            "Interval for shuffle.")
        ("read_limit,r", po::value(&read_limit)->default_value(-1),
            "Pair of reads to process")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);


    if (vm.count("help")) {
        return 0;
    } else {
    }

    if (vm.count("infile")) {
        std::cout << "Infile is set to: " << infile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: infile is not set.\n";
    }

    if (vm.count("outfile")) {
        std::cout << "Outfile is set to " << outfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outfile is not set.\n";
    }

 
    if (vm.count("top_seed")) {
        std::cout << "top_seed is set to: " << top_seed << "\n";
    } else {
        all_set = false;
        std::cout << "Error: top_seed is not set.\n";
    }

    if (vm.count("sample_count")) {
        std::cout << "sample_count is set to: " << sample_count << "\n";
    } else {
        std::cout << "Error: sample_count is not set.\n";
    }

    std::cout << "Interval is set to: " << interval << "\n";

    return all_set;

}


// Returns the number of reads in a file
unsigned long sample_fastqc::get_frag_count(std::string& infile_str) {

    fastq_reader fro(infile_str);

    std::string lword;
    unsigned long lcount = 0;
    unsigned long lfrag_count = 0;
    while (fro.getline(lword)) {
        lcount++;
    }

    lfrag_count = lcount / 4;

    return lfrag_count;

}



void sample_fastqc::initialize() {
    // get the total number of reads of the infile.
    
    // Initialize the number of reads in the input file.
    infile_frag_count = get_frag_count(infile_str);

    std::cout << "Infile_frag_count: " << infile_frag_count << "\n";
    unsigned long final_limit = infile_frag_count;
    if (read_limit > 0) {
        final_limit = read_limit;
    }
    sample_p = (sample_count * 100.0) / final_limit;

}

void sample_fastqc::main_func() {

    std::vector<unsigned long> seq_vec;
    std::set<unsigned long> write_set;
    std::mt19937_64 r_engine(top_seed);
    unsigned long long top_seed_2 = top_seed + 1;
    std::mt19937_64 r_engine_2(top_seed_2);


    std::uniform_real_distribution<double> ldist(0.0, 1.0);

    fastq_reader lreader(infile_str);
    fastq_writer lwriter(outfile_str);

    // Total number of reads loaded/read till this point.
    unsigned long total_frag_ind = 0;
   
    unsigned long frag_left = infile_frag_count; 
    unsigned long write_ind = 0;
    unsigned long frag_this_iter = 0;
    unsigned long write_lim = 0;
    unsigned long total_frag_written = 0;

    std::string lword1;
    std::string lword2;
    std::string lword3;
    std::string lword4;

    double frag_total = 0.0;

    while(true == lreader.getline(lword1)) {
        // Read other three lines
        lreader.getline(lword2);
        lreader.getline(lword3);
        lreader.getline(lword4);

        

        if (total_frag_ind % interval == 0) {

            // Do the assignment
            if (frag_left >= interval) {
                frag_this_iter = interval;
                frag_left = frag_left - interval;

                double write_lim_f = (frag_this_iter * sample_p) / 100.0; 
                write_lim = (long) write_lim_f;

                frag_total += write_lim_f - write_lim;
                if (frag_total >=1 ) {
                    double to_add = (long) frag_total;
                    
                    write_lim += to_add;
                    frag_total -= to_add;
                }

            } else {
                frag_this_iter = frag_left;
                frag_left = 0;
                unsigned long target_frag = sample_count - total_frag_written;
                write_lim = target_frag;
            }

            seq_vec.assign(frag_this_iter, 0);
            write_set.clear();
            for (unsigned long j = 0; j < frag_this_iter; j++) {
                unsigned long l_seq_val = total_frag_ind + j;
                seq_vec[j] = l_seq_val;       
            }

            // Do the shuffle
            shuffle(seq_vec.begin(), seq_vec.end(), r_engine);
            
            write_ind = 0;    
            
            for (unsigned int j = 0; j < write_lim ; j++) {
                write_set.insert(seq_vec[j]);
            }
        }

        // Write if appropriate
        if (write_ind < write_lim) {
            if (write_set.find(total_frag_ind) != write_set.end()) {

                // Write the four lines
                lwriter.putline(lword1);
                lwriter.putline(lword2);
                lwriter.putline(lword3);
                lwriter.putline(lword4);
                
                write_ind++;
                total_frag_written++;
            }
        }

        total_frag_ind++; 
        if (total_frag_ind == read_limit) {
            break;
        }

    }

    std::cout << "Total frag written: " << total_frag_written << "\n";

}

void sample_fastqc::free_vars() {

}

int main(int argc, char** argv) {
    sample_fastqc sfc;
    bool all_set = true;

    try {
        all_set = sfc.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        return 0;
    } 

    if (!all_set) {
        sfc.print_help();
        return 0;
    }

    try {
        sfc.initialize();
        sfc.main_func();
        sfc.free_vars();
    } catch(const std::runtime_error& e) {
        std::cerr << "error: "  << e.what() << "\n";
        return 1;
    }

    return 0;
}


