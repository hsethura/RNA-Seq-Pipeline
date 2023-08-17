#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <random>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <boost/program_options.hpp>


namespace po = boost::program_options;

// This class would generate a table with list of top seeds along with the
// number of fragments to sample.


class args_c {
    public:
        po::options_description desc;
        std::string infile_str;
        std::string outfile_str;
        long main_seed;
        void print_help();
        bool parse_args(int argc, char* argv[]);

};

class seed_fgenc {

    public:
        seed_fgenc(args_c args_o);
        void initialize();
        void main_func();
        void free_vars();

    private:
        std::string infile_str;
        std::string outfile_str;
        long main_seed;
        std::vector<double> sample_p_vec;
};

seed_fgenc::seed_fgenc(args_c args_o) {
    infile_str = args_o.infile_str;
    outfile_str = args_o.outfile_str;
    main_seed = args_o.main_seed;
}


void args_c::print_help() {
        std::cout << desc << "\n";
    std::cout << "Usage: seed_fgen --infile <infile> --outfile <seed_table>"
        "[ --main_seed <main_seed_number>"
        "\n\n";
}

bool args_c::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output file.")
        ("main_seed,m", po::value(&main_seed)->default_value(12345),
            "Man seed for generating the top seeds.")
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

   if (vm.count("main_seed")) {
        std::cout << "main_seed is set to: " << main_seed << "\n";
    } else {
        all_set = false;
        std::cout << "Error: main_seed is not set.\n";
    }

    return all_set;
}

void seed_fgenc::initialize() {
    
}

void seed_fgenc::main_func() {
    std::mt19937_64 r_engine(main_seed);
    std::ifstream inf(infile_str);
    std::ofstream outf(outfile_str);

    // This would be simple, just read one line from the infile, generate a 
    // random seed and append at the end of each line.

    std::string line1;

    unsigned int lcount = 0;
    while(std::getline(inf, line1)) {

        // generate one random number
        if (0 == lcount) {
            outf << line1 << "\ttop_seed\n";
        } else {
            unsigned long long lrand_num = r_engine();
            std::string outline = line1 + "\t" + std::to_string(lrand_num);
            outf << outline << "\n";
        }
        lcount++;
    }
}


void seed_fgenc::free_vars() {

}

int main(int argc, char** argv) {

    args_c args_o;
    bool all_set = true;

    try {
        all_set = args_o.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        return 0;
    }

    if (!all_set) {
        args_o.print_help();
        return 0;
    }

    seed_fgenc sfgc(args_o);
    try {
        sfgc.initialize();
        sfgc.main_func();
        sfgc.free_vars();
    } catch(const std::runtime_error& e) {
        std::cerr << "error: "  << e.what() << "\n";
        return 1;
    }

    return 0;
}
