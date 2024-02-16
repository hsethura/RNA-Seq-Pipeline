#include <iostream>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iomanip>

namespace po = boost::program_options;

class metrics_gen {

    public:
        
        void print_help();
        bool parse_args(int argc, char* argv[]);
        void init_gene_biotypes();
        void process_count_entry(std::string geneid, unsigned long gene_count);
        void load_gene_count();
        //void process_basic_metrics_file();
        void print_final_metrics();
        void load_basic_metrics();
        //void copy_basic_metrics_file();
        void print_feature_type(std::string feature_name, unsigned long count_s,
            unsigned long count_as, unsigned long aligned_frag, 
            std::ofstream& outfile); 
        double get_pct(unsigned long numer, unsigned long denom); 
        void load_basic_metrics_umi();
        bool use_umi;
    private:
        std::string infile_str;
        std::string metfile_str;
        std::string outfile_str;

        unsigned long ribosomal_rna_s = 0;
        unsigned long ribosomal_rna_as = 0;
        unsigned long protein_coding_s = 0;
        unsigned long protein_coding_as = 0;
        unsigned long pseudogene_s = 0;
        unsigned long pseudogene_as = 0;
        unsigned long long_noncoding_s = 0;
        unsigned long long_noncoding_as = 0;
        unsigned long short_noncoding_s = 0;
        unsigned long short_noncoding_as = 0;
        unsigned long uncategorized_s = 0;
        unsigned long uncategorized_as = 0;

        unsigned long check_s = 0;
        unsigned long check_as = 0;

        unsigned long total_count = 0;
        unsigned long sense_frag = 0;
        unsigned long antisense_frag = 0;
        unsigned long aligned_frag = 0;
        unsigned long unmapped_count = 0;
        unsigned long total_frags = 0;
        unsigned long s_drop_count = 0;
        unsigned long as_drop_count = 0;

        std::string umi_metrics_str;
        
        std::set<std::string> protein_codings;
        std::set<std::string> pseudogenes;
        std::set<std::string> long_noncodings;
        std::set<std::string> short_noncodings;
        po::options_description desc;
};


void metrics_gen::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: metrics_gen -i <infile> -m <metfile> -o <outfile>\n\n";
}

bool metrics_gen::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "Input count file.")
        ("metfile,m", po::value<std::string>(&metfile_str), "Metrics infile.")
        ("outfile,o", po::value<std::string>(&outfile_str), "Metrics outfile.")
        ("use_umi", po::bool_switch(&use_umi)->default_value(false), "Use UMI information with output")
        ("umi_metrics,u", po::value<std::string>(&umi_metrics_str), "File containing umi metrics")
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

    if (vm.count("metfile")) {
        std::cout << "Metrics infile is set to " << metfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Metrics infile is not set.\n";
    }

    if (vm.count("outfile")) {
        std::cout << "Metrics outfile is set to " << outfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Metrics outfile is not set.\n";
    }

    std::cout << "use_umi: " << std::boolalpha << use_umi << "\n";

    if (use_umi) {
        if (vm.count("umi_metrics")) {
            std::cout << "Umi metrics is set to: " << umi_metrics_str << "\n";
        } else {
            all_set = false;
            std::cout << "Umi metrics is not set\n";
        }
    }

    return all_set;
}


void metrics_gen::init_gene_biotypes() {
    protein_codings = {"IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
        "polymorphic_pseudogene", "protein_coding", "TR_C_gene", "TR_D_gene",
        "TR_J_gene", "TR_V_gene"};
    pseudogenes = {"IG_C_pseudogene", "IG_J_pseudogene", "IG_pseudogene", 
        "IG_V_pseudogene", "processed_pseudogene", "pseudogene", 
        "transcribed_processed_pseudogene", 
        "transcribed_unprocessed_pseudogene", "TR_J_pseudogene", 
        "TR_V_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "transcribed_unitary_pseudogene"};
    long_noncodings = {"3prime_overlapping_ncRNA", "antisense", "lincRNA",
        "non_coding", "processed_transcript", "sense_intronic", 
        "sense_overlapping"};
    short_noncodings = {"miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "rRNA",
        "scRNA", "snoRNA", "snRNA"};
}

void metrics_gen::process_count_entry(std::string geneid, 
        unsigned long gene_count) {
     
    // First check if the data is from sense or anti sense.
    bool is_AS = false;
    std::string l_geneid;
    if (0 == geneid.find("AS:")) {
        // Remove the AS: from the start of the string
        boost::regex re("^AS:");
        l_geneid = boost::regex_replace(geneid, re, "");
        is_AS = true;
        check_as += gene_count;
        
    } else {
        l_geneid = geneid;
        check_s += gene_count;
    }
    //std::cout << l_geneid << "\n";

    if (0 == l_geneid.find("NR_")) {
        if (is_AS) {
           ribosomal_rna_as += gene_count; 
        } else {
            ribosomal_rna_s += gene_count;
        }
    } else {
        // Get gene_biotype
        // Updated the regex to accomodate both old and new reference
        // after request from Lu.
        //boost::regex expr("^\\S+?:\\S+?:\\S+?:(\\S+?):");
        boost::regex expr("^\\S+?:\\S+?:\\S+?:(\\w+)");
        boost::smatch what;
        bool res = boost::regex_search(l_geneid, what, expr);
        //std::cout << std::boolalpha << "res: " << res << "\n";
        if (res) {
            // gbt is for gene_biotype
           
            std::string gbt = what[1];
            //std::cout << "gbt: " << gbt << " " << gene_count << "\n";
            if (protein_codings.find(gbt) != protein_codings.end()) {
                if (is_AS) {
                    protein_coding_as += gene_count;
                } else {
                    protein_coding_s += gene_count;
                } 
            } else if (pseudogenes.find(gbt) != pseudogenes.end()) {
                if (is_AS) {
                    pseudogene_as += gene_count;
                } else {
                    pseudogene_s += gene_count;
                }
            } else if (long_noncodings.find(gbt) != long_noncodings.end()) {
                if (is_AS) {
                    long_noncoding_as += gene_count;
                } else {
                    long_noncoding_s += gene_count;
                }
            } else if (short_noncodings.find(gbt) != short_noncodings.end()) {
                if (is_AS) {
                    short_noncoding_as += gene_count;
                } else {
                    short_noncoding_s += gene_count;
                }
            } else {
                if (is_AS) {
                    uncategorized_as += gene_count;
                } else {
                    uncategorized_s += gene_count;
                }
                //std::cout << "Uncategorized gbt: " << gbt << "\n";
            }
        }

    }
}

void metrics_gen::load_gene_count() {
    std::ifstream words(infile_str);
    std::string lstr;
    std::string header;
    boost::regex expr ("(\\S+)\\s+(\\d+)");
    boost::smatch what;
    if (!words.is_open()) {
        std::string lstr =  "The gene count file cannot be open!\n";
        throw std::runtime_error(lstr);
    } else {
        std::getline(words, header);
        while (std::getline(words, lstr)) {
            bool res = boost::regex_search(lstr, what, expr);
            if (res) {
                std::string geneid = what[1];
                std::string gene_count_str = what[2];
                unsigned long gene_count = std::stoul(gene_count_str);
                process_count_entry(geneid, gene_count);
            }
        }
    }
}

void metrics_gen::load_basic_metrics_umi() {
    std::ifstream  words(umi_metrics_str);
    std::string lstr;
    boost::regex expr("(\\S+)\\s+(\\S+)");
    boost::smatch what;
    // The variables to fill
   
    int lval_count = 0;
    std::getline(words, lstr);
    while(std::getline(words, lstr)) {
        bool res = boost::regex_search(lstr, what, expr);
        if (res) {
            std::string metrics_type = what[1];
            unsigned long lval = std::stoul(what[2]);
            if (0 == metrics_type.compare("s_drop_count")) {
                s_drop_count = lval;
                lval_count++;
                check_s += s_drop_count;
            } else if (0 == metrics_type.compare("as_drop_count")) {
                as_drop_count = lval;
                lval_count++;
                check_as += as_drop_count;
            }
        }
    }

    if (lval_count != 2) {
        std::string lstr = "Failed to extract six main metrics from the "
            "umi metrics file.";
        throw std::runtime_error(lstr);
    }

    words.close();


}

void metrics_gen::load_basic_metrics() {

    std::ifstream  words(metfile_str);
    std::string lstr;
    boost::regex expr("(\\S+)\\s+(\\S+)");
    boost::smatch what;
    // The variables to fill
   
    int lval_count = 0;
    std::getline(words, lstr);
    while(std::getline(words, lstr)) {
        bool res = boost::regex_search(lstr, what, expr);
        if (res) {
            std::string metrics_type = what[1];
            unsigned long lval = std::stoul(what[2]);
            if (0 == metrics_type.compare("total_reads")) {
                total_count = lval;
                lval_count++;
            } else if (0 == metrics_type.compare("sense_frags")) {
                sense_frag = lval;
                lval_count++;
            } else if (0 == metrics_type.compare("antisense_frags")) {
                antisense_frag = lval;
                lval_count++;
            } else if (0 == metrics_type.compare("aligned_frags")) {
                aligned_frag = lval;
                lval_count++;
            } else if (0 == metrics_type.compare("unmapped_reads")) {
                unmapped_count = lval;
                lval_count++;
            } else if (0 == metrics_type.compare("total_frags")) {
                total_frags = lval;
                lval_count++;
            }
        }
    }

    if (lval_count != 6) {
        std::string lstr = "Failed to extract six main metrics from the "
            "basic metric file.";
        throw std::runtime_error(lstr);
    }

    words.close();

}

/*
void metrics_gen::copy_basic_metrics_file() {
    // Copy the content of the metfile to the outfile fist
    // But we should also load the useful entries 
    std::ifstream  src(metfile_str, std::ios::binary);
    std::ofstream  dst(outfile_str,   std::ios::binary);
    dst << src.rdbuf();
    src.close();
    dst.close();

}
*/

double metrics_gen::get_pct(unsigned long numer, unsigned long denom) {
    if (0 == denom) {
        return 0.0;
    } else {
        double lval = (numer * 100.0) / denom;
        return lval;
    } 
}

void metrics_gen::print_feature_type(std::string feature_name, 
        unsigned long count_s, unsigned long count_as, 
        unsigned long aligned_frag, std::ofstream& outfile) {
    unsigned long count_total = count_s + count_as;
    double feature_pct_of_aligned = get_pct(count_total, aligned_frag);
    double feature_sense_pct = get_pct(count_s, count_total);
    outfile << std::fixed << std::setprecision(4);
    outfile << feature_name << "_pct_of_aligned_frags\t" << 
        feature_pct_of_aligned << "\n";
    outfile << feature_name << "_sense_pct\t" << feature_sense_pct << "\n";
}


void metrics_gen::print_final_metrics() {

    // Open the outfile in the append mode first.
    std::ofstream outfile(outfile_str);
    outfile << "metrics_type\tcount\n";
    outfile << std::fixed << std::setprecision(4);
    outfile << "total_reads\t" << total_count << "\n";

    outfile << "unmapped_reads\t" << unmapped_count << "\n";
    double unmapped_count_pct_of_total = get_pct(unmapped_count, total_count);
    outfile << "unmapped_reads_pct_of_total_reads\t" 
        << unmapped_count_pct_of_total << "\n";

    outfile << "total_frags\t" << total_frags << "\n";
    outfile << "aligned_frags\t" << aligned_frag << "\n";
    double aligned_frag_pct_of_total = get_pct(aligned_frag, total_frags);
    outfile << "aligned_frags_pct_of_total_frags\t" 
        << aligned_frag_pct_of_total << "\n";       


    outfile << "sense_frags\t" << sense_frag << "\n";
    double sense_frag_pct_of_aligned = get_pct(sense_frag, aligned_frag);
    outfile << "sense_frags_pct_of_aligned_frags\t" << 
        sense_frag_pct_of_aligned << "\n";    

    print_feature_type("ribosomal_rna_frags", ribosomal_rna_s, ribosomal_rna_as, aligned_frag, outfile);
    print_feature_type("protein_coding_frags", protein_coding_s, protein_coding_as, aligned_frag, outfile);
    print_feature_type("pseudogene_frags", pseudogene_s, pseudogene_as, aligned_frag, outfile);
    print_feature_type("long_noncoding_frags", long_noncoding_s, long_noncoding_as, aligned_frag, outfile);
    print_feature_type("short_noncoding_frags", short_noncoding_s, short_noncoding_as, aligned_frag, outfile);
    print_feature_type("uncategorized_frags", uncategorized_s, uncategorized_as, aligned_frag, outfile);
    if (use_umi) {
        print_feature_type("umi_drop_frags", s_drop_count, as_drop_count, aligned_frag, outfile);
    }

    std::cout << "ribosomal_rna_s: " << ribosomal_rna_s << "\n";
    std::cout << "ribosomal_rna_as: " << ribosomal_rna_as << "\n";
    std::cout << "protein_coding_s: " << protein_coding_s << "\n";
    std::cout << "protein_coding_as: " << protein_coding_as << "\n";
    std::cout << "pseudogene_s: " << pseudogene_s << "\n";
    std::cout << "pseudogene_as: " << pseudogene_as << "\n";
    std::cout << "long_noncoding_s: " << long_noncoding_s << "\n";
    std::cout << "long_noncoding_as: " << long_noncoding_as << "\n";
    std::cout << "short_noncoding_s: " << short_noncoding_s << "\n";
    std::cout << "sort_noncoding_as: " << short_noncoding_as << "\n";
    std::cout << "uncategorized_s: " << uncategorized_s << "\n";
    std::cout << "uncategorized_as: " << uncategorized_as << "\n";
    std::cout << "aligned_frag: " << aligned_frag << "\n";
    if (use_umi) {
        std::cout << "umi_drop_s: " << s_drop_count << "\n";
        std::cout << "umi_drop_as: " << as_drop_count << "\n";
    }
    
    

    outfile.close();

    std::cout << "check_s " << check_s << "\n";
    std::cout << "check_as " << check_as << "\n";
}


int main(int argc, char* argv[]) {
    
    metrics_gen lmgen;
    bool all_set = true;

    try {
        all_set = lmgen.parse_args(argc, argv);   
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;

    } catch (...) {
        return 0;
    }

    if (!all_set) {
        lmgen.print_help();
        return 0;
    }

    try {
        lmgen.init_gene_biotypes();
        lmgen.load_gene_count();
        lmgen.load_basic_metrics();
        if (lmgen.use_umi) {
            lmgen.load_basic_metrics_umi();
        }
        lmgen.print_final_metrics();
    }  catch(std::invalid_argument& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}


