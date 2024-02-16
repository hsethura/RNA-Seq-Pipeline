#include <iostream>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

namespace po = boost::program_options;

class count_mapper {

    public:
        
        void print_help();
        bool parse_args(int argc, char* argv[]);
        void load_trans_gene_map();
        void build_gene_count_map();
        void build_gene_umi_count_map();
        void print_outfile();
        std::string get_pre(std::string& lstr, char ldel = ':');
        bool use_umi;
    private:
        std::unordered_map<std::string, std::string> trans_gene;
        std::vector<std::string> gene_vec;
        std::set<std::string> gene_set;
        std::unordered_map<std::string, unsigned long> s_gene_count; 
        std::unordered_map<std::string, unsigned long> as_gene_count; 
        std::unordered_map<std::string, unsigned long> s_gene_umi_count; 
        std::unordered_map<std::string, unsigned long> as_gene_umi_count; 

        std::string infile_str;
        std::string outfile_str;
        std::string umi_metrics_str;
        std::string mapfile_str;
        po::options_description desc;

};

bool count_mapper::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "Input transcript count file.")
        ("mapfile,m", po::value<std::string>(&mapfile_str), "Transcript-gene mapfile.")
        ("outfile,o", po::value<std::string>(&outfile_str), "Output gene count file.")
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

    if (vm.count("outfile")) {
        std::cout << "Outfile is set to " << outfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outfile is not set.\n";
    }

    if (use_umi) {
        std::cout << "Logfile is set to " << outfile_str << ".log" << "\n";
    }

    if (vm.count("mapfile")) {
        std::cout << "Mapfile is set to " << mapfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: mapfile is not set.\n";
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

void count_mapper::load_trans_gene_map() {
    std::ifstream words(mapfile_str);
    std::string lstr;
    boost::regex expr ("(\\S+)\\s+(\\S+)");
    boost::smatch what;
    if (!words.is_open()) {
        std::string lstr =  "The mapfile cannot be open!\n";
        throw std::runtime_error(lstr);
    } else {
        while (std::getline(words, lstr)) {
            bool res = boost::regex_search(lstr, what, expr);
            if (res) { 
                std::string transcript = what[1];
                std::string gene = what[2];
                trans_gene[transcript] = gene;
                auto search = gene_set.find(gene);
                if (search == gene_set.end()) {
                   // The gene is not found in the set. So, we have to insert
                   // the gene at the end of gene_vec
                   gene_vec.push_back(gene); 
                   gene_set.insert(gene);
                }
            } else {
                std::string lstr = "parse error: transcript-to-gene map file.\n";
                throw std::runtime_error(lstr);
            }
        }
    }
    std::cout << "Loaded " << trans_gene.size() << " entries" << "\n";

}

void count_mapper::build_gene_count_map() {
    // Note that infile would have a header. Also, we need to use the first 
    // column and the last two columns. A good idea could be read the header
    // and put some assert.
    std::ifstream words(infile_str);
    std::string lstr;
    boost::smatch what;
    std::string header;
    boost::regex expr ("(\\S+)\\s+\\d+\\s+\\d+\\s+(\\d+)\\s+(\\d+)");
    if (!words.is_open()) {
        std::string lstr = "The infile cannot be open!\n";
        throw std::runtime_error(lstr);
    } else {
        // read the first line
        std::getline(words, header);
        //check_header(header);
        while (std::getline(words, lstr)) {
            bool res = boost::regex_search(lstr, what, expr);
            if (res) {
               std::string transcript = what[1];
               std::string s_str = what[2];
               unsigned long s_fragc = std::stoul(s_str);
               std::string as_str = what[3];
               unsigned long as_fragc = std::stoul(as_str);
               std::string gene = trans_gene[transcript];
               s_gene_count[gene] += s_fragc;
               as_gene_count[gene] += as_fragc;
            }
        }
    }
}

std::string count_mapper::get_pre(std::string& lstr, char ldel) {
    for (int j = 0; j < lstr.length(); j++) {
        if (lstr[j] == ldel) {
            return lstr.substr(0, j);
        }
    }
    return lstr;
}

void count_mapper::build_gene_umi_count_map() {

    std::string logfile_str = outfile_str + ".log";
    std::ofstream umi_metrics(umi_metrics_str);
    std::ofstream logfile(logfile_str);
    unsigned long s_drop_count = 0;
    unsigned long as_drop_count = 0;

    // Note that infile would have a header. Also, we need to use the first 
    // column and the last two columns. A good idea could be read the header
    // and put some assert.
    std::ifstream words(infile_str);
    std::string lstr;
    boost::smatch what;
    std::string header;
    boost::regex expr ("(\\S+)\\s+(\\S+)\\s+\\d+\\s+\\d+\\s+(\\d+)\\s+(\\d+)");
    if (!words.is_open()) {
        std::string lstr = "The infile cannot be open!\n";
        throw std::runtime_error(lstr);
    } else {
        // read the first line
        std::getline(words, header);
        //check_header(header);
        while (std::getline(words, lstr)) {
            bool res = boost::regex_search(lstr, what, expr);
            if (res) {
               std::string transcript = what[1];
               std::string umi = what[2];
               std::string s_str = what[3];
               unsigned long s_fragc = std::stoul(s_str);
               std::string as_str = what[4];
               unsigned long as_fragc = std::stoul(as_str);
               // get the gene for that transcript
               std::string gene = trans_gene[transcript];
               // Map transcript:umi => gene:umi
               std::string gene_umi = gene + ":" + umi;
               std::string ltag = "gene_" + get_pre(gene) + ", umi_" + umi + ", transcript_" + get_pre(transcript); 
               std::string ltag_s =  ltag + ", strand_S, s_frag_count: " + std::to_string(s_fragc);

               // This is where we map a transcript level UMI to 
               // a gene level UMI. If a UMI appears for a gene
               // We shall not allow any other UMI read from 
               // any other transcript of that gene.
               if (s_fragc >= 1) {
                   if (s_gene_umi_count[gene_umi] == 0) {
                       s_gene_umi_count[gene_umi] = 1;
                       s_gene_count[gene]++;
                       logfile << ltag_s << ", no conflict, counted";
                       s_drop_count += (s_fragc - 1);
                   } else {
                       logfile << ltag_s << ", conflict, not counted";
                       s_drop_count += s_fragc;
                   }
                   logfile << ", s_gene_count: " << std::to_string(s_gene_count[gene]) << "\n";
               } else {
                   //logfile << ltag_s << ", zero, not counted\n";
               }

 
               std::string ltag_as = ltag + ", strand_AS, as_frag_count: " + std::to_string(as_fragc);
               if (as_fragc >= 1) {
                   if (as_gene_umi_count[gene_umi] == 0) {
                       as_gene_umi_count[gene_umi] = 1;
                       as_gene_count[gene]++;
                       logfile << ltag_as << ", no conflict, counted";
                       as_drop_count += (as_fragc - 1);
                   } else {
                       logfile << ltag_as << ", conflict, not counted";
                       as_drop_count += as_fragc;
                   }
                   logfile << ", as_gene_count: " << std::to_string(as_gene_count[gene]) << "\n";
               } else {
                   //logfile << ltag_as << ", zero, not counted\n";
               }

               logfile << "...............\n";
            }
        }
    }
   
    // Write the log for dropping frags related to umi normalization
    umi_metrics << "metrics_type" << "\t" << "count" << "\n";
    umi_metrics << "s_drop_count" << "\t" << s_drop_count << "\n"; 
    umi_metrics << "as_drop_count" << "\t" << as_drop_count << "\n"; 
}

void count_mapper::print_outfile() {
    std::ofstream outfile(outfile_str);

    outfile << "Geneid" << "\t" << "Fragcount" << "\n";
    for (const auto& gene : gene_vec) {
        unsigned long s_count = s_gene_count[gene];
        outfile << gene << "\t" << s_count << "\n";
    }

    for (const auto& gene: gene_vec) {
        unsigned long as_count = as_gene_count[gene];
        std::string gene_as = "AS:" + gene;
        outfile << gene_as << "\t" << as_count << "\n";
    }

    outfile.close();

}

void count_mapper::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: frag_to_gene_count -i <infile> -m <mapfile> -o <outfile>\n\n";
}


int main(int argc, char* argv[]) {
    
    count_mapper lcm;
    bool all_set = true;

    try {
        all_set = lcm.parse_args(argc, argv);   
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;

    } catch (...) {
        return 0;
    }

    if (!all_set) {
        lcm.print_help();
        return 0;
    }

    try {
        lcm.load_trans_gene_map();

        if (lcm.use_umi) {
            lcm.build_gene_umi_count_map();
        } else {
            lcm.build_gene_count_map();
        }

        lcm.print_outfile();
    }  catch(std::invalid_argument& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
