#include <iostream>
#include <string>
#include <tuple>

#include <boost/program_options.hpp>

#include "fastq_reader.hpp"
#include "fastq_writer.hpp"

namespace po = boost::program_options;

class read_trimmer {
    public:

        bool parse_args(int argc, char* argv[]);
        std::tuple<bool, bool, std::string> 
            exe_keep_str_5p(std::string& instr, int lcount);    
        std::tuple<bool, bool, std::string> 
            exe_keep_str_3p(std::string& instr, int lcount);
        std::tuple<bool, bool, std::string> 
            exe_trim_str_5p(std::string& instr, int lcount);
        std::tuple<bool, bool, std::string> 
            exe_trim_str_3p(std::string& instr, int lcount);
        std::tuple<bool, bool, std::string>
            exe_trim_str_5p_3p(std::string& instr,
            unsigned int lcount_5p, unsigned int lcount_3p);
        void exe_trim_5p_3p(std::string& infile_str, 
            std::string& outfile_str);
        void exe_keep_5p(std::string& infile_str, 
            std::string& outfile_str);
        void exe_keep_3p(std::string& infile_str, 
            std::string& outfile_str);
        void initialize();
        void write_log();
        void trim_engine();
        void print_help();
        void set_suffix();
        bool has_suffix(const std::string &str, const std::string &suffix);

    private:
        std::string infile_str;
        std::string outfile_str;
        std::string logfile_str;
        int trim_5p_count;
        int trim_3p_count;
        int keep_5p_count;
        int keep_3p_count;
        bool do_trim_5p_3p = false;
        bool do_keep_5p = false;
        bool do_keep_3p = false;
        po::options_description desc;
};

class my_exception : public std::exception {
    public:
    my_exception(const std::string& msg) : msg_(msg) {}
    const char* what(); // override what to return msg_;
    private:
    std::string msg_;
};


void read_trimmer::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: combine_lanes -i <infile> -o <outfile> "
        " --l <logfile> (optional) "
        " --trim_5p <number of bases to trim on 5p side> "
        " --trim_3p <number of bases to trim on 3p side> "
        " --keep_5p <number of bases to keep from 5p side> "
        " --keep_3p <number of bases to keep from 3p side>\n\n";
        
}

bool read_trimmer::parse_args(int argc, char* argv[]) {

    bool all_set = true;
    desc.add_options()
        ("help,h", "generate help message")
        ("infile,i", po::value<std::string>(&infile_str), "path to the input file")
        ("outfile,o", po::value<std::string>(&outfile_str), "path to the output file")
        ("logfile, o", po::value<std::string>(&logfile_str)->default_value(""), "Path to the logfile")
        ("trim_5p", po::value(&trim_5p_count)->default_value(0), "number of bases to trim from 5p")
        ("trim_3p", po::value(&trim_3p_count)->default_value(0), "number of bases to trim from 3p")
        ("keep_5p", po::value(&keep_5p_count)->default_value(-1), "number of bases to keep from the 5p")
        ("keep_3p", po::value(&keep_3p_count)->default_value(-1), "number of bases to keep from the 3p")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);


    if (vm.count("help")) {
        //print_help();
        return 0;
    }else {
        //all_set = false;
    }

    if (vm.count("infile")) {
        std::cout << "The fastq infile is set to : " << infile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: The fastq infile was not set.\n";
    }

    if (vm.count("outfile")) {
        std::cout << "The output file is set to : " << outfile_str << "\n";
    } else {
        all_set = false; 
        std::cout << "Error: The output file is not set.\n";
    }

    std::cout << "logfile is set to : " << logfile_str << "\n";
    std::cout << "trim_5p is set to : " << trim_5p_count << "\n";
    std::cout << "trim_3p is set to : " << trim_3p_count << "\n";
    std::cout << "keep_5p is set to : " << keep_5p_count << "\n";
    std::cout << "keep_3p is set to : " << keep_3p_count << "\n";

    // For simplicity we shall try all the four options one after another
    
    unsigned int command_count = 0;
    if (keep_5p_count >= 0) {
        if (trim_5p_count > 0 || trim_3p_count > 0 || keep_3p_count >= 0) {
            all_set = false;
            std::cout << "Error: No other trimming options with trim_5p_count\n"; 
        } else {
            do_keep_5p = true;
            command_count++;
        }
    }

    if (keep_3p_count >= 0) {
       if (trim_5p_count > 0 || trim_3p_count > 0 || keep_5p_count >= 0) {
            all_set = false;
            std::cout << "Error: No other trimming options with trim_3p_count\n";
        } else {
            do_keep_3p = true;
            command_count++;
        }
    }

    if (trim_5p_count > 0 || trim_3p_count > 0) {
        if (keep_5p_count >= 0 || keep_3p_count >= 0) {
            all_set = false;
            std::cout << "Error: keep_5p or keep_3p not be used with tim_5p\
                 or trim_3p\n";
        } else {
            do_trim_5p_3p = true;
            command_count++;
        }
    }

    if (command_count == 0) {
        throw std::logic_error("No command was used asked!");
    } else if (command_count > 1) {
        throw std::logic_error("More than one commands were requested!");
    }

    return all_set;

}


std::tuple<bool, bool, std::string> 
read_trimmer::exe_keep_str_5p(std::string& instr, int lcount) {
    bool do_write = false;
    bool done_trim = false;
    std::string final_str;
    int instr_len = instr.length();

    if (lcount == -1) {
        std::string err_str = "Error: from exe_keep_str_5p, lcount equals to -1";
        throw std::invalid_argument(err_str);
    } else if (0 == instr.compare("")) {
        do_write = false;
        done_trim = false;
        final_str = "";
    } else {
        // If the length of the read is less than the one to keep then do not
        // write the read in the final file.
        if (lcount == 0) {
            do_write = false;
            done_trim = true;
            final_str = "";
        }
        else if (lcount > instr_len) {
            // Technically a trimming was not possible due to smaller size of 
            // the read.
            do_write = false;
            done_trim = false;
            final_str = "";
        } else {
            // 1 <= lcount <= instr_len
            do_write = true;
            done_trim = true;
            final_str = instr.substr(0, lcount);
        }
    }
    return std::make_tuple(do_write, done_trim, final_str);
}

std::tuple<bool, bool, std::string> 
read_trimmer::exe_keep_str_3p(std::string& instr, int lcount) {
    bool do_write = false;
    bool done_trim = false;
    std::string final_str;
    int instr_len = instr.length();

    if (lcount == -1) {
        std::string err_str = "Error: from exe_keep_str_3p, lcount equals to -1";
        throw std::invalid_argument(err_str);
    } else if (0 == instr.compare("")) {
        do_write = false;
        done_trim = false;
        final_str = "";
    } else {
        if (lcount == 0) {
            do_write = false;
            done_trim = true;
            final_str = "";
        } else if (lcount > instr_len) {
            do_write = false;
            done_trim = false;
            final_str = "";
        } else {
            // 1 <= lcount <= instr_len
            do_write = true;
            done_trim = true;
            int start_pos = instr_len - lcount;
            final_str = instr.substr(start_pos);
        }
    }
    return std::make_tuple(do_write, done_trim, final_str);
}

std::tuple<bool, bool, std::string> 
read_trimmer::exe_trim_str_5p(std::string& instr, int lcount) {
    bool do_write = false;
    bool done_trim = false;
    std::string final_str;
    int instr_len = instr.length();

    if (lcount == -1) {
        // This is basically illegal value, should not be passed through lcount
        std::string err_str = "Error: from exe_trim_str_5p, lcount equals to -1";
        throw std::invalid_argument(err_str);
    } else if (0 == instr.compare("")) {
        do_write = false;
        done_trim = false;
        final_str = "";
    } else {
        if (lcount == 0) {
            do_write = true;
            done_trim = false;
            final_str = instr;
        } else if (lcount > instr_len) {
            do_write = false;
            done_trim = true;
            final_str = "";
        } else {
            // 1 <= lcount <= instr_len
            do_write = true;
            done_trim = true;
            final_str = instr.substr(lcount);
        }
    }
    return std::make_tuple(do_write, done_trim, final_str);
}

std::tuple<bool, bool, std::string> 
read_trimmer::exe_trim_str_3p(std::string& instr, int lcount) {
    bool do_write = false;
    bool done_trim = false;
    std::string final_str;
    int instr_len = instr.length();

    if (lcount == -1) {
        // This is basically illegal value, should not be passed through lcount
        std::string err_str = "Error: from exe_trim_str_3p, lcount equals to -1";
        throw std::invalid_argument(err_str);
    } else if (0 == instr.compare("")) {
        do_write = false;
        done_trim = false;
        final_str = "";
    } else {
        if (lcount == 0) {
            do_write = true;
            done_trim = false;
            final_str = instr;
        } else if (lcount > instr_len) {
            do_write = false;
            done_trim = true;
            final_str = "";
        } else {
            // 1 <= lcount <= instr_len
            do_write = true;
            done_trim = true;
            int lleft = instr_len - lcount;
            final_str = instr.substr(0, lleft);
        }
    }
    return std::make_tuple(do_write, done_trim, final_str);
}

std::tuple<bool, bool, std::string> 
read_trimmer::exe_trim_str_5p_3p(std::string& instr, 
        unsigned int lcount_5p, unsigned int lcount_3p) {
    bool do_write = false;
    bool done_trim = false;
    std::string final_str;

    const auto& tuple_after_5p = exe_trim_str_5p(instr, lcount_5p);
    bool do_write_5p = std::get<0>(tuple_after_5p);
    if (do_write_5p) {
        std::string instr_after_5p = std::get<2>(tuple_after_5p);
        const auto& tuple_after_3p = exe_trim_str_3p(instr_after_5p, lcount_3p);
        bool do_write_3p = std::get<0>(tuple_after_3p);
        if (do_write_3p) {
            do_write = true;
            final_str = std::get<2>(tuple_after_3p);
        }
    }
    return std::make_tuple(do_write, done_trim, final_str);
}

void read_trimmer::exe_trim_5p_3p(std::string& infile_str, 
        std::string& outfile_str) {
    
    fastq_reader infile(infile_str);
    fastq_writer outfile(outfile_str);

    std::string inword1;
    std::string inword2;
    std::string inword3;
    std::string inword4;

    while (infile.getline(inword1)) {

        if (!infile.getline(inword2)) {break;}
        if (!infile.getline(inword3)) {break;}
        if (!infile.getline(inword4)) {break;}
        const auto& inword2_tuple = exe_trim_str_5p_3p(inword2, trim_5p_count, trim_3p_count);
        bool do_write_2 = std::get<0>(inword2_tuple);
        if (do_write_2) {
            
            const auto& inword4_tuple = exe_trim_str_5p_3p(inword4, trim_5p_count, trim_3p_count);
            std::string inword2_tr = std::get<2>(inword2_tuple);
            std::string inword4_tr = std::get<2>(inword4_tuple);
            
            bool do_write_4 = std::get<0>(inword4_tuple);
            if (!do_write_4) {
                std::string err_str = "Error: exe_trim_5p_3p, do_write_4 is\
                    false, while do_write_2 is true.";
                throw std::logic_error(err_str);
            }
            outfile.putline(inword1);
            outfile.putline(inword2_tr);
            outfile.putline(inword3);
            outfile.putline(inword4_tr);
        }
    }
}

void read_trimmer::exe_keep_5p(std::string& infile_str,
        std::string& outfile_str) {

    fastq_reader infile(infile_str);
    fastq_writer outfile(outfile_str);

    std::string inword1;
    std::string inword2;
    std::string inword3;
    std::string inword4;
    while (infile.getline(inword1)) {

        if (!infile.getline(inword2)) {break;}
        if (!infile.getline(inword3)) {break;}
        if (!infile.getline(inword4)) {break;}

        const auto& inword2_tuple = exe_keep_str_5p(inword2, keep_5p_count);
        bool do_write = std::get<0>(inword2_tuple);
        if (do_write) {
            const auto& inword4_tuple = exe_keep_str_5p(inword4, keep_5p_count);
            std::string inword2_tr = std::get<2>(inword2_tuple);
            std::string inword4_tr = std::get<2>(inword4_tuple);
            outfile.putline(inword1);
            outfile.putline(inword2_tr);
            outfile.putline(inword3);
            outfile.putline(inword4_tr);
        }
    }
}

void read_trimmer::exe_keep_3p(std::string& infile_str,
        std::string& outfile_str) {

    fastq_reader infile(infile_str);
    fastq_writer outfile(outfile_str);

    std::string inword1;
    std::string inword2;
    std::string inword3;
    std::string inword4;
    while (infile.getline(inword1)) {

        if (!infile.getline(inword2)) {break;}
        if (!infile.getline(inword3)) {break;}
        if (!infile.getline(inword4)) {break;}

        const auto& inword2_tuple = exe_keep_str_3p(inword2, keep_3p_count);
        bool do_write = std::get<0>(inword2_tuple);
        if (do_write) {
            const auto& inword4_tuple = exe_keep_str_3p(inword4, keep_3p_count);
            std::string inword2_tr = std::get<2>(inword2_tuple);
            std::string inword4_tr = std::get<2>(inword4_tuple);
            outfile.putline(inword1);
            outfile.putline(inword2_tr);
            outfile.putline(inword3);
            outfile.putline(inword4_tr);
        }
    }
}

// Obtained from http://stackoverflow.com/questions/20446201/how-to-check-if-string-ends-with-txt
bool read_trimmer::has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// The assumption is that this would be applied AFTER merging, so the suffix
// would be either _R1.fastq or _R2_fastq, however, if there is an unusual
//suffix it can be applied through the option --suffix
/*
void read_trimmer::set_suffix() {
    std::string R1_suf = "_R1.fastq";
    std::string R1_suf_gz = "_R1.fastq.gz";
    std::string R2_suf = "_R2.fastq";
    std::string R2_suf_gz = "_R2.fastq.gz";
    std::string empty_s = "";

    if (0 != suffix_str.compare(empty_s)) {
        // We already have the suffix provided by the user
        // So, do nothing
    } else if (has_suffix(infile_str, R1_suf) || 
            has_suffix(infile_str, R1_suf_gz)) {
        suffix_str = R1_suf;
    } else if (has_suffix(infile_str, R2_suf) || 
            has_suffix(infile_str, R2_suf_gz)) {
        suffix_str = R2_suf;
    } else {
        std::string err_str = "Error: suffix not extracted!";
        throw std::logic_error(err_str);
    }
}
*/

void read_trimmer::initialize() {
}

void read_trimmer::write_log() {
}

void read_trimmer::trim_engine() {

    
    if (do_keep_5p) {
        exe_keep_5p(infile_str, outfile_str);    
    } else if (do_keep_3p) {
        exe_keep_3p(infile_str, outfile_str);
    } else if (do_trim_5p_3p) {
        exe_trim_5p_3p(infile_str, outfile_str);
    }
}


int main(int argc, char* argv[]) {

    read_trimmer rto;
   
    bool all_set = true;
    try {
        all_set = rto.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        //lbs.print_help();
        return 1;

    } catch (...) {
        //lbs.print_help();
        return 0;
    }

    if (!all_set) {
        rto.print_help();
        return 0;
    }

    rto.initialize();
    try {
        rto.trim_engine();
    } catch(std::invalid_argument& e) {
        std::cerr << "error: " << e.what() << "\n";
        //lbs.print_help();
        return 1;
    }

    rto.write_log();

    return 0;
}

