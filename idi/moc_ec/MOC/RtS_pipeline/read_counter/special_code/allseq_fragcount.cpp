#include <stdio.h>
#include <iostream>
#include <string>
#include <string.h>
#include <assert.h>
#include <htslib/sam.h>
#include <map>
#include "my_exception.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class allseq_fragc {
    private:
        typedef struct flagvals_s flagvals;
    public:
        po::options_description desc;
        void initialize();
        void print_help();
        bool parse_args(int argc, char* argv[]);
        bool has_suffix(const std::string &str, 
            const std::string &suffix);
        flagvals get_flagvals(int flag);
        std::string get_refname(bam_hdr_t* lhdr, bam1_t* lread);
        void update_single_s(char* lstr);
        void update_single_as(char* lstr);
        void process_single_read(bam_hdr_t* lhdr, bam1_t* lread);
        void main_func();
        void start_allseq();
        void free_vars();
        void print_outfile();
        void print_summary();

    private:
        std::map<std::string, unsigned long> s_fragc;
        std::map<std::string, unsigned long> as_fragc;
        std::map<std::string, unsigned long> s_readc;
        std::map<std::string, unsigned long> as_readc;


        // Specifically for allseq or may be single ended dataset
        unsigned long read_s = 0;
        unsigned long read_as = 0;
        unsigned long check_var = 0;
        unsigned long total_count = 0;
        unsigned long read_unmapped = 0;

        std::string infile_str;
        std::string outfile_str;
        std::string metrics_str;

        htsFile *fp = NULL;
        bam_hdr_t *lhdr = NULL;
        bam1_t* read1 = NULL;
        bam1_t* read2 = NULL;
        bam1_t* lread = NULL;

};

struct flagvals_s{
    int is_paired;
    int is_proper_paired;
    int is_unmapped;
    int is_mate_unmapped;
    int is_reverse;
    int is_mate_reverse;
    int is_read1;
    int is_read2;
};

void allseq_fragc::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: sam_fragcount --proto <protocol> --end <paired info>"
       " --infile <sam/bam> --outfile <output file> --metrics <metrics file>"
       "\n\n";
}

bool allseq_fragc::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output file.")
        ("metrics,m", po::value<std::string>(&metrics_str), "metrics file.")
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

    if (vm.count("metrics")) {
        std::cout << "Metrics file is set to: " << metrics_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: metrics file is not set.\n";
    }

    return all_set;

}

bool allseq_fragc::has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

allseq_fragc::flagvals allseq_fragc::get_flagvals(int flag) {

    flagvals lvals;
    lvals.is_paired = flag & BAM_FPAIRED;
    lvals.is_proper_paired = flag & BAM_FPROPER_PAIR;
    lvals.is_unmapped = flag & BAM_FUNMAP;
    lvals.is_mate_unmapped = flag & BAM_FMUNMAP;
    lvals.is_reverse = flag & BAM_FREVERSE;
    lvals.is_mate_reverse = flag & BAM_FMREVERSE;
    lvals.is_read1 = flag & BAM_FREAD1;
    lvals.is_read2 = flag & BAM_FREAD2;

    return lvals;

}

std::string allseq_fragc::get_refname(bam_hdr_t* lhdr, bam1_t* lread) {
    int32_t core_tid = (lread -> core).tid;
    if (core_tid == -1) {
        throw std::runtime_error("Negative core_tid.");
    } else {
        char* refname = (lhdr -> target_name)[core_tid];
        std::string refname_str(refname);
        return refname_str;
    }
}



void sam_fragc::update_single_s(std::string lstr) {
    s_fragc[lstr]++;
    s_readc[lstr]++;

}

void sam_fragc::update_single_as(std::string lstr) {
    check_var++;
    as_fragc[lstr]++;
    as_readc[lstr]++;
}


// Note that this mapping of forward and reverse reads are based on the 
// standard Illumina protocol and using the second read for allseq.
// In that case, the read on the forward strand would be sense and the
// read on the reverse strand would be anti-sense. If things change in
// the future, we have use a different set of code.

void sam_fragc::process_single_read_allseq(bam_hdr_t* lhdr, bam1_t * lread) {
    uint16_t flag = (lread -> core).flag;
    flagvals lfvals = get_flagvals(flag);
    if (flag == 0) {
        // This is the read on the forward strand; for allseq this is the
        // sense expression.
        std::string ref_str = get_refname(lhdr, lread);
        update_single_s(ref_str);
        read_s++;
    } else if (lfvals.is_reverse) {
        // This is the read on the reverse strand; for allseq this is the 
        // anti-sense expression.
        std::string ref_str = get_refname(lhdr, lread);
        update_single_as(ref_str);
        read_as++;
    } else if (lfvals.is_unmapped){
        // The read is unmapped.
        read_unmapped++; 
    } else {
        // Unexpected situation.
        std::string lstr = "Unusual case for single read.";
        throw std::runtime_error(lstr);
    }
}


void allseq_fragc::print_outfile() {
    std::ofstream outfile(outfile_str);

    for (const auto & e : s_fragc) {
    }
}
void sam_fragc::print_outfile() {
    const char* outfile_name = outfile_str.c_str();
    FILE *fp = fopen(outfile_name, "w");
    fprintf(fp, "geneid\ts_read_count\tas_read_count\ts_frag_count\tas_frag_count\n");
    for (const auto & e : s_fragc) {
        char* gene = e.first;
        unsigned long s_read = s_readc[gene];
        unsigned long as_read = as_readc[gene];
        unsigned long s_frag = s_fragc[gene];
        unsigned long as_frag = as_fragc[gene];
        fprintf(fp, "%s\t%lu\t%lu\t%lu\t%lu\n", gene, s_read, as_read, s_frag, as_frag);
    }    
    fclose(fp);

}


void sam_fragc::start_allseq() {
    lread = bam_init1();
    while(sam_read1(fp, lhdr, lread) >= 0) {
        // For allseq process all reads separately
        process_single_read_allseq(lhdr, lread);    
        total_count++;
    }
    printf("check_var: %lu\n", check_var);
}

void sam_fragc::free_vars() {
    if (0 == end_str.compare("paired")) {
        bam_destroy1(read1);
        bam_destroy1(read2);
    } else if (0 == end_str.compare("single")) {
        bam_destroy1(read1);
    } else {
        std::string lstr = "unknown option during freeing memory: " + end_str;
        throw std::runtime_error(lstr);
    }
    bam_hdr_destroy(lhdr);
    hts_close(fp);

}

void sam_fragc::allseq_print_summary() {
    const char* outfile_name = metrics_str.c_str();
    FILE *fp = fopen(outfile_name, "w");
    fprintf(fp, "sense read:\t%lu\n", read_s);
    fprintf(fp, "anti-sense read:\t%lu\n", read_as);
    fprintf(fp, "unmapped count:\t%lu\n", read_unmapped);
    unsigned long combined_count = read_s + read_as + read_unmapped;
    assert(combined_count == total_count);
    fprintf(fp, "combined count:\t%lu\n", combined_count);
    fprintf(fp, "total count:\t%lu\n", total_count);
    fclose(fp);
}

void sam_fragc::initialize() {
    const char* format = NULL;
    if (has_suffix(infile_str, "sam")) {
        format = "r";
    } else if (has_suffix(infile_str, "bam")) {
        format = "rb";
    } else {
        std::string lstr = "File with illegal suffix: " + infile_str + "\n";
        throw std::runtime_error(lstr);
    }

    const char* infile_cstr = infile_str.c_str(); 
    fp = hts_open(infile_cstr, format);
    lhdr = sam_hdr_read(fp);

}

void sam_fragc::main_func() {
    start_allseq();
}

int main(int argc, char** argv) {

    allseq_fragc sfc;
    bool all_set = true;
    try {
        all_set = sfc.parse_args(argc, argv);   
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;

    } catch (...) {
        return 0;
    }

    if (!all_set) {
        sfc.print_help();
        return 0;
    }
    try {
        sfc.initialize();
        sfc.main_func();
        sfc.print_outfile();
        sfc.allseq_print_summary();
        sfc.free_vars();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;

}


