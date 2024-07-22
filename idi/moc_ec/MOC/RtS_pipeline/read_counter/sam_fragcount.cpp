#include <stdio.h>
#include <iostream>
#include <string>
#include <string.h>
#include <assert.h>
#include <htslib/sam.h>
#include <map>
#include "my_exception.h"
#include <boost/program_options.hpp>
#include <regex>
#include <regex.h>

namespace po = boost::program_options;

class sam_fragc {
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
        char* get_refname(bam_hdr_t* lhdr, bam1_t* lread);
        void throw_flag_error(int flag);
        void update_single_s(char* lstr);
        void update_single_as(char* lstr);
        void update_single_umi_s(char* lstr, char* umi_str);
        void update_single_umi_as(char* lstr, char* umi_str);
        void update_paired_s(char* lstr1, char* lstr2);
        void update_paired_as(char* lstr1, char* lstr2);
        void process_single_read(bam_hdr_t* lhdr, bam1_t* lread, 
            flagvals lfvals);
        void process_single_read_allseq(bam_hdr_t* lhdr, bam1_t* lread);
        void process_paired_read(bam_hdr_t* lhdr, bam1_t* lread1, 
            bam1_t* lread2, flagvals lfvals1, flagvals lfvals2);
        void process_two_reads(bam_hdr_t *lhdr, bam1_t *lread1, bam1_t *lread2);
        void rts_paired_print_summary();
        void main_func();
        void start_allseq();
        void start_rts_single();
        void start_rts_paired();
        void free_vars();
        void print_outfile();
        void print_outfile_umi();
        void allseq_print_summary();
        char* get_umi_str(bam1_t* lread);
        char* regex_match_cstr(const char* regex_str, char* full_rec, int match_len);
        char** get_geneid_umi(char* geneid_umi);
        bool use_umi;

    private:
        std::map<char*, unsigned long> s_fragc;
        std::map<char*, unsigned long> as_fragc;
        std::map<char*, unsigned long> s_readc;
        std::map<char*, unsigned long> as_readc;

        std::map<char*, unsigned long> s_fragc_umi;
        std::map<char*, unsigned long> as_fragc_umi;
        std::map<char*, unsigned long> s_readc_umi;
        std::map<char*, unsigned long> as_readc_umi;

        std::set<char*> qname_umi_set;

        unsigned long read1_s = 0;
        unsigned long read2_s = 0;
        unsigned long read1_s_mu = 0;
        unsigned long read2_s_mu = 0;
        unsigned long read1_s_i = 0;
        unsigned long read2_s_i = 0;
        unsigned long read1_as = 0;
        unsigned long read2_as = 0;
        unsigned long read1_as_mu = 0;
        unsigned long read2_as_mu = 0;
        unsigned long read1_as_i = 0;
        unsigned long read2_as_i = 0;
        unsigned long read_unmapped = 0;
        // Read pairs on same strands
        unsigned long read_both_forward = 0;
        unsigned long read_both_reverse = 0;
        unsigned long total_count = 0;
        unsigned long unequal_pair = 0;

        // Specifically for allseq or may be single ended dataset
        unsigned long read_s = 0;
        unsigned long read_as = 0;
        unsigned long frag_s = 0;
        unsigned long frag_as = 0;

        std::string infile_str;
        std::string outfile_str;
        std::string metrics_str;
        std::string proto_str;
        std::string end_str;

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


void sam_fragc::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: sam_fragcount --proto <protocol> --end <paired info>"
       " --infile <sam/bam> --outfile <output file> --metrics <metrics file>"
       " --proto <rts/allseq>"
       "\n\n";
}

bool sam_fragc::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output file.")
        ("metrics,m", po::value<std::string>(&metrics_str), "metrics file.")
        ("proto,p", po::value<std::string>(&proto_str), "protocol (rts/allseq)")
        ("end,e", po::value<std::string>(&end_str)->default_value("paired"), 
            "paired end info (single/paired)")
        ("use_umi", po::bool_switch(&use_umi)->default_value(false), "Use UMI information with output")
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

    if (vm.count("proto")) {
        std::cout << "Protocol is set to: " << proto_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: protocol is not set.\n";
    }

    if (vm.count("end")) {
        std::cout << "Paired-end info is set to: " << end_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: paired end infor is not set.\n";
    }

    std::cout << "use_umi: " << std::boolalpha << use_umi << "\n"; 
    
    

    return all_set;

}

bool sam_fragc::has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

sam_fragc::flagvals sam_fragc::get_flagvals(int flag) {

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

char* sam_fragc::get_refname(bam_hdr_t* lhdr, bam1_t* lread) {
    int32_t core_tid = (lread -> core).tid;
    char* refname = NULL;
    if (core_tid == -1) {
        throw std::runtime_error("Negative core_tid.");
    } else {
        refname = (lhdr -> target_name)[core_tid];
    }
    if (NULL == refname) {
        throw std::runtime_error("NULL refname.");
    }
    return refname;
}



void sam_fragc::throw_flag_error(int flag) {
    char* flag_str = bam_flag2str(flag);
    printf("Unusal flag_str: %s\n", flag_str);
    free(flag_str);
}

void sam_fragc::update_single_s(char* lstr) {
    s_fragc[lstr]++;
    s_readc[lstr]++;
    frag_s++;

}

void sam_fragc::update_single_as(char* lstr) {
    as_fragc[lstr]++;
    as_readc[lstr]++;
    frag_as++;
}

void sam_fragc::update_single_umi_s(char* lstr, char* umi_str) {
    char* lstr_umi = new char[2 + strlen(lstr)+ strlen(umi_str)];
    sprintf(lstr_umi, "%s:%s", lstr, umi_str);
    s_fragc_umi[lstr_umi]++;
    s_readc_umi[lstr_umi]++;

    if (qname_umi_set.count(lstr_umi) == 0) {
        qname_umi_set.insert(lstr_umi);
    } else {
        delete[] lstr_umi;
    }

    frag_s++;
        
    
}

void sam_fragc::update_single_umi_as(char* lstr, char* umi_str) {
    char* lstr_umi = new char[2 + strlen(lstr)+ strlen(umi_str)];
    sprintf(lstr_umi, "%s:%s", lstr, umi_str);
    as_fragc_umi[lstr_umi]++;
    as_readc_umi[lstr_umi]++;
    if (qname_umi_set.count(lstr_umi) == 0) {
        qname_umi_set.insert(lstr_umi);
    } else {
        delete[] lstr_umi;
    }
    frag_as++;
}

void sam_fragc::update_paired_s(char* lstr1, char* lstr2) {
    assert(0 == strcmp(lstr1, lstr2));
    s_fragc[lstr1]++;
    s_readc[lstr1] += 2;
    frag_s++;
}

void sam_fragc::update_paired_as(char* lstr1, char* lstr2) {
    
    assert(0 == strcmp(lstr1, lstr2));
    as_fragc[lstr1]++;
    as_readc[lstr1] += 2;
    frag_as++;
}


char* sam_fragc::regex_match_cstr(const char* regex_str, char* full_rec,
    int match_len) {
    regex_t regex;
    int nmatch = 2;
    regmatch_t pmatch[2];

    int reti = regcomp(&regex, regex_str, REG_EXTENDED);
    if (reti) {
        fprintf(stderr, "Could not compile regex\n");
        return NULL;
    }


    reti = regexec(&regex, full_rec, 2, pmatch, 0);
    if (!reti) {
        char* result = new char[match_len + 1];
        int len = pmatch[1].rm_eo - pmatch[1].rm_so;
        memcpy(result, full_rec + pmatch[1].rm_so, len);
        result[len] = 0;
        regfree(&regex);
        return result;
    }
    return NULL;
}


char* sam_fragc::get_umi_str(bam1_t* lread) {

    // Get the query_name
    char* qname_src = bam_get_qname(lread);
    int qname_len = strlen(qname_src) + 1;
    char* qname = new char[qname_len];
    memcpy(qname, qname_src, qname_len);

    // Get the UMI string by regex

    const char* regex_str = "^\\S+?umi_(\\w{6}).*$";
    int match_len = 6;
    char* umi_str = regex_match_cstr(regex_str, qname, match_len);

    if (qname != nullptr) {
        delete[] qname;
    }

    if (!umi_str) {
        std::string err_str = "umi str not found, qname: " + std::string(qname);
        throw std::runtime_error(err_str);
    }


    return (umi_str);

}


// Note that this mapping of forward and reverse reads are based on the 
// standard Illumina protocol and using the second read for allseq.
// In that case, the read on the forward strand would be sense and the
// read on the reverse strand would be anti-sense. If things change in
// the future, we have use a different set of code.

void sam_fragc::process_single_read_allseq(bam_hdr_t* lhdr, bam1_t * lread) {
    uint16_t flag = (lread -> core).flag;
    flagvals lfvals = get_flagvals(flag);
    // Get umi str
    if (flag == 0) {
        // This is the read on the forward strand; for allseq this is the
        // sense expression.
        char* ref_str = get_refname(lhdr, lread);
        if (use_umi) {
            char* umi_str = get_umi_str(lread);
            update_single_umi_as(ref_str, umi_str);
            if (umi_str != nullptr) {
                delete[] umi_str;
            }
        } else {
            update_single_as(ref_str);
        }
        read_as++;
    } else if (lfvals.is_reverse) {
        // This is the read on the reverse strand; for allseq this is the 
        // anti-sense expression.
        char* ref_str = get_refname(lhdr, lread);
        if (use_umi) {
            char* umi_str = get_umi_str(lread);
            update_single_umi_s(ref_str, umi_str);
            if (umi_str != nullptr) {
                delete[] umi_str;
            }
        } else {
            update_single_s(ref_str);
        }
        read_s++;
    } else if (lfvals.is_unmapped){
        // The read is unmapped.
        read_unmapped++; 
      
    } else {
        // Unexpected situation.
        std::string lstr = "Unusual case for single read.";
        throw std::runtime_error(lstr);
    }
}


// These are the cases when the mate is unmapped

void sam_fragc::process_single_read(bam_hdr_t* lhdr, bam1_t * lread, flagvals lfvals) {

    assert(lfvals.is_mate_unmapped);
    char* ref_str = get_refname(lhdr, lread);

    if (lfvals.is_read1) {
        if (lfvals.is_reverse) {
            // read1 (sense)
            update_single_s(ref_str);
            read1_s_mu++;
        } else {
            // read 1 (anti sense)
            update_single_as(ref_str);
            read1_as_mu++;
        }
    } else if (lfvals.is_read2) {
        if (lfvals.is_reverse) {
            // read 2 (anti sense)
            update_single_as(ref_str);
            read2_as_mu++;
        } else {
            // read 2 (sense)
            update_single_s(ref_str);
            read2_s_mu++;
        }
    } else {
        std::string lstr = "Unusual case for single read.";
        throw std::runtime_error(lstr);
    }
}

void sam_fragc::process_paired_read(bam_hdr_t* lhdr, bam1_t* lread1, bam1_t* lread2,
        flagvals lfvals1, flagvals lfvals2) {

    // Put asserts so that mate_unmapped is not set for any
    // of the reads
    assert(!lfvals1.is_mate_unmapped);
    assert(!lfvals2.is_mate_unmapped);

    char* ref_str1 = get_refname(lhdr, lread1);
    char* ref_str2 = get_refname(lhdr, lread2);

    if (lfvals1.is_read1 && !lfvals1.is_reverse && !lfvals1.is_mate_reverse &&
        lfvals2.is_read2 && !lfvals2.is_reverse && !lfvals2.is_mate_reverse) {

        // Both on forward strand
        read_both_forward += 2;

    } else if (lfvals1.is_read1 && lfvals1.is_reverse && lfvals1.is_mate_reverse &&
        lfvals2.is_read2 && lfvals2.is_reverse && lfvals2.is_mate_reverse) {
        // Both on forward strand
        read_both_reverse += 2;

    }  else if (lfvals1.is_read1 && lfvals1.is_reverse && 
        lfvals2.is_read2 && lfvals2.is_mate_reverse) {

        // Sense pair
        // Note that we are not differentiating between proper and improper
        // reads, but we are keepin a note on them.

        if (0 == strcmp(ref_str1, ref_str2)) {
            update_paired_s(ref_str1, ref_str2);
        } else {
            update_single_s(ref_str1);
            update_single_s(ref_str2);
            unequal_pair++;
        }
        
        if (lfvals1.is_proper_paired) {
            read1_s++;
        } else {
            read1_s_i++;
        }

        if (lfvals2.is_proper_paired) {
            read2_s++;
        } else {
            read2_s_i++;
        }
    } else if (lfvals1.is_read1 && lfvals1.is_mate_reverse &&
        lfvals2.is_read2 && lfvals2.is_reverse) {

        // Anti sense pair
         if (0 == strcmp(ref_str1, ref_str2)) {
            update_paired_as(ref_str1, ref_str2);
        } else {
            update_single_as(ref_str1);
            update_single_as(ref_str2);
            unequal_pair++;
        }

        if (lfvals1.is_proper_paired) {
            read1_as++;
        } else {
            read1_as_i++;
        }

        if (lfvals2.is_proper_paired) {
            read2_as++;
        } else {
            read2_as_i++;
        }
    } else {
        throw std::runtime_error("Unlikely condition!");
    }
}

void sam_fragc::process_two_reads(bam_hdr_t *lhdr, bam1_t *lread1, bam1_t *lread2) {
    uint16_t flag1 = (lread1 -> core).flag;
    uint16_t flag2 = (lread2 -> core).flag;

    flagvals lfvals1 = get_flagvals(flag1);
    flagvals lfvals2 = get_flagvals(flag2);

    if (lfvals1.is_unmapped && lfvals2.is_unmapped) {
        // Both of them are unmapped; so not tracked. However we can keep 
        // two variables to track it.
        read_unmapped += 2;
    } else if (!lfvals1.is_unmapped && !lfvals2.is_unmapped) {
        // Both of them are mapped
        // identify the first and second
        if (lfvals1.is_read1 && lfvals2.is_read2) {
            process_paired_read(lhdr, lread1, lread2, lfvals1, lfvals2);  
        } else {
            // illegal situation
            std::string lstr = "Both of the reads are labeled as read1 or " 
                " both of the reads are labeld as read 2";
            throw std::runtime_error(lstr);
        }

    } else if (lfvals2.is_unmapped) {
        // fvals1 is mapped
        if (lfvals1.is_read1) {
            process_single_read(lhdr, lread1, lfvals1);
            read_unmapped++;
        } else {
            std::string lstr = "First read is not read 1.";
            throw std::runtime_error(lstr);
        }
        
    } else if (lfvals1.is_unmapped) {
        // fvals2 is mapped
        if (lfvals2.is_read2) {
            process_single_read(lhdr, lread2, lfvals2);
            read_unmapped++;
        } else {
            std::string lstr = "Second read is not read 2.";
            throw std::runtime_error(lstr);
        }

    }
}



void sam_fragc::print_outfile() {
    const char* outfile_name = outfile_str.c_str();
    FILE *fp = fopen(outfile_name, "w");
    fprintf(fp, "transcript\ts_read_count\tas_read_count\ts_frag_count\tas_frag_count\n");
    for (int32_t pos = 0; pos < (lhdr->n_targets); pos++) {
        char* gene = (lhdr -> target_name)[pos];
        unsigned long s_read = s_readc[gene];
        unsigned long as_read = as_readc[gene];
        unsigned long s_frag = s_fragc[gene];
        unsigned long as_frag = as_fragc[gene];
        fprintf(fp, "%s\t%lu\t%lu\t%lu\t%lu\n", gene, s_read, as_read, s_frag, as_frag);
    }    
    fclose(fp);

}

char** sam_fragc::get_geneid_umi(char* geneid_umi) {
    char ldel = ':';
    int lpos =-1;
    unsigned int llen = strlen(geneid_umi);
    for (lpos = llen -1; lpos >= 0; lpos--) {
        if (geneid_umi[lpos] == ldel) {
            break;
        } 
    }
    // Now lpos is supposed to contain the position of the last ":"
    char* geneid = new char[lpos + 1];
    unsigned int umi_len = llen - lpos -1;
    char* umi = new char[ umi_len + 1];
    memcpy(geneid, geneid_umi, lpos);
    geneid[lpos] = 0;
    memcpy(umi, &geneid_umi[lpos + 1], umi_len);
    umi[umi_len] = 0;
    
    char** retarr = new char* [2];
    retarr[0] = geneid;
    retarr[1] = umi;
    return (retarr);
}


void sam_fragc::print_outfile_umi() {
    const char* outfile_name = outfile_str.c_str();
    FILE *fp = fopen(outfile_name, "w");
    fprintf(fp, "transcript\tumi\ts_read_count\tas_read_count\ts_frag_count\tas_frag_count\n");

    for( auto it = qname_umi_set.begin(); it != qname_umi_set.end(); ++it ) {
        char* gene_umi = *it;
        unsigned long s_read = s_readc_umi[gene_umi];
        unsigned long as_read = as_readc_umi[gene_umi];
        unsigned long s_frag = s_fragc_umi[gene_umi];
        unsigned long as_frag = as_fragc_umi[gene_umi];
        char** res = get_geneid_umi(gene_umi);
        char* geneid = res[0];
        char* umi = res[1];
        fprintf(fp, "%s\t%s\t%lu\t%lu\t%lu\t%lu\n", geneid, umi, s_read, as_read, s_frag, as_frag);
        delete[] geneid;
        delete[] umi;
        delete[] res;
        
    }    
        
   fclose(fp);

}


void sam_fragc::start_rts_paired() {
    read1 = bam_init1();
    read2 = bam_init1();
    while(sam_read1(fp, lhdr, read1) >= 0) {
        if (sam_read1(fp, lhdr, read2) < 0) {
            std::string lstr = "Failed to read second read.\n";
            throw std::runtime_error(lstr);
        }
        process_two_reads(lhdr, read1, read2);
        total_count += 2;
    }

}

void sam_fragc::start_rts_single() {
    /*
    lread = bam_init1();
    uint16_t flag = (lread -> core).flag;
    flagvals lfvals = get_flagvals(flag);
    while(sam_read1(fp, lhdr, lread) >= 0) {
        // For allseq process all reads separately
        process_single_read_rts(lhdr, lread, lfvals, flag);    
        total_count++;
    }
    */
}

void sam_fragc::start_allseq() {
    lread = bam_init1();
    while(sam_read1(fp, lhdr, lread) >= 0) {
        // For allseq process all reads separately
        process_single_read_allseq(lhdr, lread);    
        total_count++;
    }
}

void sam_fragc::free_vars() {
    if (0 == proto_str.compare("rts")) {
        if (0 == end_str.compare("paired")) {
            bam_destroy1(read1);
            bam_destroy1(read2);
        } else if (0 == end_str.compare("single")) {
            bam_destroy1(read1);
        } else {
            std::string lstr = "unknown option during freeing memory: " + end_str;
            throw std::runtime_error(lstr);
        }
    } else if (0 == proto_str.compare("allseq")) {
        bam_destroy1(lread);
    } else {
        std::string lstr = "Invalid option: " + proto_str;
        throw std::runtime_error(lstr);
    }

    bam_hdr_destroy(lhdr);
    hts_close(fp);

    
    for( auto it = qname_umi_set.begin(); it != qname_umi_set.end(); ++it ) {
        char* qname_umi = *it;
        if (qname_umi != nullptr) {
            delete[] qname_umi;
        }
    }

}

void sam_fragc::allseq_print_summary() {
    const char* outfile_name = metrics_str.c_str();
    FILE *fp = fopen(outfile_name, "w");

    fprintf(fp, "metrics_type\tcount\n");
    fprintf(fp, "total_reads\t%lu\n", total_count);
    fprintf(fp, "unmapped_reads\t%lu\n", read_unmapped);
    unsigned long aligned_read = read_s + read_as;
    //fprintf(fp, "sense_read:\t%lu\n", read_s);
    //fprintf(fp, "anti_sense_read:\t%lu\n", read_as);
    //fprintf(fp, "aligned_read:\t%lu\n", aligned_read);
    unsigned long total_frags = total_count;
    fprintf(fp, "total_frags\t%lu\n", total_frags);
    
    unsigned long aligned_frag = frag_s + frag_as;
    fprintf(fp, "sense_frags\t%lu\n", frag_s);
    fprintf(fp, "antisense_frags\t%lu\n", frag_as);
    fprintf(fp, "aligned_frags\t%lu\n", aligned_frag);

    unsigned long combined_count = read_s + read_as + read_unmapped;
    if (combined_count != total_count) {
        std::cout << "combined_count: " << combined_count << "\n";
        std::cout << "total_count: " << total_count << "\n";
    }
    assert(combined_count == total_count);
    //fprintf(fp, "combined_count:\t%lu\n", combined_count);
    fclose(fp);
}

void sam_fragc::rts_paired_print_summary() {
    const char* outfile_name = metrics_str.c_str();
    FILE *fp = fopen(outfile_name, "w");
    fprintf(fp, "metrics_type\tcount\n");
    fprintf(fp, "sense_read_1\t%lu\n", read1_s);
    fprintf(fp, "sense_read_2\t%lu\n", read2_s);
    fprintf(fp, "sense_read_1_improper\t%lu\n", read1_s_i);
    fprintf(fp, "sense_read_2_improper\t%lu\n", read2_s_i);
    fprintf(fp, "sense_read_1_mate_unmapped\t%lu\n", read1_s_mu);
    fprintf(fp, "sense_read_2_mate_unmapped\t%lu\n", read2_s_mu);
    fprintf(fp, "antisense_read_1\t%lu\n", read1_as);
    fprintf(fp, "antisense_read_2\t%lu\n", read2_as);
    fprintf(fp, "antisense_read_1_improper\t%lu\n", read1_as_i);
    fprintf(fp, "antisense_read_2_improper\t%lu\n", read2_as_i);
    fprintf(fp, "antisense_read_1_mate_unmapped\t%lu\n", read1_as_mu);
    fprintf(fp, "antisense_read_2_mate_unmapped\t%lu\n", read2_as_mu);
    fprintf(fp, "reads_both_forward_strand\t%lu\n", read_both_forward);
    fprintf(fp, "reads_both_reverse_strand\t%lu\n", read_both_reverse);

    fprintf(fp, "total_reads\t%lu\n", total_count);
    fprintf(fp, "unmapped_reads\t%lu\n", read_unmapped);
    unsigned long total_frags = total_count /2;
    fprintf(fp, "total_frags\t%lu\n", total_frags);
    fprintf(fp, "sense_frags\t%lu\n", frag_s);
    fprintf(fp, "antisense_frags\t%lu\n", frag_as);
    unsigned long aligned_frag = frag_s + frag_as;
    fprintf(fp, "aligned_frags\t%lu\n", aligned_frag);

    unsigned long combined_count = read1_s + read2_s + read1_s_i + read2_s_i + 
        read1_s_mu + read2_s_mu + read1_as + read2_as + read1_as_i + 
        read2_as_i + read1_as_mu + read2_as_mu + read_both_forward + 
        read_both_reverse + read_unmapped;
    
    //printf("combined count:\t%lu\n", combined_count);
    //printf("unequal pair:\t%lu\n", unequal_pair);
    assert(combined_count == total_count);
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
    std::string rts_str("rts");
    std::string single_str("single_str");
    std::string paired_str("paired");
    std::string allseq_str("allseq_str");

    if (0 == proto_str.compare("rts")) {
        if (0 == end_str.compare("single")) {
            start_rts_single();
        } else if (0 == end_str.compare("paired")){
            start_rts_paired();
            rts_paired_print_summary();
        } else {
            std::string lstr = "Invalid rts end info: " + end_str;
            throw std::runtime_error(lstr);
        }
    } else if (0 == proto_str.compare("allseq")) {
        start_allseq();
        allseq_print_summary();
    } else {
        std::string lstr = "Invalid option: " + proto_str;
        throw std::runtime_error(lstr);
    }
}

int main(int argc, char** argv) {

    sam_fragc sfc;
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
        if (sfc.use_umi) {
            sfc.print_outfile_umi();
        } else {
            sfc.print_outfile();
        }
        sfc.free_vars();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;

}


