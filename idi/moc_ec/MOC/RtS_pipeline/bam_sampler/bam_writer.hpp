#ifndef _BAM_WRITER_HPP
#define _BAM_WRITER_HPP

#include <htslib/sam.h>


class bam_writer {

    public:

    bam_writer(std::string& outfile_str, bam_hdr_t* lhdr1) {
        const char* format = NULL;
        if (has_suffix(outfile_str, "sam")) {
            format = "w";
        } else if (has_suffix(outfile_str, "bam")) {
            format = "wb";
        } else {
            std::string lstr = "File with illegal suffix: " + outfile_str + "\n";
            throw std::runtime_error(lstr);
        }

        const char* outfile_cstr = outfile_str.c_str();
        printf("%s\t%s\n", outfile_cstr, format);
        if (!(fp = sam_open(outfile_cstr, format))) {
            std::cout << "Error in opening sam file" << "\n";
        } else {
            std::cout << "Successfully created the outfile: " << outfile_str << "\n";
        }
        /*
        if(!(lhdr = bam_hdr_init())) {
            std::cout << "Error in bam header init" << "\n";
        } else {
        
            std::cout << "Successfully created bam header" << "\n";
        }
        */
        lhdr = bam_hdr_dup(lhdr1);
        if (sam_hdr_write(fp, lhdr) < 0 ) {
            std::cout << "Error in sam_hdr_write." << "\n";
        } else {
            std::cout << "Wrote the header successfully." << "\n";
        }
        lread_w = bam_init1();



    }

    size_t get_m(size_t var) {
        size_t lvar = (size_t)exp2(ceil(log2(var)));
        return lvar;
    }

    void write_record(const std::string& record) {
        kstring_t str;
        str.s = (char*) record.c_str();
        str.l = record.length();
        str.m = get_m(str.l);
        
        if (sam_parse1(&str, lhdr, lread_w) < 0) {
            std::cout << "Problem with sam_parse1" << "\n";
        } 

        // Now write it into the file
        if (sam_write1(fp, lhdr, lread_w) < 0) {
            std::cout << "Problem with sam_write1" << "\n";
        }  

    }

bool has_suffix(const std::string &str, const std::string &suf)
    {
        return str.size() >= suf.size() &&
           str.compare(str.size() - suf.size(), suf.size(), suf) == 0;
    }

    ~bam_writer() {
        bam_destroy1(lread_w);
        bam_hdr_destroy(lhdr);
        sam_close(fp); // clean up 
    }

    private:
    std::string outfile_str;
    htsFile *fp = NULL;
    bam_hdr_t *lhdr = NULL;
    bam1_t* lread_w = NULL;

      

};
#endif
