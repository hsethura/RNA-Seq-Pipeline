#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib> 
#include <cstring>
#include <set>
#include <cctype>
#include <memory>
#include <set>

#include "BKTree.h"
#include "fastq_reader.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

struct classcomp {
	bool operator() (const double lhs, const double rhs) const {
		return rhs > lhs;
	}
};


class bc_splitter {

	public:
	//bc_splitter();
	bool parse_args(int argc, char* argv[]);
	int distance(std::string source, std::string target, bool remove_last = false);	
	unsigned long updateMaps(std::string& barcode_str,
    	std::string& lword1, std::string& lword2,
    	std::string& lword3, std::string& lword4,
    	unsigned long totalcap);	
	void writeMapsToFile();
	void split_engine();
	void write_log();
	BKTree<std::string>& getTree();
	void initialize();
	void print_help();
	bool isAlpha(const std::string &str);
	bool isNumber(const std::string& str);
    bool has_suffix(const std::string &str, const std::string &suffix);
    std::unique_ptr<bio::filtering_istream> get_instream(std::string infile_str);

	private:
	int cutoff;
	std::string ltype;
	std::string dict_file;
	std::string file_str;
	std::string file2_str;
	std::string prefix_str;
	std::string outdirpath;
	int barcode_start;
	int barcode_size;
	int umi_start;
	int umi_size;
	int allowed_MB;
	std::string bc_used_file;
	std::string bc_all_file;
	//std::map<std::string, std::vector<std::string>> lQueueMap;
    //std::map<std::string, std::vector<std::string>> rQueueMap;
	std::map<std::string, std::vector<std::unique_ptr<std::string>>> lQueueMap;
    
	std::map<std::string, std::vector<std::unique_ptr<std::string>>> rQueueMap;
	std::set<std::string> barcode_set;
	std::set<std::string> outfile_set;
	std::map<std::string, unsigned long> zero_dist_map;
	std::map<std::string, unsigned long> one_dist_map;
	std::map<std::string, unsigned long> higher_dist_map;
    BKTree<std::string> tree;
	po::options_description desc;
	std::map<int, int> distmap;
	std::multimap<double, std::string, classcomp> bar_map;
	std::map<int, std::string> barseq_map;
	std::set<std::string> used_barcodes;
    std::set<std::string> all_nodes;

	unsigned long totalcap = 0;
	unsigned long  match_total = 0;
	unsigned long ambiguous_total = 0;
    unsigned long no_match_total = 0;

    // Value of keep_last would be set by the command line option
    bool keep_last;

};

class my_exception : public std::exception {
	public:
 	my_exception(const std::string& msg) : msg_(msg) {}
	const char* what(); // override what to return msg_;
	private:
    std::string msg_;
};


void bc_splitter::print_help() {
    std::cout << desc << "\n";
	std::cout << "Usage: bc_splitter_rts -d <dict_file> --file <file>"
			" -p <prefix_str> -o <outdir>\n\n";
}


BKTree<std::string>& bc_splitter::getTree() {
	return tree;
}


void bc_splitter::initialize() {
	std::ifstream iff(dict_file);
    boost::archive::text_iarchive iar(iff);
        
    iar >> tree;
    all_nodes = tree.get_nodes();

	struct stat st = {0};

	if (stat(outdirpath.c_str(), &st) == -1) {
		mkdir(outdirpath.c_str(), 0755);
	}
}

bool bc_splitter::isAlpha(const std::string& str) {
    for(int i = 0; i < str.size(); i++)
        if(!isalpha(str[i]))
            return false;
    return true;
}

bool bc_splitter::isNumber(const std::string& str) {
	for (int i = 0; i < str.size(); i++)
		if (!isdigit(str[i]))
			return false;
	return true;
}



bool bc_splitter::parse_args(int argc, char* argv[]) {
	bool all_set = true;
	desc.add_options()
		("help,h", "produce help message")
		("dict-file,d", po::value<std::string>(&dict_file), "Dictionary file")
		("file,f", po::value<std::string>(&file_str), "Read file")
		("prefix,p", po::value<std::string>(&prefix_str), "Prefix string")
		("outdir,o", po::value<std::string>(&outdirpath), "Output directory")	
        ("keep_last,k", "Optional/Do use last base of barcode (RNATag-Seq)")
		("mismatch,m", po::value(&cutoff)->default_value(1), 
			"Optional/Maximum allowed mismatches.")
		("allowed-mb", po::value(&allowed_MB)->default_value(2048),
			"Optional/Estimated memory requirement in MB.")
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

	barcode_start = 0;
	barcode_size = 9;

	std::cout << "Max mismatch is set to " << cutoff << ".\n";
	std::cout << "Barcode-start is set to " << barcode_start << ".\n";
	std::cout << "Barcode-size is set to " << barcode_size << ".\n";



	if (vm.count("file")) {
		std::cout << "Fastq file is set to: " << file_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Fastq file is not set.\n";
	}


	if (vm.count("prefix")) {
		std::cout << "Prefix string is set to: " << prefix_str << ".\n";
	} else {
		std::cout << "Error: Prefix string is not set.\n";
	}

	if (vm.count("dict-file")) {
		std::cout << "Dict_file is set to " << dict_file << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Dict_file is not set.\n";
	}

	if (vm.count("outdir")) {
		std::cout << "Outdir is set to: " << outdirpath << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Outdir is not set.\n";
	}
    keep_last = vm.count("keep_last");
    printf("keep_last set to: %s\n", keep_last ? "true" : "false");
	
	return all_set;
}

int bc_splitter::distance(std::string source, std::string target, bool remove_last) {

    const int n = source.length();
    const int m = target.length();
    if (n == 0) {
        throw std::invalid_argument("Length of source is zero.");
    }
    if (m == 0) {
        throw std::invalid_argument("Length of target is zero.");
    }

    if (m !=n ) {
        throw std::invalid_argument("Source and target have different length");
	}

    int ldist = 0;

    int n1 = n;
    if (remove_last)
    {
        n1--;
    }

    for (int j = 0; j < n1; j++) {
        if (source[j] != target[j]) {
            ldist++;
        }
    }

    return ldist;
}	


unsigned long 
bc_splitter::updateMaps(std::string& barcode_str, 
	std::string& lword1, std::string& lword2_ori, 
	std::string& lword3, std::string& lword4_ori, 
	unsigned long totalcap) {

	// Get the capacity of the the eight strings as we have to insert
	// them into the two maps.

	// Trim lword2 and lword4. Trim the first 9 bases

	int lword2_len = lword2_ori.length();
	int lword4_len = lword4_ori.length();
	int lword2_lefts = lword2_len - barcode_size;
	int lword4_lefts = lword4_len - barcode_size;

	std::string lword2 = lword2_ori.substr(barcode_size, lword2_lefts);
	std::string lword4 = lword4_ori.substr(barcode_size, lword4_lefts);
	
	int lcap1 = lword1.capacity();
	int lcap2 = lword2.capacity();
	int lcap3 = lword3.capacity();
	int lcap4 = lword4.capacity();

	
	int allcap = lcap1 + lcap2 + lcap3 + lcap4 ;

	totalcap += allcap;

	lQueueMap[barcode_str].push_back(std::make_unique<std::string>(lword1));
	lQueueMap[barcode_str].push_back(std::make_unique<std::string>(lword2));
	lQueueMap[barcode_str].push_back(std::make_unique<std::string>(lword3));
	lQueueMap[barcode_str].push_back(std::make_unique<std::string>(lword4));


	return totalcap;
}

void bc_splitter::writeMapsToFile() {
	
	for (auto& kv : lQueueMap) {
	    std::string barcode = kv.first;

        const std::string file1 = outdirpath + "/" + prefix_str + "_" + barcode + "_R.fastq";

        std::ofstream ofs1;
        std::ofstream ofs2;
        if (outfile_set.count(barcode) > 0) {
	        ofs1 = std::ofstream(file1, std::ofstream::out|std::ofstream::app);
        } else {

            ofs1 = std::ofstream(file1, std::ofstream::out|std::ofstream::trunc);
            outfile_set.insert(barcode);
        }

        // Dump the content of the two maps to the two files
		std::vector<std::unique_ptr<std::string>> & valSet1 = lQueueMap[barcode];

        //std::vector<std::string> & valSet1 = lQueueMap[barcode];
        //std::vector<std::string> & valSet2 = rQueueMap[barcode];

        for (auto const& kv1 : valSet1) {
 	        std::string val1 = *kv1;
            ofs1 << val1 << '\n';
        }

        ofs1.close();

	}

	lQueueMap.clear();

} 

		
void bc_splitter::split_engine() {

	// The words starting with lword is for read 1
	std::string lword1;
	std::string lword2;
	std::string lword3;
	std::string lword4;


	for (int j = 0; j <= cutoff; j++) {
		distmap[j] = 0;
	}   
	const unsigned long MB_SIZE = 1024 * 1024;
	unsigned long total_allowed = allowed_MB * MB_SIZE;

	//const std::string logfile_detailed = outdirpath + "/logfile_detailed.txt";
	//std::ofstream log_detailed(logfile_detailed);

	// We shall start reading the first line. The assumption is that second line 
	// contains the barcode.

    fastq_reader file1(file_str);

    std::cout << "Here we are too!\n";

	while (file1.getline(lword1)) {

		if (!file1.getline(lword2)) {break;}
		if (!file1.getline(lword3)) {break;}
		if (!file1.getline(lword4)) {break;}
		

		// The barcode stays at the second line of each four lines of first
		// read file.
        
		// We want to the 9th bases of the barcode, since it is not useful.
		std::string barcode_str = lword2.substr(barcode_start, barcode_size);	
		std::vector<std::string> results;
		
        bool remove_last = !keep_last;
        //printf("remove_last set to: %s\n", remove_last ? "true" : "false");
		results = tree.find(barcode_str, cutoff, remove_last);

		// calculate the minimum dIstance between the target and references

		int smallest_dist = cutoff + 1;

		// The smallest_barcode initialization could technically be anything
		// since we overwrite this variable.

		std::string smallest_barcode = "JJJJJJJJ";

		std::vector<int> dist_vec;
		for (auto const& val : results) {
			int ldist = distance(val, barcode_str, remove_last);
			if (ldist < smallest_dist) {
				smallest_dist = ldist;
				smallest_barcode = val;
			}

			dist_vec.push_back(ldist);
		}

		int smallest_count = 0;
		for (auto const& temp_dist : dist_vec) {
			if (temp_dist == smallest_dist) {
				smallest_count++;
			}
		}

		//std::cout << "actual_barcode: " << barcode_str << 
		//	", smallest barcode: " <<  smallest_barcode <<  
		//	", sallest dist: " << smallest_dist << 
		//	", smallest_count: " << smallest_count << "\n";

		// So the smallest dist has to be unique, otherwise we shall put 
		// 	them in a fil called unknow.

		std::string write_barcode;

		if (smallest_count == 1) {
			write_barcode = smallest_barcode;
			if (smallest_dist == 0) {
				zero_dist_map[write_barcode]++;
			} else if (smallest_dist == 1) {                        
				one_dist_map[write_barcode]++;
			} else {
				higher_dist_map[write_barcode]++;
			}
			match_total++;
				
		} else if (smallest_count > 0) {
			write_barcode = "ambiguous";
			ambiguous_total++;
		} else {
			write_barcode = "no_match";
			no_match_total++;
		}
		barcode_set.insert(write_barcode);
 
		totalcap = updateMaps(write_barcode, lword1, lword2, lword3, lword4, totalcap);
			
		//std::cout << "total cap: " << totalcap << "\n";

		if (totalcap > total_allowed) {
			writeMapsToFile();
			// Write all the data in the respective files sequentially
			totalcap = 0;
		}

		distmap[smallest_dist]++;
	}

	// final writing to the files
	writeMapsToFile();

	//log_detailed.close();
}

void bc_splitter::write_log() {

	// Writing the logs
	const std::string logfile1 = outdirpath + "/" + prefix_str + "_frequency_logfile.txt";
	std::ofstream log_freq(logfile1);

	std::setprecision(2);

	unsigned long total_reads = match_total + ambiguous_total + no_match_total;

	double ambiguous_percent = ((double) ambiguous_total / (double) total_reads) * 100;
	double no_match_percent = ((double) no_match_total / (double) total_reads) * 100;

	log_freq << "Total reads: " << total_reads << "\n..................\n";

	log_freq << "Ambiguous:\n";
	log_freq << ".................." << "\n";
	log_freq << "Total ambiguous reads: " << ambiguous_total << " (" << ambiguous_percent << "%)\n\n";

	log_freq << "No match:\n";
    log_freq << ".................." << "\n";
    log_freq << "Total non-match reads: " << no_match_total << " (" << no_match_percent << "%)\n\n";


	for (const auto& lbarcode : all_nodes) {

	
		unsigned long total_correct_count = 0;
		double zero_dist_percent = 0;
		double one_dist_percent = 0;
		double higher_dist_percent = 0;
		double barcode_read_percent = 0;

        if (barcode_set.count(lbarcode) > 0) {  

		    unsigned long zero_dist_count = zero_dist_map[lbarcode];
		    unsigned long one_dist_count = one_dist_map[lbarcode];
		    unsigned long higher_dist_count = higher_dist_map[lbarcode];

		    total_correct_count = zero_dist_count + one_dist_count + higher_dist_count;
		    zero_dist_percent = ((double)zero_dist_count / (double)total_correct_count) * 100.0;
		    one_dist_percent = ((double)one_dist_count / (double)total_correct_count) * 100.0;
		    higher_dist_percent = 100 - zero_dist_percent - one_dist_percent;
		    barcode_read_percent = ((double)total_correct_count / (double)total_reads) * 100;  
        }

  		
		log_freq << "Barcode: " << lbarcode << "\n";
		log_freq << ".................." << "\n";
		log_freq << "Zero base mismatch: " << zero_dist_percent << "%\n";
		log_freq << "One base mismatch: " << one_dist_percent << "%\n";
		log_freq << "Total read for this barcode: " << total_correct_count << 
			" (percent of total reads: " << barcode_read_percent << "%)\n";
		log_freq << "\n";
		bar_map.insert(std::pair<double, std::string>(barcode_read_percent, lbarcode));

	}

	log_freq.close();

	std::cout << std::fixed;
    std::cout << std::setprecision(4);	

	for (auto& entry : bar_map) {
		
		std::cout << entry.second << ": " << entry.first << "%\n";
	}

	std::cout << "Ambiguous: " << ambiguous_percent << "%\n";
	std::cout << "No-match: " << no_match_percent << "%\n";
}



int main(int argc, char* argv[]) { 

	bc_splitter lbs;

	// To avoid performance hit but keeping a large number of files open 
	// and writing them in parallel, we have decided to maintain a queue
	// of lines according to barcodes. So, it reality, there would be
	// n number of queues where n is the number of barcodes.
	// Once the total size of the queue exceeds a predefined limit such as
	// 2 GB, we shall write the queu to files. We shall open those files
	// at that time sequentially, update them and close them sequentially. 
	// This is a sort of lazy update, I hope that it would provide good 
	// performence.
	
	bool all_set = true;
	try {
		all_set = lbs.parse_args(argc, argv);	
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << "\n";
        //lbs.print_help();
        return 1;

	} catch (...) {
		//lbs.print_help();
		return 0;
	}

 	if (!all_set) {
		lbs.print_help();
		return 0;
	}

	lbs.initialize();
	try {
		lbs.split_engine();
	} catch(std::invalid_argument& e) {
        std::cerr << "error: " << e.what() << "\n";
		//lbs.print_help();
		return 1;
    }

	lbs.write_log();
        
    return 0;
}



