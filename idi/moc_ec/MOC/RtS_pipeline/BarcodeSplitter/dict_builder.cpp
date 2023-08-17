#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib> 
#include <cstring>

#include "BKTree.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class dict_builder{

	public:
	bool parse_args(int argc, char* argv[]);
	bool build_data();
	void save_data();
	std::string& get_type();
	void print_help();
	dict_builder();
	~dict_builder();

	private:
	BKTree<std::string> tree;
	std::string infile;
	std::string outfile;
	std::string ltype;
    po::options_description desc;

};

class my_exception : public std::exception {
    public:
    my_exception(const std::string& msg) : msg_(msg) {}
    const char* what(); // override what to return msg_;
    private:
    std::string msg_;
};



void dict_builder::print_help() {
	std::cout << desc << "\n";
    std::cout << "Usage: dict_builder -i <infile> -o <outfile>\n\n";
}

dict_builder::dict_builder() {
	tree = BKTree<std::string>();
}

std::string& dict_builder::get_type() {
	return ltype;
}

dict_builder::~dict_builder() {
	// Nothing yet	
}

bool dict_builder::build_data() {
	std::ifstream words(infile);
	std::string lstr;
	boost::regex expr ("(\\w+)\\s+(\\w+)");
	boost::smatch what;
	if (!words.is_open()) {
    	std::cerr << "The infile cannot be open!\n";
		return false;
	} else {
		
		while (std::getline(words, lstr)) {
			bool res = boost::regex_search(lstr, what, expr);
			if (res) {
				std::string seq = what[1];
				std::string barcode = what[2];
				tree.insert(barcode);
			} else {
				throw my_exception("Problem in the parsing the barcode lines.\n");
			}
		}
	}
	std::cout<< "Loaded " <<tree.size()<< " entries" << std::endl;
	return true;
}

void dict_builder::save_data() {
    
	std::ofstream ofs(outfile);
    boost::archive::text_oarchive oa(ofs);
    oa << tree;
	ofs.close();
	return;
}

bool dict_builder::parse_args(int argc, char* argv[]) {
	 bool all_set = true;
	 bool none_set = true;

    desc.add_options()
        ("help,h", "produce help mesage")
        ("infile,i", po::value<std::string>(&infile), "Input file")
        ("outfile,o", po::value<std::string>(&outfile), "Output file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
		//print_help();
        return 0;
    } else {
		//all_set =false;
	}

    if (!vm.count("infile")) {
        std::cout << "Infile not set.\n";
        all_set = false;
		
    } else {
		
        std::cout << "Infile is set to: " << infile << ".\n";
    }

    if (!vm.count("outfile")) {
        std::cout << "Outfile not set.\n";
        all_set = false;
    } else {
        std::cout << "Outfile is set to: " << outfile << ".\n";
    }


	return all_set;
}

int main(int argc, char* argv[]) { 
	dict_builder ldict = dict_builder();
	bool all_set = true;
	
	po::options_description desc("Allowed options");
	try {

		all_set = ldict.parse_args(argc, argv);
	} catch(std::exception& e) {
       	std::cerr << "error: " << e.what() << "\n";
		//ldict.print_help();
       	return 1;
   	}
   	catch(...) {
       	std::cerr << "Exception of unknown type!\n";
		//ldict.print_help();
   	}

	if (!all_set) {
		ldict.print_help();
		return 0;
	}

	ldict.build_data();

	ldict.save_data();
        
    return 0;
}

