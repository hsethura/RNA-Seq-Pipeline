#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <cstdlib>

using namespace std;

#define MAX_LINES 10000000000

struct SIPHTInfo {	string name;
	string UpORF_num;
	string UpORF_name;
	string UpORF_dir;
	int UpORF_dist;
	int start;
	int end;
	string dir;
	int DnORF_dist;
	string DnORF_dir;
	string DnORF_num;
	string DnORF_name;
	string can_type;
	int score;
	float expect;
	char term_type;
	int num_BLAST;
	string partners;
	string QRNA;
	string homol;
	int num_para;
	char WRP;
	string paral;
	string synt;
	char UTR;
	string ribo;
	string TFBS;
};
 

void BUFF_OF (unsigned long int, unsigned long int, string);
int FILE_LINE_COUNT (string);
void PRINT_ERROR_FILE (string, int);

//------------- BUFF_OF -----------------------------

void 
BUFF_OF (unsigned long int max, unsigned long int count, string name)
{
	
	if(count== max)
	{
		cout << "BUFFER OVERFLOW IN FUNCTION " << name << ".\n";
		cerr << "BUFFER OVERFLOW IN FUNCTION " << name << ".\n";
		exit(1);
	}
}
//------------- FILE_LINE_COUNT -----------------------------

int 
FILE_LINE_COUNT (string file)
{
	ifstream InFile;
	unsigned long int count=0;
	int line_len=1000000;
	char input[line_len];
	string s;
	
	InFile.open(file.c_str());
	
	if (!InFile)
	{
		cerr << "Error opening " << file << ".\n";
		exit(1);
	}
	
	while(!InFile.eof())
	{
		count++;
		BUFF_OF (MAX_LINES, count, "FILE_LINE_COUNT");  // Report buffer overflow
		InFile.getline (input, line_len);		
		s=input;
		//cerr << input << s << " " << file << endl;
		if (s.empty())			
			count--;
	}
	
	InFile.close();
	return count;
}

//------------- PRINT_ERROR_FILE -----------------------------

void 
PRINT_ERROR_FILE (string message, int number)
{
	ofstream OutFile;
	
	if(number==0)
		OutFile.open("SIPHTseq_err.txt");
	else		
		OutFile.open("SIPHTseq_err.txt", ios::app);

	if (!OutFile)
	{
		cout << "Error opening 'SIPHTseq_err.txt'.\n";
		exit(1);
	}
	
	OutFile << message << endl;
	
	OutFile.close();
}

