#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cctype> 
#include <sstream>
#include <cstdlib>

#include "SIPHTseq_info.h"

using namespace std;
/*---------------------------------------------------------------------------
 This program takes in a file with start, stop coords and dir of transripts
 and outputs an array of numbers (1 per line) corresponding to the number of 
 of times that coordinates from  0-largest coord are represented 
 The first column corresponds to the + strand, the second to the - strand
 ----------------------------------------------------------------------------*/


static string usage = "Usage: Coords_to_histogram <path to _coord.txt file> < path to .gff file > <min abundance (default 0)>  <strand specificity> < path to distribution file > [window size (default 1)]\n " ;


struct feature 
{
	string acc; 
	string source;
	string type;
	long int start;
	long int end;
	char strand;
	string tag;
	double counts_sense;
	double counts_antisense;
	float read_overlap_sense;
	float read_overlap_antisense;
	string coord_name;
	char discard;
};

int read_features (string, struct feature[], int);

int main(int argc, char *argv[])
{	
	
	string InGFF=argv[1];
	//string InTAGS=argv[2];

	// Determine number of features
	int numFeatures = FILE_LINE_COUNT(InGFF);
	//int numTags = FILE_LINE_COUNT(InTAGS);

	struct feature f[numFeatures];
	//string tag[numTags];
	int count, n; 
	char input[10000];
	float total_counts=0;

	ifstream InFile1,InFile2;
	ofstream OutFile1;

	InFile1.open(InGFF.c_str());
	if (!InFile1)
	{
		cout << "Error opening " << InGFF << ".\n";
		exit(1);
	} 

// 	InFile2.open(InTAGS.c_str());
// 	if (!InFile2)
// 	{
// 		cout << "Error opening " << InTAGS << ".\n";
// 		exit(1);
// 	} 

// 	OutFile1.open(OutMet.c_str(), ios_base::app);
// 	if (!OutFile1)  
// 	{
// 		cerr << "Error opening file '"<<OutMet<<"'.\n\n";
// 		exit(1);
// 	} 
	
	
	// Read in tags
// 	for (count=0; count < numTags; count++)
// 	{
// 		InFile2 >> tag[count];
// 		cout << tag[count] << endl;
// 	}

	// Read in features
	cerr << "Reading in " << numFeatures << " genomic features..." << endl;
	read_features (InGFF, f, numFeatures);	
	

	for (count=0; count < numFeatures; count++)
	{	
		for (n=count+1; n < numFeatures; n++)
		{	
			if(f[count].acc==f[n].acc && f[count].start==f[n].start && f[count].end==f[n].end && f[count].strand==f[n].strand)
			{
				if(f[count].type=="gene")
					f[count].type=f[n].type;				
				if(f[count].type=="gene")
					f[count].type=f[n].type;
				if(f[n].type=="gene")
					f[n].type=f[count].type;
				if(f[n].type=="exon") 
					f[n].type=f[count].type;
				if((f[n].type=="pseudogene" && f[count].type=="CDS") || (f[count].type=="pseudogene" && f[n].type=="CDS"))
					f[count].type="pseudogene";

				f[count].tag=f[count].tag+";"+f[n].tag;
				f[n].discard='Y';
			}
		}
	}
	for (count=0; count < numFeatures; count++)
	{	
				if(f[count].discard!='Y')
					cout << f[count].acc << "	" << f[count].source << "	" << f[count].type << "	" << f[count].start << "	" << f[count].end << "	.	" << f[count].strand << "	.	" << f[count].tag << endl;
	}

/*	
	// Read in coords, determine overlap with features, and increment f.counts
	cerr << "Reading in " << numCoords << " coordinates, determining overlap with features, and incrementing counts..." << endl;
		
	for(count=0; count < numCoords; count++)
	{
		if(count % 100000 == 0)
			cerr << "Read in "<<count/1e6<< " million total_reads so far"<< endl;
		
		InFile1 >> c.acc >> c.start >> c.end >> c.strand >> c.pair >> c.name;		
		
		if((c.end-c.start) >= min_frag_len && (paired_only=="N" || (paired_only=="Y" && c.pair=='Y')))
		{	
			overlap (c.acc, c.start, c.end, c.strand, f, numFeatures, c.name);
			coords_counted++;
		}
		
		//InFile1.getline (input, 100);
		//cerr << c.start << " " <<  c.name<< endl;
	}
	
	for (count=0; count < numFeatures; count++)
	{	
		cout << f[count].tag  << ":" << f[count].acc << "	" << f[count].counts_sense  << "	" << f[count].coord_name << "	" << "S" << endl;
		cout << f[count].tag  << ":" << f[count].acc<< "	" << f[count].counts_antisense  << "	" << f[count].coord_name << "	" << "AS" << endl;
		
		//cerr << f[count].acc << endl;
		for(n=0; n < numACC; n++)
		{	
			if(ACC[n]==f[count].acc)
				ACC_reads[n]=ACC_reads[n]+total_counts+f[count].counts_sense+f[count].counts_antisense;
				
			
		}
		//total_counts=total_counts+f[count].counts_sense+f[count].counts_antisense;

	}
	
	OutFile1 << endl;
	OutFile1 << "@Total_fragments_mapped_to_all_reps: " << numCoords << endl;
	for(n=0; n < numACC; n++)
		OutFile1 << "@Total_counts_for_rep_"<< ACC[n] << ":\t" << ACC_reads[n] << endl;
*/
	return 0;
}

//###########################
int read_features (string FeatIn, struct feature f[], int numFeatures)
{

	int count, n, i, k;
	string s;
	char input[10000];
	
	ifstream InFile1, InFile2;
	
	InFile1.open(FeatIn.c_str());
	if (!InFile1)
	{	
		cout << "Error opening file " << FeatIn << "\n.";
		exit(1);
	}
	
	//cout << numFeatures << endl;
		
	// Read in genomic features
	for (count=0; count < numFeatures; count++)
	{
		
		f[count].discard='N';
		InFile1 >> f[count].acc;		
		
		if(f[count].acc[0] == '#')
		{	
			InFile1.getline(input, 100000);
		}
		else
		{	
			InFile1 >> f[count].source >> f[count].type >> f[count].start >> f[count].end >> s >> f[count].strand >> s >> f[count].tag;
		}
				
		if(InFile1.eof())
			break;
		
	}
		
	InFile1.close();
	return 0;
}


