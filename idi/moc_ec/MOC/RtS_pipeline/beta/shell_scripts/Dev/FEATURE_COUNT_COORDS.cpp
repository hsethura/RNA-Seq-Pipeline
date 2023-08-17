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

struct COORDS 
{
	string acc;
	float start;
	float end;
	char strand;
	char pair;
	string name;
	float q;
};

struct feature 
{
	string acc; 
	string type;
	float start;
	float end;
	char strand;
	string tag;
	double counts_sense;
	double counts_antisense;
	float read_overlap_sense;
	float read_overlap_antisense;
	string coord_name;
	
};

int read_features (string, struct feature[], int);
int overlap (string, float, float, char, feature[], int, string);

int main(int argc, char *argv[])
{	
	
	string InCoords=argv[1];
	string InFeat=argv[2];
	string OutMet=argv[3];
	string ACC_file=argv[4];			
	string paired_only="N";
	
	if(argc == 6)
	{
		string paired_only=argv[5];
	}
	
	// Determine number of features
	int numFeatures = FILE_LINE_COUNT(InFeat);
	// Determine number of coordinates
	unsigned long int numCoords = FILE_LINE_COUNT(InCoords);
	unsigned long int coords_counted=0;
	
	struct feature f[numFeatures];
	struct COORDS c;
	int count, n;
	char input[10000];
	float total_counts=0;
	int	min_frag_len=10;
	
	
	
	ifstream InFile1, InFile3;
	ofstream OutFile1;

	InFile1.open(InCoords.c_str());
	if (!InFile1)
	{
		cout << "Error opening " << InCoords << ".\n";
		exit(1);
	}
	InFile3.open(ACC_file.c_str());
	if (!InFile3)
	{
		cout << "Error opening " << ACC_file << ".\n";
		exit(1);
	}

	OutFile1.open(OutMet.c_str(), ios_base::app);
	if (!OutFile1)  
	{
		cerr << "Error opening file '"<<OutMet<<"'.\n\n";
		exit(1);
	} 
	
	
	// Read in accessions
	int numACC = FILE_LINE_COUNT(ACC_file);
	
	string ACC[numACC];
	float ACC_reads[numACC];
	
	for(count=0; count < numACC; count++)
	{
		InFile3 >> ACC[count];
		ACC_reads[count]=0;
	}

	// Read in features
	cerr << "Reading in " << numFeatures << " genomic features..." << endl;
	read_features (InFeat, f, numFeatures);	
	
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
		
		InFile1 >> f[count].acc;		
		
		if(f[count].acc[0] == '#')
			InFile1.getline(input, 100000);
		else
			InFile1 >> s >> f[count].type >> f[count].start >> f[count].end >> s >> f[count].strand >> s >> f[count].tag;
		//cerr << "FEATURE " << f[count].acc << "	" << f[count].start << "	" << f[count].end << "	" << f[count].strand << "	" << f[count].tag << endl;
		
		f[count].counts_sense=0;
		f[count].counts_antisense=0;
		f[count].coord_name="";
		
		if(InFile1.eof())
			break;
		
	}
		
	InFile1.close();
	return 0;
}
//######################################
int overlap (string cacc, float cs, float ce, char cd, struct feature f[], int numFeatures, string cn)
{
	int clen;
	int flen;
	int count;
	float total_OL=0;
	float total_coord_q=0;
	
	for (count=0; count < numFeatures; count++)
	{
		if(cacc==f[count].acc)
		{
			f[count].read_overlap_sense=0;
			f[count].read_overlap_antisense=0;
			if(f[count].start > ce)
				break;	
			if(f[count].end  < cs)
				continue;	
			//cout << "	" << f[count].acc << "	" << f[count].start << "	" << f[count].end << "	" << f[count].strand << "	" << f[count].tag << endl;
		
			if(cd==f[count].strand)
			{	
				if (ce >= f[count].start && ce <= f[count].end && cs <= f[count].start)						// End of read is within feature
					f[count].read_overlap_sense=ce-f[count].start;
				if (cs >= f[count].start && cs <= f[count].end && ce >= f[count].end)						// Start of read is within feature
					f[count].read_overlap_sense=f[count].end-cs;
				if (cs > f[count].start && ce < f[count].end)										// Entire read is contained w/in feature
					f[count].read_overlap_sense=ce-cs;
				if (cs <= f[count].start && ce >= f[count].end)										// Entire feature is within read
					f[count].read_overlap_sense=f[count].end-f[count].start;
			
				if(f[count].read_overlap_sense > 0)
					total_OL=total_OL+f[count].read_overlap_sense;
			}
			if(cd!=f[count].strand)
			{	
				if (ce >= f[count].start && ce <= f[count].end && cs <= f[count].start)						// End of read is within feature
					f[count].read_overlap_antisense=ce-f[count].start;
				if (cs >= f[count].start && cs <= f[count].end && ce >= f[count].end)						// Start of read is within feature
					f[count].read_overlap_antisense=f[count].end-cs;
				if (cs > f[count].start && ce < f[count].end)										// Entire read is contained w/in feature
					f[count].read_overlap_antisense=ce-cs;
				if (cs <= f[count].start && ce >= f[count].end)										// Entire feature is within read
					f[count].read_overlap_antisense=f[count].end-f[count].start;
			
				if(f[count].read_overlap_antisense > 0)
					total_OL=total_OL+f[count].read_overlap_antisense;
			}
		}
	}

	for (count=0; count < numFeatures; count++)
	{
		if(cacc==f[count].acc)
		{
			if(f[count].read_overlap_sense > 0)
			{	
				f[count].counts_sense=f[count].counts_sense+f[count].read_overlap_sense/total_OL;
				total_coord_q=total_coord_q+(f[count].read_overlap_sense/total_OL);
			}
			if(f[count].read_overlap_antisense > 0)
			{	
				f[count].counts_antisense=f[count].counts_antisense+f[count].read_overlap_antisense/total_OL;
				total_coord_q=total_coord_q+(f[count].read_overlap_antisense/total_OL);
			}
			if(f[count].start > ce)
				break;	
			if(f[count].end  < cs)
				continue;	
		}
		
	}
	
	//if(total_coord_q == 0)
		//cerr << total_coord_q << "	" << cn << endl;
	//if(total_coord_q != 0)
		//cerr << total_coord_q << "	" << cn << endl;

	return 0;	
}

