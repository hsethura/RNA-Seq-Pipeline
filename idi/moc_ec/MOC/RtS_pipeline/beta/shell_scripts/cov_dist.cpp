#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <cstdlib>

#include "SIPHTseq_info.h" 

using namespace std;

struct GENE
{
	string name;
	int	start;
	int end;
	char dir;
	float len; 
};

struct hist
{
	double p;
	double m;
};

static string usage = "Usage: coverage_bin <path to GENE file> <path to HISTO file> < # of bins >";

int main(int argc, char *argv[])
{
	if(argc < 2)
	{
		cout << usage;
		exit(1);
	}
	
	string InGENE=argv[1];			
	string InHISTO=argv[2];			
	int BIN_NUM=atoi(argv[3]);			
	int ADD_5=atoi(argv[4]);	
	int ADD_3=atoi(argv[5]);	

	int count, n, k, i, r;
	float bin_len=0;
	double total_abund;
	double min_total=0;
	double bin_abund;
	char input[10000];
	int num=0;
	
	ifstream InFile1, InFile2;

	InFile1.open(InGENE.c_str());

	if (!InFile1)
	{
		cout << "Error opening " << InGENE << ".\n";
		exit(1);
	}
	InFile2.open(InHISTO.c_str());
	
	if (!InFile2)
	{
		cout << "Error opening " << InHISTO << ".\n";
		exit(1);
	}
	
	int GENE_NUM = FILE_LINE_COUNT(InGENE);
	int HISTO_SIZE = FILE_LINE_COUNT(InHISTO);
	
	struct GENE c[GENE_NUM];
	
	double *p = (double *) calloc (HISTO_SIZE, sizeof(double));
	double *m = (double *) calloc (HISTO_SIZE, sizeof(double));
	double temp;
	
	float b[BIN_NUM];
	
	//------------------- read in coords  ----------------------------//
	
	//GENE_NUM=3;
	
	for(count=0; count < GENE_NUM; count++)
	{
		InFile1 >> c[count].name >> c[count].start >> c[count].end >> c[count].dir;
		
		if(c[count].dir=='+')
		{
			c[count].start=c[count].start-ADD_5;
			c[count].end=c[count].start+ADD_3;
		}
		if(c[count].dir=='-')
		{
			c[count].start=c[count].start-ADD_3;
			c[count].end=c[count].start+ADD_5;
		}

		
		if(InFile1.eof())
		{	
			break;
		}
		c[count].len=c[count].end-c[count].start+1;
		
		//cerr <<  c[count].name << "\t" << c[count].start << "\t" << c[count].dir << "\t" << count << "\t" << GENE_NUM << endl;
	}
	
	//------------------- input historgram -----------------------------//
	for(count=1; count < HISTO_SIZE+1; count++)
	{
		InFile2 >> p[count] >> m[count];
		//cerr << p[count] << " " << m[count]<< " " << count <<  endl;

		temp=m[count];
		m[count]=-(temp);
		if(InFile2.eof())
			break;
	} 
	cerr << "Abundance at " << count << " positions read in" << endl;
	
	
	//------------------- calculate bin coverage for each gene -----------------------------//

	
	for(count=0; count < GENE_NUM; count++)
	{
		
		//------------------- initialize bin coverage array -----------------------------/
		//cerr << "initializing bin coverage array for " << c[count].name  << endl;
		for(n=0; n < BIN_NUM; n++)
			b[n]=0;
		
		bin_len=(c[count].len/BIN_NUM);

		//------------------- calculate total coverage of gene -----------------------------//
		//cerr << "calculating total coverage of gene " << endl;

		total_abund=0;
		
		
		for(k=c[count].start; k < c[count].end+1; k++)
		{
			if(c[count].dir=='+')
				total_abund=total_abund+p[k];
			if(c[count].dir=='-')
				total_abund=total_abund+m[k];
		}
		
		//cerr << total_abund << endl;
		
		//------------------- calculate bin coverages -----------------------------//
		//cerr << " calculating bin coverages " << endl;
		if(total_abund > min_total)
		{
			num++;
			i=0;
			r=0;
			bin_abund=0;
			
			//cerr << bin_len << "\t" << BIN_NUM << endl;
			cout << c[count].name << "\t" << total_abund << "\t";		
			if(c[count].dir=='+')
			{
				for(k=c[count].start; k < c[count].end+1; k++)
				{
					bin_abund=bin_abund+p[k];
					r++;
					if(r >= ((i+1)*bin_len))
					{
						b[i]=bin_abund/total_abund*100;
						cout << b[i] << "\t";
						bin_abund=0;
						i++;
					}
					if(k == c[count].end && i < BIN_NUM)
					{
						b[i]=bin_abund/total_abund*100;
						cout << b[i] << "\t";
						bin_abund=0;
						i++;
					}
					
				}
			}
			if(c[count].dir=='-')
			{
				for(k=c[count].end; k > c[count].start-1; k--)
				{
					bin_abund=bin_abund+m[k];
					r++;
					if(r >= ((i+1)*bin_len))
					{
						b[i]=bin_abund/total_abund*100;
						cout << b[i] << "\t";
						bin_abund=0;
						i++;
					}
					if(k == c[count].start && i < BIN_NUM)
					{
						b[i]=bin_abund/total_abund*100;
						cout << b[i] << "\t";
						bin_abund=0;
						i++;
					}
				}
			}
			
			//cerr << c[count].name << "\t" << c[count].start << "\t" << c[count].end  << "\t" << total_abund << endl;
			//for(n=0; n < BIN_NUM; n++)
				//cerr << n << "|" << b[n] << "\t";
			cout << endl;
			
		}
	}
	cerr << "A total of  " << num << " genes have been processed..." << endl;

	return 0;
}





