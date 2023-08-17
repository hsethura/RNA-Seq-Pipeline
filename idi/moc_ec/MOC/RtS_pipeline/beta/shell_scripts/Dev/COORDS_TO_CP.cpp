#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <cstdlib>

using namespace std;
/*---------------------------------------------------------------------------
 This program takes in a file with start, stop coords and dir of transripts
 and outputs an array of numbers (1 per line) corresponding to the number of 
 of times that coordinates from  0-largest coord are represented 
 The first column corresponds to the + strand, the second to the - strand
 ----------------------------------------------------------------------------*/

#define NUM_COORDS 500000000
#define REP_SIZE 10000000

static string usage = "Usage: Coords_to_histogram <path to _coord.txt file> < path to .gff file > <min abundance (default 0)>  <strand specificity> < path to distribution file > [window size (default 1)]\n " ;

struct COORDS 
{
	string name;
	string acc;
	float start;
	float stop;
	char dir;
	char pair;
};


void coords_to_his(string, char, int, string);


int main(int argc, char *argv[])
{
	
	if(argc == 1)
	{
		cout << usage;
		exit(1);
	}
	
	string InCoords=argv[1];
	string ACC=argv[2];
	char s_spec=argv[3][0];			// set strand specificity flag 
	
	cerr << "Strand specific histogram? " << s_spec <<  "\n";
	int w=1;
	
	if(argc >= 8)
		w=atoi(argv[7]);		// set window size
	if(w%2 != 0)
		w=w-1;
	coords_to_his(InCoords, s_spec, w, ACC);

	return 0;
}

//******************************************************************

void 
coords_to_his(string InCoords, char s_spec, int w, string ACC)
{
	struct COORDS c;
	
	int count;
	int array=100000;
	char input[array];
	float temp1;
	float temp2;
	int k=0;
	int n;
	double max_coord=0;
	double total_abundance=0;
	
	double reads=0;
	double pair=0;
	double mapped=0;
	
	string x;
	char r;

	ifstream InFile1, InFile2;
	ofstream OutFile1;

	InFile1.open(InCoords.c_str());
	
	if (!InFile1)
	{
		cout << "Error opening " << InCoords << ".\n";
		exit(1);
	}

	double *p = (double *) calloc (REP_SIZE, sizeof(double));
	double *m = (double *) calloc (REP_SIZE, sizeof(double));

	
	//------------------- initialize abundance arrays -----------------------------//
	for(count=0; count < REP_SIZE; count++)
	{
		p[count]=0;
		m[count]=0;
	}
	
	//------------------- read in coords and increment abundance arrays -----------------------------//
	
	cerr << "Reading in coords from " << InCoords << endl;

	for(count=0; count < NUM_COORDS; count++)
	{
		InFile1 >> c.acc >> c.start >> c.stop >> c.dir >> c.pair;		
		//cerr << c.start << endl;
		InFile1.getline (input, 100);		
		
		if (c.acc==ACC)
		{
			if(c.pair=='Y')
			{	
				pair++;
				reads+=2;
			}
			else
				reads++;

			if(InFile1.eof())
			{	
				break;
			}

			if(c.stop > 0 && c.start > 0)
			{
				if (c.stop  > max_coord)
					max_coord = c.stop;
			
				total_abundance=total_abundance+(c.stop-c.start+1);
	
				if (c.dir == '+')
					for(k=c.start; k < c.stop+1; k++)
						p[k-1]+=1;
				if (c.dir == '-')
					for(k=c.start; k < c.stop+1; k++)
						m[k-1]+=1;
			
			}
			if(count == NUM_COORDS-1)
			{
		
				cerr << "TOO MANY COORDS!\n";
				exit(1);
			
			}
		}
	}
cerr << count << " coordinate pairs read in with maximum coord " << max_coord << endl;


//------------------- adjust data for window size -----------------------------//
	if(w > 0)	
	{	
		for(k=(w/2); k < max_coord-(w/2); k++)
		{
			temp1=0;
			temp2=0;
			for(n=(k-w/2); n < k+w/2; n++)
			{
				temp1=temp1+p[n];
				temp2=temp2+m[n];
			}
			p[k]=temp1/w;
			m[k]=temp2/w;
		}
	}
//------------------- output abundance arrays -----------------------------//
	cerr << "Writing out CP... " << endl;
	
	for(k=0; k < max_coord-(w/2); k++)
	{
		
		if(s_spec=='Y')
		{
			if(m[k] > 0  )
				m[k]=-m[k];
			cout << p[k] << " " << m[k] << endl;
		}
		
		
		if(s_spec=='N')
		{	
			p[k]=p[k]+m[k];
			cout << p[k] << endl;
		}
	}

}



//******************************************************************
