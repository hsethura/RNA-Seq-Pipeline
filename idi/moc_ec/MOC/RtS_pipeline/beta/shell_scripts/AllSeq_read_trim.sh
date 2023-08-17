#!/bin/bash

FASTQ=$1
OUTFILE=$2
MAX_A=$3
MIN_LEN=$4

awk -v max_A=$MAX_A -v MIN_LEN=$MIN_LEN '{
    
    if(NR % 4 == 1)
    {	
        name1=$1
        name2=$2
    }

    if(NR % 4 == 2)
    {
    
        len=length($1)

        last_end = 0
        A_count = 0
        for (i = 1; i <= len; i++) {
            lchar = substr($1, i, 1)
            if (lchar == "A") {
                A_count++
                if (A_count == max_A) {
                    # Time to trim
                    break
                }
            } else {
                # Reset
                A_count = 0
                last_end = i
            }
        }

        
        if(last_end >= MIN_LEN)
		{	
			if(len != last_end)
				print name1":trimmed_"last_end + 1"_"(len - last_end)" "name2
			else 
				print name1":untrimmed_0_0 "name2

			printf "%s\n", substr($1, 1, last_end)
		}
    }
	
	if(last_end >= MIN_LEN)
	{	
		if(NR % 4 == 3) {
			print $1
		}

		if(NR % 4 == 0)
		{	
			printf "%s\n", substr($1, 1, last_end)
		}
	}

}' $FASTQ > $OUTFILE

echo $OUTFILE

			

