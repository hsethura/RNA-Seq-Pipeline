#!/bin/sh


FASTQ=$1


awk '{
				if(NR % 4 == 1 || NR % 4 == 3)
					print $0
				
				if(NR % 4 == 2)
				{
					
					len=length($1)
					for(i=1; i < len+1; i++)
					{
						y=substr($1, i, 1)
						if(y != "G")
						{
							p=i
							break
						}
					
					}
					t=substr($1, p, len+1)
					print t
				}
				if(NR % 4 == 0)	
				{
					t=substr($1, p, len+1)
					print t
				}
							
		}' $FASTQ 
		
		
