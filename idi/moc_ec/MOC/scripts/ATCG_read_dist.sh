#!/bin/sh


for FILE in $@ff
do
	printf "%s\t" $FILE
	cat $FILE | head -100000 | awk -v totalA=0 -v totalC=0 -v totalG=0 -v totalT=0   '{

																										totalA=0
																										totalC=0
																										totalG=0
																										totalT=0
																										if(NR % 4 == 2)
																										{
																											l=length($1)
																											for(i=1; i<l+1; i++)
																											{
																												t=substr($1, i, 1)
																												if(t=="A")
																													totalA++
																												if(t=="C")
																													totalC++
																												if(t=="G")
																													totalG++
																												if(t=="T")
																													totalT++
																											}
																										
																											printf "%s\t%s\t%s\t%s\n", totalA/l*100,totalC/l*100,totalG/l*100,totalT/l*100
																										
																										}
																									}' > temp.txt
																			
																																												
	awk -v totalA=0 -v totalC=0 -v totalG=0 -v totalT=0 '{
														
														if($1 > 80)
															totalA++
														if($2 > 80)
															totalC++
														if($3 > 80)
															totalG++
														if($4 > 80)
															totalT++
														printf "%s\t%s\t%s\t%s\n", totalA/NR*100,totalC/NR*100,totalG/NR*100,totalT/NR*100
														
														
														}' temp.txt	| tail -1 
														
done					