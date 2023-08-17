#!/bin/sh

BLAST=$1
COORD=$2


cat $BLAST | grep -v "#" | awk '{

									if($9 < $10)
									{	
										DIR="+"
										S=$9
										E=$10
									}
									if($9 > $10)
									{	
										DIR="-"
										S=$10
										E=$9
									}
										print $2, S, E, DIR, $1 
								}' > $COORD