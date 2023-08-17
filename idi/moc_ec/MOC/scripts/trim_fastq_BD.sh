
#!/bin/sh

INPUT_FILE=$1
TRIM_FILE=$2
TRIM_LEN=$3


	
zcat $INPUT_FILE | awk '{
								# qual
								if(NR%4==0)
								{
									y=substr($0,'$TRIM_LEN',10000)
									print y
								}
								
								# name
								if(NR%4==1)
									print $0
								
								# seq
								if(NR%4==2)
								{
									y=substr($0,'$TRIM_LEN',10000)
									print y
								}
								# +
								if(NR%4==3)
									print $0
							
							}' > $TRIM_FILE
							
							"Gzipping $TRIM_FILE"
							gzip $TRIM_FILE
							


