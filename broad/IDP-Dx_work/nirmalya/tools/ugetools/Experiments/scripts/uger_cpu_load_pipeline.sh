#! /bin/bash
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

set -e
set -o pipefail

progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [ -d duration-minutes ] [ -c num-floats ] [ -r random-seed ] [ -s sort-chunk-size ] [-t]

Run a pipeline to measure UGER CPU utilization:
generate_random_numbers | sqrt(abs(number) | sort_chunks_of_numbers | sum numbers

-d duration-minutes: generate random numbers for this many minutes.  Default: no time limit.
-c num-floats:       generate this many numbers.  Default: no limit.
-r random-seed:      RNG seed.  Default: no seed.
-s sort-chunk-size:  Sort numbers in chunks of this size.  Default: program default.
-t:                  Use taskset to force processes to different CPUs

EOF
}

generate_args=
sort_args=
taskset=0
while getopts ":d:c:r:s:t" options; do
  case $options in
    d ) generate_args="$generate_args --duration $OPTARG";;
    c ) generate_args="$generate_args --count $OPTARG";;
    r ) generate_args="$generate_args --seed $OPTARG";;
    s ) sort_args="--chunk $OPTARG";;
    t ) taskset=1;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;

  esac
done
shift $(($OPTIND - 1))

if [ $# != 0 ]
then echo "ERROR: Extra arguments on command line" >&2
     usage
     exit 1
fi

script_dir=`dirname $0`

if (( $taskset ))
then T1='taskset 0x1'
     T2='taskset 0x2'
     T3='taskset 0x4'
     T4='taskset 0x8'
else T1=
     T2=
     T3=
     T4=
fi
$T1 $script_dir/generate_random_numbers.py $generate_args | \
    $T2 $script_dir/square_root.py | \
    $T3 $script_dir/sort_float_chunks.py $sort_args | \
    $T4 $script_dir/sum_floats.py

echo "Pipeline done"

