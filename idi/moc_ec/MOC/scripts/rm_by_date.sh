#!/bin/sh

DIR=$1
max_period=$2


now=`date +'%s'`

#find $DIR -printf '%p %C@\n' | awk '{split($2, ar, "."); period=('$now'-ar[1])/(60*60*24); if(period >= '$max_period')print $1, period}'

find $DIR -type d -empty