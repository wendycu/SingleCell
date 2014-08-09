#!/bin/bash
#$ -cwd


list=$1

f1=`head -1 $list`
output=$list.counts.txt

cat $f1  > $output

for f in `tail -n +2 $list`
do
    cut -f2 $f | paste $output - > $output.temp
#     paste $output $f.c2 > $output.temp
    mv $output.temp $output 
done



