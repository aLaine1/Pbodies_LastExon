#!/bin/bash

bed_file=$1
join_depth=$2
split_depth=$3
vectorPB=$4
vectorCYTO=$5

while read p; do
  chr=$(echo $p |cut -f1 -d' ')
  coord1=$(echo $p |cut -f2 -d' ')
  coord2=$(echo $p |cut -f3 -d' ')
  name=$(echo $p |cut -f4 -d' ')
  echo "$name" >> $split_depth
  cat $join_depth | awk -v chr=$chr -v coord1=$coord1 -v coord2=$coord2 -v name=$name '$1 == chr && $2 >= coord1 && $2 <= coord2 {a=($3+$4+$5)/3;b=($6+$7+$8)/3; print a"\t"b}' >> $split_depth
done < $bed_file

sed -i 's/,/\./g' $split_depth

vectorUp=""
vectorDown=""
while read p; do
  if [[ $p == ENST* ]];
  then
    echo $vectorUp >> $vectorPB
    echo $vectorDown >> $vectorCYTO
    vectorUp=$p
    vectorDown=$p
  else
    Up=$(echo $p |cut -f1 -d' ')
    Down=$(echo $p |cut -f2 -d' ')
    vectorUp=$vectorUp" "$Up
    vectorDown=$vectorDown" "$Down
  fi
done < $split_depth
