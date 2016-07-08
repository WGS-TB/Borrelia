#!/bin/bash
# align velvet out-puts to reference strain
echo Hellow world!
cd /home/ali/Desktop/needleman_wunsch-0.3.5


for i in  {5..9};
do
  for j in {1..2};
do
    n=$((10-i))	
   a='AB.'$i'.'$n'_'$j

  b='./needleman_wunsch --printscores --colour --file /home/ali/Downloads/velvet_1.2.10/'$a'/contigs.fa'
  echo $a
  eval $b
done
done
