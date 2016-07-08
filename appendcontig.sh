#!/bin/bash
# append velvet out-put to  correspoding reference for running needleman aligner
echo Hellow world!
cd /home/ali/Downloads/velvet_1.2.10

num="1"
for i in  {5..9};
do
  for j in {1..2};
do
    n=$((10-i))	
   a='AB.'$i'.'$n'_'$j
  
      
   if (( j = 1))
   then
	b='cat /home/ali/Desktop/reference/ref1 >> /home/ali/Downloads/velvet_1.2.10/'$a'/contigs.fa'
        eval $b
   else
	b='cat /home/ali/Desktop/reference/ref2 >> /home/ali/Downloads/velvet_1.2.10/'$a'/contigs.fa'
        eval $b
   fi

 

done
done
