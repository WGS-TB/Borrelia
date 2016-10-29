#!/bin/bash


////////////////////////////////////////////////////////////
directory_rawdata="$1table.txt"

#$1 is the directory of the input folder 
ls $1 > $directory_rawdata
readarray names < $directory_rawdata
nlen=${#names[@]}-1

for (( i=0; i<${nlen}; i++ ));
do

   a=${names[$i]} 
   b=${a:0:$len-1}
    
    
   dir2="$1$b/zztable"
   echo $dir2
   ls "$1$b" > $dir2


   readarray names_f < $dir2
    
   n1len=${#names_f[@]}-1


  
 
   j=0
      
	 q="Rscript --vanilla Wrapper.R $1$b/${names_f[$j]}  $1$b/${names_f[$j+1]} $2 $1$b/leftOver"
         eval $q

  
  

done




/////////////////////////////////////////////////////////////


for (( i=0; i<${nlen}; i++ ));
do

   a=${names[$i]} 
   b=${a:0:$len-1}
    
    
   dir2="$1$b/leftOver/zztable"
   echo $dir2
   ls "$1$b/leftOver" > $dir2


   readarray names_c < $dir2
    
   n1len=${#names_c[@]}-1


  
   j=0

   q="./interleave_fastq.sh $1$b/leftOver/${names_c[$j]}  $1$b/leftOver/${names_c[$j+1]}  > $1$b/leftOver/interleave"
   eval $q
   

  q="VelvetOptimiser.pl --d $1$b/leftOver/Results/assembly_result -s 25 -e 51 --k 'n50' --c 'tbp' --f '-fastq -shortPaired $1$b/leftOver/interleave' "
  eval $q

   q="Rscript LongContig.R $1$b/leftOver/Results/assembly_result/"
   eval $q
  

done


