#!/bin/bash
# call vevet assember for all AB mixed sample 5.5
echo Hellow world!
cd /home/ali/Downloads/velvet_1.2.10
for i in  {5..9};
do
  for j in {1..2};
do
    n=$((10-i))	
   a='AB.'$i'.'$n'_'$j
   b='AB.'$i'.'$n'_R1.fq\ '$j
   d='AB.'$i'.'$n'_R2.fq\ '$j

echo $b
echo $d
 
q='./velveth '$a' 21 -shortPaired -fastq -separate '$b' '$d
eval $q


p='./velvetg '$a' -ins_length 200 -ins_length_sd 30 -exp_cov auto -cov_cutoff auto  -scaffolding yes'
eval $p
done
done
