#!/bin/bash
# Split a PLINK file by chromosome chunks and send 23 jobs (23 chromosomes) and then recombine results
# Syntax: bash parallel_prune.sh <PLINK_filename> <output_filename>
# Example: bash parallel_prune.sh 06_duplicated_snps_rm 09_pruned
# suggested r^2 be 0.1
PLINK_filename=$1
output_filename=$2
r2=$3
export PATH=$PATH:/home/Bowei/ADHD_clinical
# Split files into 23 chunks

for i in {1..23}
do
 nohup plink --noweb --bfile ${PLINK_filename} --chr ${i} --make-bed --out ${PLINK_filename}_chr${i} &
done
for i in {1..23}
do
 waiton $(findpid "${PLINK_filename} --chr ${i} --make-bed")
done


# Send the pruning jobs
let i=1 a=23
    while [ $i -le $a ]
      do
        nohup plink --noweb --bfile ${PLINK_filename}_chr${i} --indep-pairwise 1500 100 ${r2} --out ${output_filename}_chr${i} &
        let i=$i+1
      done
for i in {1..23}
do
  waiton $(findpid "${PLINK_filename}_chr${i} --indep-pairwise")
done

cat ${output_filename}_chr*prune.in > ${output_filename}.prune.in
plink --noweb --bfile ${PLINK_filename} --extract ${output_filename}.prune.in --make-bed --out ${output_filename}

rm ${PLINK_filename}_chr*
rm ${output_filename}_chr*
