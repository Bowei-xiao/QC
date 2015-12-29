# !/bin/bash

#QC script
# would assume this is the binary format
plinkOri=$1
# We want noPedID.txt to indicate whom we want to throw away
prefix=1
scriptADS=/hpf/projects/arnold/bowei/ADHD_ClinicalData
export PATH=$PATH:${scriptADS}
#get Missingness
# remove Individuals that we couldn't find pedigree ID
plink --noweb --bfile ${plinkOri} --remove noPedID.txt --make-bed --out ${prefix}_${plinkOri}_withPID 
plink --noweb --bfile ${prefix}_${plinkOri}_withPID --missing --out ${prefix}_${plinkOri}
#Fix the blank; change delimiter to single blank
awk '{gsub(" ?","",$1)}1' ${prefix}_${plinkOri}.imiss > ${prefix}_${plinkOri}_blankFix.imiss
awk '{gsub(" ?","",$1)}1' ${prefix}_${plinkOri}.lmiss > ${prefix}_${plinkOri}_blankFix.lmiss
# Filter Inds with missingness >= 3%
awk '$6 <= 0.03' ${prefix}_${plinkOri}_blankFix.imiss | cut -f1,2 > ${prefix}_${plinkOri}IndLowMiss.txt
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
plink --noweb --bfile ${plinkOri} --keep ${prefix_old}_${plinkOri}IndLowMiss.txt --make-bed --out ${prefix}_${plinkOri}_highMissIndRMV
# remove SNPs with call rate less than 95%
awk '$5 <= 0.05' ${prefix_old}_${plinkOri}_blankFix.lmiss | cut -f2 > ${prefix}_${plinkOri}_highCallSNP.txt
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
plink --noweb --bfile ${prefix_old}_${plinkOri}_highMissIndRMV --extract ${prefix_old}_${plinkOri}_highCallSNP.txt --make-bed --out ${prefix}_${plinkOri}_highMissIndSNPRMV

#Heterozygosity check
#pruning:
#prefix_old=$prefix
#prefix=$(( $prefix + 1 ))
plink --noweb --bfile ${prefix}_${plinkOri}_highMissIndSNPRMV --maf 0.05 --make-bed --out ${prefix}_${plinkOri}_commonSNP
cp ${scriptADS}/waiton .
cp ${scriptADS}/findpid .
bash ${scriptADS}/parallel_prune.sh  ${prefix}_${plinkOri}_commonSNP ${prefix}_${plinkOri}_commonSNPAfterPrune 0.2

# Calculate heterozygosity for autosomal and X separately
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonSNPAfterPrune --het --out ${prefix}_${plinkOri}_autoHet
# Chromosome X:
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonSNPAfterPrune --chr 23 --make-bed --out ${prefix}_${plinkOri}_chrX
sed -i 's/^23/1/g' ${prefix}_${plinkOri}_chrX.bim
plink --noweb --bfile ${prefix}_${plinkOri}_chrX --het --out ${prefix}_${plinkOri}_sexHet
# Check Heteozygousity and removes individuals that failed it (with plots)
Rscript ${scriptADS}/1118HetCheck.R ${prefix}_${plinkOri}
# remove people that failed Heteozygousity check
plink --noweb --bfile ${prefix_old}_${plinkOri}_highMissIndSNPRMV --remove ${prefix}_${plinkOri}_indNeedremove.txt --make-bed --out ${prefix_old}_${plinkOri}_Hetremove 
# only keep common SNPs (MAF>=0.05) since we don't have a lot of ppl
plink --noweb --bfile ${prefix_old}_${plinkOri}_Hetremove --maf 0.05 --make-bed --out ${prefix_old}_${plinkOri}_commonHetremove
# Sex check
prefix=$(( $prefix + 1 )) #prefix changed to 5, but prefix_old is still 3
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonHetremove --check-sex --out ${prefix}_${plinkOri}_sexCheck
#impute sex
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonHetremove --make-bed --impute-sex --out ${prefix}_${plinkOri}_sexImpute

plink --noweb --bfile ${prefix}_${plinkOri}_sexImpute --check-sex --out ${prefix}_${plinkOri}_sexCheck2
# Removing heterozygous SNPs in chromosome X:
cut -f3 ${prefix}_${plinkOri}_sexCheck2.hh > ${prefix}_${plinkOri}_SNPtoRemoveHH
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
plink --noweb --bfile ${prefix_old}_${plinkOri}_sexImpute --exclude ${prefix_old}_${plinkOri}_SNPtoRemoveHH --make-bed --out ${prefix}_${plinkOri}_hetSNPremoved

# Hardy-Weinberg Equilibrium
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
plink --noweb --bfile ${prefix_old}_${plinkOri}_hetSNPremoved --hardy --make-bed --out  ${prefix}_${plinkOri}_HWE

# summarize the final result files:
cp ${prefix}_${plinkOri}_HWE.bim ${plinkOri}_afterQC.bim
cp ${prefix}_${plinkOri}_HWE.bed ${plinkOri}_afterQC.bed
cp ${prefix}_${plinkOri}_HWE.fam ${plinkOri}_afterQC.fam 
# recording individuals that were removed in the QC
awk '$6 > 0.03 && $1 + 0 == $1' 1_${plinkOri}_blankFix.imiss | cut -d' ' -f1,2 > TMP_lowRate.txt
yes 'Low call Rate' | head -n `wc -l TMP_lowRate.txt | cut -d' ' -f1` > TMP_lowRateReason.txt
yes 'Failed heteo check' | head -n `wc -l 4_${plinkOri}_indNeedremove.txt | cut -d' ' -f1` > TMP_heteoReason.txt
yes 'No pedigree ID found' | head -n `wc -l noPedID.txt | cut -d' ' -f1` > TMP_noPedReason.txt
paste TMP_lowRate.txt TMP_lowRateReason.txt > TMP_2.txt
paste noPedID.txt TMP_noPedReason.txt > TMP_1.txt
paste 4_${plinkOri}_indNeedremove.txt TMP_heteoReason.txt > TMP_3.txt
cat TMP_1.txt TMP_2.txt TMP_3.txt > ${plinkOri}_indRemovedInQC.txt

rm TMP_*



