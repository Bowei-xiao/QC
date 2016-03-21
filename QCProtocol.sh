# !/bin/bash

#QC script
# Syntax: bash QCProtocol.sh ${PlinkFileName} 0.03 0.05
# would assume this is the binary plink format
plinkOri=$1
indRate=$2 #Individual missing rate, I normally use 0.03 for Ilumina; and 0.05 for Affy6.0
snpRate=$3 #SNP missing rate, I normally use 0.05
prefix=1
# We want noPedID.txt to indicate whom we want to throw away
# If there is no such thing, then just skip this step
# remove Individuals that we couldn't find pedigree ID
# It will looks for files noPedID.txt (should be in usual plink two columns format)
if [ -e "noPedID.txt" ]
then
plink --noweb --bfile ${plinkOri} --remove noPedID.txt --make-bed --out ${prefix}_${plinkOri}_withPID 
else
cp ${plinkOri}.fam ${prefix}_${plinkOri}_withPID.fam
cp ${plinkOri}.bim ${prefix}_${plinkOri}_withPID.bim
cp ${plinkOri}.bed ${prefix}_${plinkOri}_withPID.bed
fi
# Firstly removing chromosome 24 and above
awk '$1 >23 && $1!="X"{print $2}' ${prefix}_${plinkOri}_withPID.bim > ${prefix}_${plinkOri}_notUseChr.txt 
plink --noweb --bfile ${prefix}_${plinkOri}_withPID --exclude ${prefix}_${plinkOri}_notUseChr.txt --make-bed --out ${prefix}_${plinkOri}Chr123

#missing rate check
plink --noweb --bfile ${prefix}_${plinkOri}Chr123 --missing --out ${prefix}_${plinkOri}
#Fix the blank; change delimiter to single blank
awk '{gsub(" ?","",$1)}1' ${prefix}_${plinkOri}.imiss > ${prefix}_${plinkOri}_blankFix.imiss
awk '{gsub(" ?","",$1)}1' ${prefix}_${plinkOri}.lmiss > ${prefix}_${plinkOri}_blankFix.lmiss
# Filter Inds with missingness >= the preset rate
awk -v miss="$indRate" '$6 <= miss' ${prefix}_${plinkOri}_blankFix.imiss | cut -f1,2 > ${prefix}_${plinkOri}IndLowMiss.txt
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
###prefix=2 pre_old=1 
plink --noweb --bfile ${plinkOri} --keep ${prefix_old}_${plinkOri}IndLowMiss.txt --make-bed --out ${prefix}_${plinkOri}_highMissIndRMV
# remove SNPs with call rate that are more than pre-set rate
awk -v miss="$snpRate" '$5 <= miss' ${prefix_old}_${plinkOri}_blankFix.lmiss | cut -f2 > ${prefix}_${plinkOri}_highCallSNP.txt
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
###prefix=3 pre_old=2
plink --noweb --bfile ${prefix_old}_${plinkOri}_highMissIndRMV --extract ${prefix_old}_${plinkOri}_highCallSNP.txt --make-bed --out ${prefix}_${plinkOri}_highMissIndSNPRMV

#Heterozygosity check
#First do maf filter + prune the data
plink --noweb --bfile ${prefix}_${plinkOri}_highMissIndSNPRMV --maf 0.05 --make-bed --out ${prefix}_${plinkOri}_commonSNP

# This is an old way to do it. Only use this part with plink1.07. We could directly prune it with plink 1.9
if false; then
bash ${scriptADS}/parallel_prune.sh  ${prefix}_${plinkOri}_commonSNP ${prefix}_${plinkOri}_commonSNPAfterPrune 0.2
fi
# For plink 1.9; we can directly use this
plink --noweb --bfile ${prefix}_${plinkOri}_commonSNP --indep-pairwise 1500 100 ${r2} --out ${prefix}_${plinkOri}_commonSNPAfterPrune 
plink --noweb --bfile ${prefix}_${plinkOri}_commonSNP --extract ${prefix}_${plinkOri}_commonSNPAfterPrune.prune.in --make-bed --out ${prefix}_${plinkOri}_commonSNPAfterPrune

# Calculate heterozygosity for autosomal and X separately
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
###prefix=4 pre_old=3
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonSNPAfterPrune --het --out ${prefix}_${plinkOri}_autoHet
# Chromosome X:
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonSNPAfterPrune --chr 23 --make-bed --out ${prefix}_${plinkOri}_chrX
sed -i 's/^23/1/g' ${prefix}_${plinkOri}_chrX.bim
plink --noweb --bfile ${prefix}_${plinkOri}_chrX --het --out ${prefix}_${plinkOri}_sexHet
# Check Heteozygousity and mark individuals that failed it (with plots)
Rscript ${scriptADS}/1118HetCheck.R ${prefix}_${plinkOri}

# Sex check
prefix=$(( $prefix + 1 )) #prefix changed to 5, but prefix_old is still 3
###prefix=5 pre_old=3
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonSNP --check-sex --out ${prefix}_${plinkOri}_sexCheck
#impute sex
plink --noweb --bfile ${prefix_old}_${plinkOri}_commonSNP --make-bed --impute-sex --out ${prefix}_${plinkOri}_sexImpute

plink --noweb --bfile ${prefix}_${plinkOri}_sexImpute --check-sex --out ${prefix}_${plinkOri}_sexCheck2
# Removing heterozygous SNPs in chromosome X:
cut -f3 ${prefix}_${plinkOri}_sexCheck2.hh > ${prefix}_${plinkOri}_SNPtoRemoveHH
prefix_old=$prefix
prefix=$(( $prefix + 1 ))
###prefix=6 pre_old=5
plink --noweb --bfile ${prefix_old}_${plinkOri}_sexImpute --exclude ${prefix_old}_${plinkOri}_SNPtoRemoveHH --make-bed --out ${prefix}_${plinkOri}_hetSNPremoved

# Hardy-Weinberg Equilibrium (Mark SNPs but do not remove)
plink --noweb --bfile ${prefix}_${plinkOri}_hetSNPremoved --hardy --out ${prefix}_${plinkOri}_hwe

# summarize the final result files:
cp ${prefix}_${plinkOri}_hetSNPremoved.bim ${plinkOri}_afterQC.bim
cp ${prefix}_${plinkOri}_hetSNPremoved.bed ${plinkOri}_afterQC.bed
cp ${prefix}_${plinkOri}_hetSNPremoved.fam ${plinkOri}_afterQC.fam 
# recording individuals that were removed or marked in the QC
# $1+=0 == $1 used to remove header, it's not perfect, can only deal with header strating with Characers not numbers(But ours is fine since it's "FID")
awk -v miss="$indRate" '$6 > miss && $1 + 0 == $1' 1_${plinkOri}_blankFix.imiss | cut -d' ' -f1,2 > TMP_lowRate.txt
yes 'Low call Rate(Removed)' | head -n `wc -l TMP_lowRate.txt | cut -d' ' -f1` > TMP_lowRateReason.txt
yes 'Failed Heterozygosity check(not Removed)' | head -n `wc -l 4_${plinkOri}_indNeedremove.txt | cut -d' ' -f1` > TMP_FailedHetReason.txt
paste TMP_lowRate.txt TMP_lowRateReason.txt > TMP_2.txt
paste 4_${plinkOri}_indNeedremove.txt TMP_FailedHetReason.txt > TMP_3.txt
if [ -e "noPedID.txt" ]
then
yes 'No pedigree ID found' | head -n `wc -l noPedID.txt | cut -d' ' -f1` > TMP_noPedReason.txt
paste noPedID.txt TMP_noPedReason.txt > TMP_1.txt
cat TMP_1.txt TMP_2.txt TMP_3.txt > ${plinkOri}_indMarkedInQC.txt
else 
cat TMP_2.txt TMP_3.txt > ${plinkOri}_indMarkedInQC.txt
fi
# Recording SNPs that were removed/marked in the QC(rare SNPs not recorded)
yes 'Not on Chr1-23(Removed)' | head -n `wc -l 1_${plinkOri}_notUseChr.txt | cut -d' ' -f1` > TMP_SNPnotON.txt
paste ${prefix}_notUseChr.txt TMP_SNPnotON.txt > TMP_SNP1.txt
awk -v miss="$snpRate" '$5 > miss && $1 + 0 == $1' 1_${plinkOri}_blankFix.lmiss | cut -d' ' -f2 > TMP_SNPlowRate.txt
yes 'Low call Rate(Removed)' | head -n `wc -l TMP_SNPlowRate.txt | cut -d' ' -f1` > TMP_SNPlowRateReason.txt
paste TMP_SNPlowRate.txt TMP_SNPlowRateReason.txt > TMP_SNP2.txt
yes 'Heterozygous SNPs on ChrX for Male(Removed)' | head -n `wc -l 5_${plinkOri}_SNPtoRemoveHH | cut -d' ' -f1` > TMP_SNPHetReason.txt
paste 5_${plinkOri}_SNPtoRemoveHH TMP_SNPHetReason.txt > TMP_SNP3.txt
awk '$3 == "ALL" && $9 <= 0.0001{print $2}' ${prefix}_${plinkOri}_hwe.hwe >TMP_hweSNP.txt
yes 'Failed HWE(Not Removed)' | head -n `wc -l TMP_hweSNP.txt | cut -d' ' -f1` > TMP_SNPhweReason.txt
paste TMP_hweSNP.txt TMP_SNPhweReason.txt > TMP_SNP4.txt
cat TMP_1.txt TMP_2.txt TMP_3.txt TMP_4.txt> ${plinkOri}_SNPMarkedInQC.txt


rm TMP_*
mkdir ${plinkOri})intermediateSteps
mv [1-7]_* > ${plinkOri})intermediateSteps


