
#!/bin/bash


awk -F, '{ if (NR>1) { print $1 }}' Betas.csv > rsidlist.txt

awk -F, '{ if (NR>1) { print sprintf("%02d", $2)":"$3"-"$3 }}' Betas.csv > chrposlist.txt

cmd=""
for i in {1..22}
do
  bgenix -g ukb_imp_chr${i}_v3.bgen -incl-rsids rsidlist.txt -incl-range chrposlist.txt > chr_${i}.bgen
  cmd=$cmd"chr_${i}.bgen "
done

# Combine the .bgen files for each chromosome into one
cat-bgen -g  $cmd -og initial_chr.bgen -clobber
# Write index file .bgen.bgi
bgenix -g initial_chr.bgen -index -clobber

# Remove the individual chromosome files
for i in {1..22}
do
  rm chr_${i}.bgen
done

# Import the betas into the sqlite database as a table called Betas
sqlite3 initial_chr.bgen.bgi "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator "," initial_chr.bgen.bgi ".import Betas.csv Betas"

sqlite3 initial_chr.bgen.bgi "DROP TABLE IF EXISTS Joined;"
# And inner join it to the index table (Variants), making a new table (Joined)
# By joining on alleles as well as chromosome and position 
# we can ensure only the relevant alleles from any multi-allelic SNPs are retained
sqlite3 -header -csv initial_chr.bgen.bgi \
"CREATE TABLE Joined AS 
  SELECT Variant.*, Betas.chr_name, Betas.Beta FROM Variant INNER JOIN Betas 
    ON Variant.chromosome = printf('%02d', Betas.chr_name) 
    AND Variant.position = Betas.chr_position 
    AND Variant.allele1 = Betas.noneffect_allele 
    AND Variant.allele2 = Betas.effect_allele 
  UNION 
  SELECT Variant.*, Betas.chr_name, -Betas.Beta FROM Variant INNER JOIN Betas 
    ON Variant.chromosome = printf('%02d', Betas.chr_name) 
    AND Variant.position = Betas.chr_position 
    AND Variant.allele1 = Betas.effect_allele AND 
    Variant.allele2 = Betas.noneffect_allele;"

# Filter the .bgen file to include only the alleles specified in the Betas for each SNP 
bgenix -g initial_chr.bgen -table Joined  > single_allelic.bgen

# And produce an index file for the new .bgen
bgenix -g single_allelic.bgen -index


plink2 --bgen single_allelic.bgen ref-first \
--hard-call-threshold 0.1 \
--sample ukbA_imp_chrN_v3_sP.sample \
--memory 15000 \
--set-all-var-ids @:#_\$r_\$a \
--freq \
--make-pgen \
--out raw 


awk '/^[^#]/ { if( $5>0.4 && $5<0.6 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' \
raw.afreq > exclrsIDs_ambiguous.txt

########################################################
#----------------------- SNP QC -----------------------#
########################################################

for i in {1..22}
do
  awk -v chr=$i 'BEGIN {FS="\t"; OFS="\t"} { print chr,$0,chr":"$3"_"$4"_"$5 }' ukb_mfi_chr${i}_v3.txt
done > ukb_mfi_all_v3.tsv


plink2 --pfile raw \
--memory 15000 \
--exclude exclrsIDs_ambiguous.txt \
--extract-col-cond ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 \
--maf 0.005 \
--write-snplist \
--make-pgen \
--out snpQC 

########################################################
#--------------------- Sample QC ----------------------#
########################################################

plink2 --pfile raw \
--memory 15000 \
--extract snpQC.snplist \
--keep-fam usedinpca.txt \
--write-samples \
--out sampleQC 

########################################################
#----------------- Calculate the PRS ------------------#
########################################################

sqlite3 -separator " " -list initial_chr.bgen.bgi \
"SELECT chr_name || ':' || position || '_' ||allele1 || '_' || allele2, allele2, Beta FROM Joined;" > score.txt


plink2 --pfile raw \
--memory 15000 \
--extract snpQC.snplist \
--keep sampleQC.id \
--score score.txt no-mean-imputation \
--out PRS

