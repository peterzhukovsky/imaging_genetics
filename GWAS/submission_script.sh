#!/usr/bin/env bash

#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=13
#BSATCH --mem-per-cpu=5
#SBATCH --time=20:00:00
#SBATCH --export=ALL
#SBATCH --job-name="gwas_pz"
#SBATCH --output=/KIMEL/tigrlab/scratch/pzhukovsky/biobank/genetics/CT/lbl_reg_%j.out
#SBATCH --error=/KIMEL/tigrlab/scratch/pzhukovsky/biobank/genetics/CT/lbl_reg_%j.err


module load bio/PLINK2
#PLINK/1.9b_6.15-x86_64
roi=X_SUBJECT_X

export MDIR=/external/rprshnas01/kcni/mwainberg
plink2 --out $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/king \
--keep $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/samples_to_include.txt \
--pfile $MDIR/ukb/genomic_data/full_genome_merged/ukb_final \
--make-king triangle bin \
--make-king-table \
--threads 12

chromosomes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"  
for i in $chromosomes; do 
export MDIR=/external/rprshnas01/kcni/mwainberg 
plink2 --out $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/chr${i}_gwas_sumstats \
--keep $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/samples_to_include.txt \
--pfile $MDIR/ukb/genomic_data_imp/pgen_white/chr${i} \
--pheno $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/phenotype.phe \
--pheno-name Outcome \
--pheno-quantile-normalize \
--covar $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/phenotype.phe \
--covar-name age sex site TIV PC1-PC10 \
--covar-variance-standardize \
--glm skip-invalid-pheno firth-fallback hide-covar omit-ref no-x-sex \
--king-cutoff $SCRATCH/biobank/genetics/EURnonWB/EUR/CT/$roi/king 0.0884 \
--threads 12 
done  



wait
