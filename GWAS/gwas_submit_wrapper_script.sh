module load bio/PLINK2
#PLINK/1.9b_6.15-x86_64
module load bio/METAL


for i in `seq 2 34`; do 
cp $SCRATCH/biobank/genetics/EURnonWB/WB/sbatch_gwas.sh $SCRATCH/biobank/genetics/EURnonWB/WB/sbatch_gwas${i}.sh
sed -i -e "s/X_SUBJECT_X/${i}/g" $SCRATCH/biobank/genetics/EURnonWB/WB/sbatch_gwas${i}.sh
done

for i in `seq 2 34`; do 
sbatch --mem=80G -n 1 -c 13 $SCRATCH/biobank/genetics/EURnonWB/WB/sbatch_gwas${i}.sh
done
