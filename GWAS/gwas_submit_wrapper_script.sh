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


### METAL section

for i in `seq 1 34`; do 
cp $SCRATCH/biobank/genetics/EURnonWB/METAL/metal_script.txt $SCRATCH/biobank/genetics/EURnonWB/METAL/metal_script${i}.txt
sed -i -e "s/X_roi_X/${i}/g" $SCRATCH/biobank/genetics/EURnonWB/METAL/metal_script${i}.txt
done

for i in `seq 1 34`; do 
metal $SCRATCH/biobank/genetics/EURnonWB/METAL/metal_script${i}.txt
done

for i in `seq 1 34`; do 
gzip METAANALYSIS_CT${i}_WB_EUR_1.txt
done
