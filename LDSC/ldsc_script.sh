## installing ldsc
# python -m virtualenv venvldsc
source venvldsc/bin/activate
#follow ldsc page https://github.com/bulik/ldsc
conda activate ldsc
module load bio/ldsc

sed -i -e "s/LogOR/beta/g" ./rg_META/previous_studies/PGC_UKB_depression_genome-wide.txt
./ldsc/munge_sumstats.py \
--out ./rg_META/previous_studies/dep \
--merge-alleles w_hm3.snplist \
--N 807553 \
--sumstats ./rg_META/previous_studies/PGC_UKB_depression_genome-wide.txt \
--snp MarkerName \
--a1 A1 \
--a2 A2 \
--p P

sed -i '1,73d' ./rg_META/previous_studies/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv

./ldsc/munge_sumstats.py \
--out ./rg_META/scz \
--merge-alleles w_hm3.snplist \
--N 130644 \
--sumstats ./rg_META/previous_studies/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv \
--snp ID \
--a1 A1 \
--a2 A2 \
--p PVAL


./ldsc/munge_sumstats.py \
--out ./rg_META/alz \
--merge-alleles w_hm3.snplist \
--N 74004 \
--sumstats ./rg_META/previous_studies/PGCALZ2sumstatsExcluding23andMe.txt \
--snp PosGRCh37 \
--a1 testedAllele \
--a2 otherAllele \
--p p


./ldsc/munge_sumstats.py \
--out ./rg_META/sesa \
--merge-alleles w_hm3.snplist \
--N 346978 \
--sumstats ./rg_META/previous_studies/SESA_neuro_clus_sumstats.txt \
--snp SNP \
--a1 A1 \
--a2 A2 \
--ignore Z \
--p P
 
./ldsc/munge_sumstats.py \
--out ./rg_META/adhd \
--merge-alleles w_hm3.snplist \
--N 53293 \
--sumstats ./rg_META/previous_studies/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.txt \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p P

./ldsc/munge_sumstats.py \
--out ./rg_META/int \
--merge-alleles w_hm3.snplist \
--N 225955 \
--sumstats ./rg_META/previous_studies/SavageJansen_2018_intelligence_metaanalysis.txt \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p P

./ldsc/munge_sumstats.py \
--out ./rg_META/insom \
--merge-alleles w_hm3.snplist \
--N-col N \
--sumstats ./rg_META/previous_studies/Insomnia_sumstats_Jansenetal.txt \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p P

./ldsc/munge_sumstats.py \
--out ./rg_META/ad19 \
--merge-alleles w_hm3.snplist \
--N-col Neff \
--sumstats ./rg_META/previous_studies/AD_sumstats_Jansenetal_2019sept.txt \
--snp SNP \
--a1 A1 \
--a2 A2 \
--ignore Z \
--p P


##### CT rgs
for i in `seq 2 34`; do
./ldsc/munge_sumstats.py \
--out ./rg_META/META_CT_sumstats/ct${i} \
--merge-alleles w_hm3.snplist \
--N 34000 \
--sumstats ./rg_META/METAANALYSIS_CT${i}_WB_EUR_1.txt \
--snp MarkerName \
--a1 Allele1 \
--a2 Allele2 \
--p P-value
done

for i in `seq 1 34`; do
echo -ne "rg_META/META_CT_sumstats/ct${i}.sumstats.gz,"
done


# run the regional rgs with AD stats
for i in `seq 1 34`; do
./ldsc/ldsc.py \
--rg rg_META/META_CT_sumstats/ct${i}.sumstats.gz,rg_META/adhd.sumstats.gz,rg_META/scz.sumstats.gz,rg_META/dep.sumstats.gz,rg_META/sesa.sumstats.gz,rg_META/int.sumstats.gz,rg_META/insom.sumstats.gz,rg_META/ad19.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ./rg_META/META_regionalCT_rgs/AD_ct${i}
done

cd /external/rprshnas01/tigrlab/scratch/pzhukovsky/biobank/genetics/LDSC/rg_META/META_regionalCT_rgs
for i in `seq 1 34`; do
tail AD_ct${i}.log -n 10 | head -n 7 >> sumstats.txt
done
