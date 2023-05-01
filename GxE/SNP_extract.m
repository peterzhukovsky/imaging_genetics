
%% Extracting SNPS with main effects for variant prioritization for interaction testing
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_WMH_META\GenomicRiskLoci.txt'));
snps.Trait=repmat({'WMH'}, [length(snps.rsID),1]);snps.TraitNum=100*ones(length(snps.rsID),1); allsnps=snps;
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_hippo_META\GenomicRiskLoci.txt'));
snps.Trait=repmat({'Hippo'}, [length(snps.rsID),1]); snps.TraitNum=101*ones(length(snps.rsID),1); allsnps=vertcat(allsnps, snps);
regions_short={'GlobalMeanMean thickness';'STS';'cACC';'cMFG';'cuneus';'entorhinal';'fusiform';'iparietal';'itemporal';'isthmuscingulate';'loccipital';'lorbitofrontal';'lingual';'morbitofrontal';'mtemporal';'phippocampal';'pocentral';'popercularis';'porbitalis';'ptriangularis';'pcalcarine';'pcentral';'pcingulate';'precentral';'precuneus';'rACC';'rMFG';'SupF';'SupPar';'SupT';'SupMarg';'fPole';'trTemp';'insula'};
for i=1:34
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_CT', num2str(i),'_META\GenomicRiskLoci.txt'));
snps.Trait=repmat(regions_short(i), [length(snps.rsID),1]); snps.TraitNum=i*ones(length(snps.rsID),1);
allsnps=vertcat(allsnps, snps);
end


for chr=1:length(unique(allsnps.chr))
    i=unique(allsnps.chr); i=i(chr);
rsids=allsnps.rsID(allsnps.chr==i);writetable(table(rsids), strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\snpextract\snpids_sign', num2str(i), '.txt'), 'WriteVariableNames',0);%,  'FileType', 'text', 'Delimiter', '\t');
commands{chr, :}=strcat('plink2 --pfile $MDIR/ukb/genomics/genomic_data_imp/pgen_white/chr', num2str(i), ' --extract <(cat snpids_sign', num2str(i),'.txt) --export A --out $SCRATCH/biobank/genetics/EURnonWB/METAL/snpextract/chr', num2str(i), 'SNPSshort --keep $SCRATCH/biobank/genetics/samples_to_extract_snps.txt');
end


for chr=1:length(unique(allsnps.chr))
    i=unique(allsnps.chr); i=i(chr);
try movefile(strcat('chr', num2str(i), 'SNPSshort.raw'), strcat('chr', num2str(i), 'SNPS.raw.txt')); end
dta=readtable(strcat('chr', num2str(i), 'SNPS.raw.txt'));
if chr==1
alldta=dta;
else
    dta(:,1:6)=[];
alldta=[alldta, dta];
end
end
alldta=sortrows(alldta,1);
datanames_rs=readtable('D:\Canada_2020\UK_biobank\reports\subjects_rs.txt');
datanames_fs=readtable('D:\Canada_2020\UK_biobank\reports\subjects_fs.txt');
datanames_all=unique(vertcat(datanames_fs.Var1, datanames_rs.Var1));
CT_SNPs=table;
for i=1:length(datanames_all)
    try
    if sum(datanames_all(i)==alldta.FID)>0
        CT_SNPs(i,:)=alldta((datanames_all(i)==alldta.FID),:); 
    else 
        CT_SNPs(i,:)=NaN;
    end
    catch
    end
end

save('allsnps.mat','CT_SNPs')

