
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

%%% WMH interaction testing


clear SNP_P SNP_Cardio_WMH_P SNP_Cardio_WMH_T
for snp=1:length(snps.Properties.VariableNames)
Outcome=sqrt(AD_ordered.x25781_2_0); %Outcome=(AD_ordered.x25019_2_0+ AD_ordered.x25020_2_0)/2; Outcome=clean(Outcome, 4);
APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 

mdl = fitlm(T,'Outcome~APOE*cardio+sex+age+site+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_Cardio_WMH_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
SNP_Cardio_WMH_T(1,snp)=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
SNP_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));

mdl = fitlm(T,'Outcome~APOE*Dep+sex+age+site+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_Dep_int_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));

mdl = fitlm(T,'Outcome~APOE*AlcFreq+sex+age+site+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_AlcFreq_int_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:AlcFreq'));

mdl = fitlm(T,'Outcome~APOE*BMI+sex+age+site+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
% plotInteraction(mdl,'BMI','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_BMI_int_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:BMI'));

end
sum(mafdr(SNP_Cardio_WMH_P,'BHFDR', 'true')<0.1)
