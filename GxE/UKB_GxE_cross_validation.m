%%% load data from UKB_GxE.m

%% %%%%%%%%%%%%%%%%%
%% cross-validation of interactions
%% %%%%%%%%%%%%%%%%
clear SNP_Cardio_int_P SNP_Cardio_int_T SNP_Cardio_int_SE SNP_P
SNP_Cardio_int_P=zeros(360, 22); SNP_Cardio_int_P(SNP_Cardio_int_P==0)=NaN;
SNP_Dep_int_P=zeros(360, 22); SNP_Dep_int_P(SNP_Dep_int_P==0)=NaN;
%%%% step 1 regen the maps for each fold
folds=1:3584:length(T.APOE);folds(11)=35846;
for i=1:10; ixfold(:,i)=zeros(1,35846); ixfold(folds(i):folds(i+1),i)=1; end
for i=1:360%[2:34, 36:68]%38:68
    %tmp=sortrows(unique(allsnps.rsID(allsnps.TraitNum==1)));
    tmp=sortrows(unique(allsnps.rsID(allsnps.TraitNum<100)));
    for j=1:length(tmp); ixx(:,j)=~cellfun('isempty',strfind(CT_SNPs.Properties.VariableNames, tmp{j})); end; ixx=sum(ixx')';
    snps=CT_SNPs(:,ixx'==1);
    %snps=CT_SNPs.rs9926320_G;%rs7875050_G;
    
for snp=1:length(snps.Properties.VariableNames)
Outcome=ct_ordered(:,i); Outcome=clean(Outcome, 4); 
APOE=snps(:, snp); APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
%APOE=snps; APOE=round(APOE);APOE(CT_SNPs.FID==0)=NaN;
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
parfor f=1:10
mdl = fitlm(T(ixfold(:,f)==0,:),'Outcome~APOE*cardio+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_Cardio_int_P(i,snp, f)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
%SNP_Cardio_int_T(i,snp)=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
SNP_Cardio_int_beta(i,snp, f)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
SNP_Cardio_int_SE(i,snp, f)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
%SNP_P(i,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));

mdl = fitlm(T(ixfold(:,f)==0,:),'Outcome~APOE*Dep+sex+age+site+age*sex+age^2+TIV+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_Dep_int_P(i,snp,f)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
%SNP_Dep_int_T(i,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
SNP_Dep_int_beta(i,snp,f)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
SNP_Dep_int_SE(i,snp,f)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
end
end
end
%%%% step 2 run the interactions
folds=1:3584:35846;fold(11)=35846;
for snp=1:length(snps.Properties.VariableNames)
    for f=1:10
    fdr=SNP_Cardio_int_P(:,snp,f);%mafdr(SNP_Cardio_int_P(:,snp,f), 'BHFDR', 'false');    %label_names_all(fdr<0.1)
    fdr2=SNP_Dep_int_P(:,snp,f);%mafdr(SNP_Dep_int_P(:,snp,f), 'BHFDR', 'false');    %label_names_all(fdr<0.1)
    if sum(fdr<0.05)==0; fdr=zeros(1,360); end
    if sum(fdr2<0.05)==0; fdr2=zeros(1,360); end
    try
    Outcome=nanmean(ct_ordered(:,fdr<0.05)')'; if sum(fdr<0.05)==1; Outcome=ct_ordered(:,fdr<0.05); end 
    APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    T=T(folds(f):folds(f+1),:);
    mdl = fitlm(T,'Outcome~APOE*cardio+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    SNP_CrossVal_Cardio_int_P(f,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
    SNP_CrossVal_Cardio_int_beta(f,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
    SNP_CrossVal_Cardio_int_SE(f,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
    catch 
    SNP_CrossVal_Cardio_int_P(f,snp)=NaN;    SNP_CrossVal_Cardio_int_beta(f,snp)=NaN;   SNP_CrossVal_Cardio_int_SE(f,snp)=NaN;
    end
    
         
    try
    Outcome=nanmean(ct_ordered(:,fdr2<0.05)')'; if sum(fdr2<0.05)==1; Outcome=ct_ordered(:,fdr2<0.05); end 
    APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    T=T(folds(f):folds(f+1),:);
    mdl = fitlm(T,'Outcome~APOE*Dep+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    SNP_Comb_Dep_int_P(f,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
    SNP_Comb_Dep_int_beta(f,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
    SNP_Comb_Dep_int_SE(f,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
    catch
    SNP_Comb_Dep_int_P(f,snp)=NaN;SNP_Comb_Dep_int_beta(f,snp)=NaN;SNP_Comb_Dep_int_SE(f,snp)=NaN;
    end;
end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for snp=1:220; [sp, fp]=meta_analyzep(SNP_CrossVal_Cardio_int_P(:,snp)); combinedp(snp,1:2)=[sp, fp];end
snps.Properties.VariableNames(combinedp(:,2)<0.01)
combinedp((combinedp(:,2)<0.01),2)
for snp=1:220; [sp, fp]=meta_analyzep(SNP_Comb_Dep_int_P(:,snp)); combinedp(snp,1:2)=[sp, fp];end
sum(mafdr(combinedp(:,2), 'BHFDR', 'true')<0.1)
combinedp(:,3)=mafdr(combinedp(:,2), 'BHFDR', 'true');
