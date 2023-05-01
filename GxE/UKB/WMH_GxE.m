
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

%% %%% variance explained
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
mdl = fitlm(T,'Outcome~APOE+sex+age+site+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10')
r1=mdl.Rsquared.Ordinary;
mdl = fitlm(T,'Outcome~sex+age+site+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10')
(r1-mdl.Rsquared.Ordinary)*100
