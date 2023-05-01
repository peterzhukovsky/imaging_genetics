
%% cog analysis by snp
%%% what are the y variables: 
%%% word pairs (PAL),  
Y=cognitive.x20197_2_0;

for snp=1:length(snps.Properties.VariableNames)
    Outcome=cognitive_ordered.x20197_2_0; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    mdl = fitlm(T,'Outcome~APOE+sex+age+site+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    SNP_PAL_int_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
    SNP_PAL_int_beta(1,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
    SNP_PAL_int_SE(1,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
end;

[SNP_PAL_int_Pfdr]=mafdr(SNP_DSST_int_P', 'BHFDR', false);

t=table(snps.Properties.VariableNames', SNP_PAL_int_beta', SNP_PAL_int_SE',SNP_PAL_int_P', 'VariableNames', {'tmp', 'beta','se', 'p'});
t=innerjoin(tmp, t);
