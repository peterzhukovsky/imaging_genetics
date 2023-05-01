
%% rerunning interactions in 360 ROI space
clear SNP_Cardio_int_P SNP_Cardio_int_T SNP_Cardio_int_SE SNP_P
SNP_Cardio_int_P=zeros(360, 22); SNP_Cardio_int_P(SNP_Cardio_int_P==0)=NaN;
SNP_Dep_int_P=zeros(360, 22); SNP_Dep_int_P(SNP_Dep_int_P==0)=NaN;
for i=1:360%
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

mdl = fitlm(T,'Outcome~APOE*cardio+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_Cardio_int_P(i,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
%SNP_Cardio_int_T(i,snp)=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
SNP_Cardio_int_beta(i,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
SNP_Cardio_int_SE(i,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
%SNP_P(i,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));

mdl = fitlm(T,'Outcome~APOE*Dep+sex+age+site+age*sex+age^2+TIV+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
SNP_Dep_int_P(i,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
%SNP_Dep_int_T(i,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
SNP_Dep_int_beta(i,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
SNP_Dep_int_SE(i,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));

%mdl = fitlm(T,'Outcome~APOE*AlcFreq+sex+age+HMotion+site+age*sex+age^2+TIV+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
%SNP_AlcFreq_int_P(i,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:AlcFreq'));

%mdl = fitlm(T,'Outcome~APOE*BMI+sex+age+HMotion+site+age*sex+age^2+TIV+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
% plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %hold on; scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
%SNP_BMI_int_P(i,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:BMI'));

end
end

%snp2 fdr<0.1
%""lh.L_V1_ROI.label" "lh.L_MST_ROI.label" "lh.L_V2_ROI.label" "lh.L_V3_ROI.label" "lh.L_V4_ROI.label" "lh.L_FEF_ROI.label" "lh.L_55b_ROI.label" "lh.L_V3A_ROI.label" "lh.L_RSC_ROI.label" "lh.L_V3B_ROI.label" "lh.L_MT_ROI.label" "lh.L_7m_ROI.label" "lh.L_d23ab_ROI.label" "lh.L_5m_ROI.label" "lh.L_24dd_ROI.label" "lh.L_SCEF_ROI.label" "lh.L_1_ROI.label" "lh.L_6d_ROI.label" "lh.L_8Ad_ROI.label" "lh.L_8BL_ROI.label" "lh.L_8C_ROI.label" "lh.L_IFJa_ROI.label" "lh.L_9-46d_ROI.label" "lh.L_47s_ROI.label" "lh.L_OP4_ROI.label" "lh.L_AVI_ROI.label" "lh.L_ProS_ROI.label" "lh.L_IP1_ROI.label" "lh.L_IP0_ROI.label" "lh.L_PGi_ROI.label" "lh.L_PGs_ROI.label" "lh.L_31pd_ROI.label" "lh.L_LBelt_ROI.label"
%snp12
%"lh.L_45_ROI.label" "lh.L_PHA1_ROI.label" "lh.L_STSda_ROI.label" "lh.L_STSvp_ROI.label" "lh.L_pOFC_ROI.label"
% snp rs9926320_G

snp=strcmp(snps.Properties.VariableNames,'rs11126806_C');  %rs11126806_C %rs77690628_T %rs9926320_G %rs983960_C
fdr=mafdr(SNP_Cardio_int_P(:,snp), 'BHFDR', 'true');
label_names_all(fdr<0.1)
Outcome=nanmean(ct_ordered(:,fdr<0.1)')'; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
%T(T.cardio==1,:)=[];
mdl = fitlm(T,'Outcome~APOE*cardio+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
mdl
r1=mdl.Rsquared.Ordinary;
mdl = fitlm(T,'Outcome~APOE+cardio+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
(r1-mdl.Rsquared.Ordinary)*100

%figure;plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
pred=mdl.predict;
[means,predci] = grpstats(pred(T.cardio==0), T.APOE(T.cardio==0), {'mean','sem'},'Alpha',0.05);
figure(1); hold off; errorbar(means,3.291.*predci); hold on
[means,predci] = grpstats(pred(T.cardio==1), T.APOE(T.cardio==1), {'mean','sem'},'Alpha',0.05);
figure(1); errorbar(means,3.291.*predci); hold off

figure(1);violinplot(pred(T.cardio==0), T.APOE(T.cardio==0)); %ylim([2.5 2.8]); 
figure(2);violinplot(pred(T.cardio==1), T.APOE(T.cardio==1)); %ylim([2.5 2.8]);
figure(3);violinplot(pred, T.APOE); %ylim([2.5 2.8]);

snp=strcmp(snps.Properties.VariableNames,'rs7575796_A');  
fdr=mafdr(SNP_Dep_int_P(:,snp), 'BHFDR', 'true');
label_names_all(fdr<0.1)
Outcome=nanmean(ct_ordered(:,fdr<0.1)')'; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 

mdl = fitlm(T,'Outcome~Dep*APOE+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
mdl
%figure;plotInteraction(mdl,'cardio','APOE','predictions'); ylabel('Cortical Thickness'); %scatter(T.APOE(T.cardio==0), T.Outcome(T.cardio==0),2, 'filled', 'r'); hold on; scatter(T.APOE(T.cardio==1), T.Outcome(T.cardio==1),2, 'filled', 'b');
pred=mdl.predict; 
[means,predci] = grpstats(pred(T.Dep==0), T.APOE(T.Dep==0), {'mean','sem'},'Alpha',0.05);
figure(1); hold off; errorbar(means,3.291.*predci); hold on
[means,predci] = grpstats(pred(T.Dep==1), T.APOE(T.Dep==1), {'mean','sem'},'Alpha',0.05);
figure(1); errorbar(means,3.291.*predci); hold off

pred=mdl.predict; figure(1);violinplot(pred(T.Dep==0), T.APOE(T.Dep==0)); %ylim([2.5 2.8]); 
figure(2);violinplot(pred(T.Dep==1), T.APOE(T.Dep==1)); %ylim([2.5 2.8]);

%%
tmp=sortrows(unique(allsnps.rsID(allsnps.TraitNum<100)));
for j=1:length(tmp); ixx(:,j)=~cellfun('isempty',strfind(CT_SNPs.Properties.VariableNames, tmp{j})); end; ixx=sum(ixx')';
snps=CT_SNPs(:,ixx'==1);
    
for snp=1:length(snps.Properties.VariableNames)
    fdr=mafdr(SNP_Cardio_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)
    fdr2=mafdr(SNP_Dep_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)

    try
    Outcome=nanmean(ct_ordered(:,fdr<0.1)')'; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    mdl = fitlm(T,'Outcome~APOE*cardio+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    SNP_Comb_Cardio_int_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
    SNP_Comb_Cardio_int_beta(1,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
    SNP_Comb_Cardio_int_SE(1,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:cardio'));
    
    catch 
    SNP_Comb_Cardio_int_P(1,snp)=NaN;
    SNP_Comb_Cardio_int_beta(1,snp)=NaN;
    SNP_Comb_Cardio_int_SE(1,snp)=NaN;
    end
    
         
    try
    Outcome=nanmean(ct_ordered(:,fdr2<0.1)')'; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    mdl = fitlm(T,'Outcome~APOE*Dep+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    SNP_Comb_Dep_int_P(1,snp)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
    SNP_Comb_Dep_int_beta(1,snp)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
    SNP_Comb_Dep_int_SE(1,snp)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Properties.RowNames,'APOE:Dep_1'));
    catch
    SNP_Comb_Dep_int_P(1,snp)=NaN;
    SNP_Comb_Dep_int_beta(1,snp)=NaN;
    SNP_Comb_Dep_int_SE(1,snp)=NaN;
    
end;end

figure;histogram(log10(SNP_Comb_Cardio_int_P), 8)
snps.Properties.VariableNames(SNP_Comb_Cardio_int_P<0.05/220)
figure;histogram(log10(SNP_Comb_Dep_int_P), 8)
snps.Properties.VariableNames(SNP_Comb_Dep_int_P<0.05/220)
SNP_Comb_Dep_int_P(SNP_Comb_Dep_int_P<0.05/220)
        %plot all snps with sign interactions for cardio
figure;imagesc((log10(SNP_Comb_Cardio_int_P(SNP_Comb_Cardio_int_P<0.05/220)))); colormap bone; 
set(gca, 'XTick', (1:length(SNP_Comb_Cardio_int_P<0.05/220)), 'XTickLabel', snps.Properties.VariableNames(SNP_Comb_Cardio_int_P<0.05/220), 'XTickLabelRotation',90);
rosmapSNP_Comb_Cardio_int_P=[0.1	0.1	0.1	0.0006	0.1	0.1	0.0137	0.0034	0.0106	0.1	0.0152	0.1	0.0011	0.018	0.018	0.018	0.018	0.1	0.003	0.0006	0.1	0.0206	0.0007	0.1];
figure;imagesc(log10(rosmapSNP_Comb_Cardio_int_P)); colormap bone; 
        %plot all snps with sign interactions for Dep
figure;imagesc((log10(SNP_Comb_Dep_int_P(SNP_Comb_Dep_int_P<0.05/220)))); colormap bone; set(gca, 'XTick', (1:length(SNP_Comb_Dep_int_P<0.05/220)), 'XTickLabel', snps.Properties.VariableNames(SNP_Comb_Cardio_int_P<0.05/220), 'XTickLabelRotation',90);

%%% plotting the t-stats for each subgroup individually
for snp=1:length(snps.Properties.VariableNames)
    if SNP_Comb_Cardio_int_P(snp)<0.05/220
    fdr=mafdr(SNP_Cardio_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)
    fdr2=mafdr(SNP_Dep_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)
    
    Outcome=nanmean(ct_ordered(:,fdr<0.1)')'; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    T(T.cardio==1,:)=[];
    mdl = fitlm(T,'Outcome~APOE+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    tmpC1=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    T(T.cardio==0,:)=[];
    mdl = fitlm(T,'Outcome~APOE+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    tmpC0=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
    
    figure(33); scatter(tmpC0,tmpC1, 'filled','k'); hold on; text(tmpC0+0.05,tmpC1+rand(1,1)/10, snps.Properties.VariableNames(snp)); hold on
    end 
end

for snp=1:length(snps.Properties.VariableNames)
    if SNP_Comb_Dep_int_P(snp)<0.05/220
    fdr=mafdr(SNP_Dep_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)
        
    Outcome=nanmean(ct_ordered(:,fdr<0.1)')'; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    T(T.Dep==1,:)=[];
    mdl = fitlm(T,'Outcome~APOE+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    tmpC1=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
    T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    T(T.Dep==0,:)=[];
    mdl = fitlm(T,'Outcome~APOE+sex+age+site+age*sex+TIV+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10'); %HMotion
    tmpC0=mdl.Coefficients.tStat(strcmp(mdl.Coefficients.Properties.RowNames,'APOE'));
    figure(36); scatter(tmpC0,tmpC1, 'filled','k'); hold on; text(tmpC0+0.05,tmpC1+rand(1,1)/10, snps.Properties.VariableNames(snp)); hold on
    end 
end;clear tmp1 tmp0
%%%%% combine the sign results into a table 

sumstats_cardio=table;  sumstats_dep=table;  
for snp=1:length(snps.Properties.VariableNames)
    if SNP_Comb_Cardio_int_P(snp)<0.05/220
        tmp=table;
        fdr=mafdr(SNP_Cardio_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)
        tmp.beta=SNP_Cardio_int_beta(fdr<0.1, snp);
        tmp.SE=SNP_Cardio_int_SE(fdr<0.1, snp);
        tmp.P=SNP_Cardio_int_P(fdr<0.1, snp);
        tmp.Pfdr=fdr(fdr<0.1);
        tmp.region=label_names_all(fdr<0.1);
        s=snps.Properties.VariableNames(snp);
        tmp.snpname=repmat(s, [length(tmp.P), 1]);
        tmp.Overallbeta=repmat(SNP_Comb_Cardio_int_beta(snp), [length(tmp.P), 1]);
        tmp.OverallSE=repmat(SNP_Comb_Cardio_int_SE(snp), [length(tmp.P), 1]);
        tmp.OverallP=repmat(SNP_Comb_Cardio_int_P(snp), [length(tmp.P), 1]);
        sumstats_cardio=vertcat(sumstats_cardio, tmp);
    end
    if SNP_Comb_Dep_int_P(snp)<0.05/220
        tmp=table;
        fdr=mafdr(SNP_Dep_int_P(:,snp), 'BHFDR', 'true');    %label_names_all(fdr<0.1)
        tmp.beta=SNP_Dep_int_beta(fdr<0.1, snp);
        tmp.SE=SNP_Dep_int_SE(fdr<0.1, snp);
        tmp.P=SNP_Dep_int_P(fdr<0.1, snp);
        tmp.Pfdr=fdr(fdr<0.1);
        tmp.region=label_names_all(fdr<0.1);
        s=snps.Properties.VariableNames(snp);
        tmp.snpname=repmat(s, [length(tmp.P), 1]);
        tmp.Overallbeta=repmat(SNP_Comb_Dep_int_beta(snp), [length(tmp.P), 1]);
        tmp.OverallSE=repmat(SNP_Comb_Dep_int_SE(snp), [length(tmp.P), 1]);
        tmp.OverallP=repmat(SNP_Comb_Dep_int_P(snp), [length(tmp.P), 1]);
        sumstats_dep=vertcat(sumstats_dep, tmp);
    end
end
clear s tmp
