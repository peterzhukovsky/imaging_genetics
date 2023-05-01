clear
cd D:\Canada_2020\CLSA

clsa_base=readtable('D:\Canada_2020\CLSA\data\clsa_dru_bl_cog1.csv'); clsa_base=sortrows(clsa_base, 2);
clsa_snps=readtable('D:\Canada_2020\CLSA\data\CLSA_dosage.csv');clsa_snps=sortrows(clsa_snps, 2);
clsa_prs=readtable('D:\Canada_2020\CLSA\data\clsa_ct_prs.csv');clsa_prs(26623:end,:) = [];clsa_prs.ADM_GWAS_COM=(clsa_prs.FID);
clsa_snps=innerjoin(clsa_snps, clsa_prs);
clsa_snps(:,'ADM_GWAS_COM') = [];clsa_snps(:,'ADM_GWAS3') = [];
clsa_base.CCC_HBP_COM(clsa_base.CCC_HBP_COM==2)=0;
clsa_base.CCC_HEART_COM(clsa_base.CCC_HEART_COM>3)=NaN;clsa_base.CCC_HBP_COM(clsa_base.CCC_HBP_COM>3)=NaN; clsa_base.CCC_ANGI_COM(clsa_base.CCC_ANGI_COM>3)=NaN; clsa_base.CCC_AMI_COM(clsa_base.CCC_AMI_COM>3)=NaN;
clsa_base.DEP=clsa_base.DEP_CESD10_COM>10;
clsa_base.STP_INTFR_RATIO_COM(clsa_base.STP_INTFR_RATIO_COM>10)=NaN;


clsa_joint=innerjoin(clsa_base, clsa_snps); %clsa_joint.cardio=clsa_joint.CCC_HEART==1 | clsa_joint.CCC_HBP==1 | clsa_joint.CCC_AMI==1 | clsa_joint.CCC_ANGI==1;
clsa_joint.ED_UDR11_COM(clsa_joint.ED_UDR11_COM>40)=NaN;

regions_short={'GlobalMeanMean thickness';'STS';'cACC';'cMFG';'cuneus';'entorhinal';'fusiform';'iparietal';'itemporal';'isthmuscingulate';'loccipital';'lorbitofrontal';'lingual';'morbitofrontal';'mtemporal';'phippocampal';'pocentral';'popercularis';'porbitalis';'ptriangularis';'pcalcarine';'pcentral';'pcingulate';'precentral';'precuneus';'rACC';'rMFG';'SupF';'SupPar';'SupT';'SupMarg';'fPole';'trTemp';'insula'};

%clsa_joint.cardio=clsa_joint.CCC_HEART_COM==1 | clsa_joint.CCC_HBP_COM==1 | clsa_joint.CCC_AMI_COM==1 | clsa_joint.CCC_ANGI_COM==1;

%%
% for i in `seq 1 34`; do head -n 26622 CT${i}.all_score > CT${i}.all_score.txt; done
% for i in `seq 1 34`; do sed -ie 's/[ ][ ]*/ /g' CT${i}.all_score.txt ; done
clear snp_*
for i=1:34; i
prs=readtable(strcat('D:\Canada_2020\CLSA\data\PRS\CT', num2str(i),'.all_score.txt'), 'Delimiter', ' ');
prs.Properties.VariableNames={'FID', 'IID', 'Pt_5e08', 'Pt_5e07', 'Pt_5e06', 'Pt_5e05', 'Pt_0p0005', 'Pt_0p005', 'Pt_0p05', 'Pt_0p1', 'Pt_0p5', 'Pt_1', 'Missing'};
prs.ADM_GWAS_COM=prs.FID;
T=innerjoin(clsa_joint, prs);
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_5e08 + AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_main(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Pt_5e08'),:);
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_5e08*CCC_HBP_COM+ AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');%+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
    snp_stats_int(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'CCC_HBP_COM:Pt_5e08'),:);
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_5e08*DEP+ AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');%+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10
    snp_stats_int_DEP(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'DEP_1:Pt_5e08'),:);
    
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_5e05 + AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_main5e05(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Pt_5e05'),:);
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_0p005 + AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_main0p005(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Pt_0p005'),:);
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_5e07 + AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_main5e07(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Pt_5e07'),:);
    mdl=fitlm(T, 'COG_CONSTR_EF4_COM ~ Pt_5e06 + AGE_NMBR_COM + ED_UDR11_COM + SEX_ASK_COM + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_main5e06(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Pt_5e06'),:);
    
end

regions_short(mafdr(snp_stats_main.pValue, 'BHFDR', true)<0.1)
regions_short(mafdr(snp_stats_int.pValue, 'BHFDR', true)<0.1)
regions_short(mafdr(snp_stats_int_DEP.pValue, 'BHFDR', true)<0.1)

snp_stats_main(1,:)=[];snp_stats_int(1,:)=[];snp_stats_int_DEP(1,:)=[];
snp_stats_main.region=regions_short; snp_stats_main.pfdr=mafdr(snp_stats_main.pValue, 'BHFDR', true);
snp_stats_int.pfdr=mafdr(snp_stats_int.pValue,'BHFDR', false);
snp_stats_int_DEP.pfdr=mafdr(snp_stats_int_DEP.pValue, 'BHFDR', true);

i=3 %19 %20
figure; ix=T.hypert==0; X=prs_thick(ix,i); Y=clsa_joint.COG_CONSTR_MEM_COM(ix); %figure; scatter(X,Y);lsline;
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'b-','LineWidth',3); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'b-','LineWidth',1);
ix=T.hypert==1; X=prs_thick(ix,i); Y=clsa_joint.COG_CONSTR_MEM_COM(ix);  %figure; scatter(X,Y);lsline;
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'r-','LineWidth',3); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'r-','LineWidth',1);

i=29 %19 %20 %8 %29
figure; ix=T.hypert==0; X=prs_thick(ix,i); Y=clsa_joint.COG_CONSTR_EF4_COM(ix); %figure; scatter(X,Y);lsline;
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'b-','LineWidth',3); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'b-','LineWidth',1);
ix=T.hypert==1; X=prs_thick(ix,i); Y=clsa_joint.COG_CONSTR_EF4_COM(ix);  %figure; scatter(X,Y);lsline;
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'r-','LineWidth',3); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'r-','LineWidth',1);

i=3 %19 %20
T=table(clsa_joint.CCC_HBP_COM, clsa_joint.DEP, clsa_joint.entity_id,clsa_joint.COG_CONSTR_MEM_COM,clsa_joint.AGE_NMBR_COM,clsa_joint.ED_UDR11_COM,clsa_joint.SEX_ASK_COM, clsa_snps.PC1,clsa_snps.PC2,clsa_snps.PC3,clsa_snps.PC4,clsa_snps.PC5,clsa_snps.PC6,clsa_snps.PC7,clsa_snps.PC8,clsa_snps.PC9,clsa_snps.PC10, (prs_thick(:,i)),...
    'VariableNames', {'hypert', 'DEP', 'entity_id_clsa_long','COG_REYII_SCORE','AGE_BL','ED_UDR11','SEX_ASK','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','SNP'});
mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP+ AGE_BL + ED_UDR11 + SEX_ASK + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
figure; X=prs_thick(:,i); Y=clsa_joint.COG_CONSTR_MEM_COM;  %scatter(X,Y);lsline;
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'k-','LineWidth',3); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'k-','LineWidth',1);
figure; X=clsa_joint.ED_UDR11_COM; Y=clsa_joint.COG_CONSTR_MEM_COM; %scatter(X,Y);lsline;
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'k-','LineWidth',3); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'k-','LineWidth',1);


for i=1:33; i
    T=table(clsa_joint.CCC_HBP_COM, clsa_joint.DEP, clsa_joint.entity_id,clsa_joint.COG_REYII_SCORE_COM,clsa_joint.AGE_NMBR_COM,clsa_joint.ED_UDR11_COM,clsa_joint.SEX_ASK_COM, clsa_snps.PC1,clsa_snps.PC2,clsa_snps.PC3,clsa_snps.PC4,clsa_snps.PC5,clsa_snps.PC6,clsa_snps.PC7,clsa_snps.PC8,clsa_snps.PC9,clsa_snps.PC10, (prs_thick(:,i)),...
        'VariableNames', {'hypert', 'DEP', 'entity_id_clsa_long','COG_REYII_SCORE','AGE_BL','ED_UDR11','SEX_ASK','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','SNP'});
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP + AGE_BL + ED_UDR11 + SEX_ASK');
    snp_stats_main(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'SNP'),:);
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP*hypert + AGE_BL+ED_UDR11+SEX_ASK');
    snp_stats_int(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'hypert:SNP'),:);
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP*DEP + AGE_BL+ED_UDR11+SEX_ASK');
    snp_stats_int_DEP(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'DEP_1:SNP'),:);
end
regions_short(snp_stats_main.pValue<0.0017)
regions_short(snp_stats_int.pValue<0.05)
regions_short(snp_stats_int_DEP.pValue<0.01)


%% viz PRS interactions and main effects!
i=30 %1 5 19 30
corr(clsa_joint.COG_REYI_NORMED_ZSCORE_COM,prs_thick(:,i), 'rows','pairwise')
figure; scatter(prs_thick(:,i)*-1,clsa_joint.COG_REYII_SCORE_COM)
figure; scatter(prs_thick(:,i)*-1,clsa_joint.COG_REYI_NORMED_ZSCORE_COM);lsline

corr(clsa_joint.COG_REYI_NORMED_ZSCORE_COM,prs_thick(:,i), 'rows','pairwise')
corr(clsa_joint.COG_CONSTR_OVERALLCOG6_COM,prs_thick(:,i), 'rows','pairwise')
figure; scatter(prs_thick(:,i)*-1,clsa_joint.COG_CONSTR_OVERALLCOG6_COM);lsline

i=18 % 8 iparietal 16 phippo 18 popercularis interaction
T=table(clsa_joint.CCC_HBP_COM, clsa_joint.DEP, clsa_joint.entity_id,clsa_joint.COG_REYII_SCORE_COM,clsa_joint.AGE_NMBR_COM,clsa_joint.ED_UDR11_COM,clsa_joint.SEX_ASK_COM, clsa_snps.PC1,clsa_snps.PC2,clsa_snps.PC3,clsa_snps.PC4,clsa_snps.PC5,clsa_snps.PC6,clsa_snps.PC7,clsa_snps.PC8,clsa_snps.PC9,clsa_snps.PC10, (prs_thick(:,i)),...
    'VariableNames', {'hypert', 'DEP', 'entity_id_clsa_long','COG_REYII_SCORE','AGE_BL','ED_UDR11','SEX_ASK','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','SNP'});
mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP*hypert + AGE_BL+ED_UDR11+SEX_ASK')

figure; ix=T.hypert==2; %scatter(prs_thick(ix,i)*-1,clsa_joint.COG_REYII_NORMED_ZSCORE_COM(ix));lsline ;xlim([-4 2]);
X=prs_thick(ix,i)*-1; Y=clsa_joint.COG_REYII_NORMED_ZSCORE_COM(ix); 
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'b-'); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'b-','LineWidth',2);
ix=T.hypert==1; %scatter(prs_thick(ix,i)*-1,clsa_joint.COG_REYII_NORMED_ZSCORE_COM(ix));lsline ;xlim([-4 2]);
X=prs_thick(ix,i)*-1; Y=clsa_joint.COG_REYII_NORMED_ZSCORE_COM(ix); 
m=fitlm(X,Y);  [ypred,ci] = predict(m,X);hold on; plot(X,ypred,'r-'); hold on; tmp=[X, ci]; tmp=sortrows(tmp,1); plot(tmp(:,1),tmp(:,2:3),'r-','LineWidth',2);

%% viz {'rs11126806'}

 T=table(clsa_joint.CCC_HBP_COM, clsa_joint.DEP, clsa_joint.entity_id,clsa_joint.COG_REYII_SCORE_COM,clsa_joint.AGE_NMBR_COM,clsa_joint.ED_UDR11_COM,clsa_joint.SEX_ASK_COM, clsa_snps.PC1,clsa_snps.PC2,clsa_snps.PC3,clsa_snps.PC4,clsa_snps.PC5,clsa_snps.PC6,clsa_snps.PC7,clsa_snps.PC8,clsa_snps.PC9,clsa_snps.PC10, round(clsa_snps.rs11126806),...
        'VariableNames', {'hypert', 'DEP', 'entity_id_clsa_long','COG_REYII_SCORE','AGE_BL','ED_UDR11','SEX_ASK','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','SNP'});
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP*hypert + AGE_BL + ED_UDR11 + SEX_ASK+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 ')%+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10

pred=mdl.predict; 
[means,predci] = grpstats(pred(T.hypert==0), T.SNP(T.hypert==0), {'mean','sem'},'Alpha',0.05);
figure(1); hold off; errorbar(means,3.291.*predci); hold on
[means,predci] = grpstats(pred(T.hypert==1), T.SNP(T.hypert==1), {'mean','sem'},'Alpha',0.05);
figure(1); errorbar(means,3.291.*predci); hold off

pred=clsa_joint.COG_REYI_NORMED_ZSCORE_COM; figure(1);violinplot(pred(T.hypert==1), T.SNP(T.hypert==1)); %ylim([0 10]); 
figure(2);violinplot(pred(T.hypert==2), T.SNP(T.hypert==2)); %ylim([0 10]);

