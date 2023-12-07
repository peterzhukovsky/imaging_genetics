clear
cd D:\Canada_2020\CLSA

clsa_base=readtable('D:\Canada_2020\CLSA\data\clsa_dru_bl_cog1.csv');
clsa_snps=readtable('D:\Canada_2020\CLSA\data\CLSA_dosage.csv');
clsa_snps(:,'ADM_GWAS3') = [];clsa_snps(30098:end,:) = [];
clsa_base.CCC_HEART_COM(clsa_base.CCC_HEART_COM>3)=NaN;clsa_base.CCC_HBP_COM(clsa_base.CCC_HBP_COM>3)=NaN; clsa_base.CCC_ANGI_COM(clsa_base.CCC_ANGI_COM>3)=NaN; clsa_base.CCC_AMI_COM(clsa_base.CCC_AMI_COM>3)=NaN;
clsa_base.CCC_HBP_COM(clsa_base.CCC_HBP_COM==2)=0;
clsa_base.DEP=clsa_base.DEP_CESD10_COM>10;
clsa_base.STP_INTFR_RATIO_COM(clsa_base.STP_INTFR_RATIO_COM>10)=NaN;

clsa_joint=outerjoin(clsa_base, clsa_snps);
clsa_joint.cardio=clsa_joint.CCC_HEART_COM==1 | clsa_joint.CCC_HBP_COM==1 | clsa_joint.CCC_AMI_COM==1 | clsa_joint.CCC_ANGI_COM==1;
%histogram(clsa_joint.AGE_NMBR_COM)


%% 
%% USE THIS MODEL!!! Overall memory/ overall cognition + Covary for PCs - dont covary for age sex education
%%
%%
clear snp_*
for i=1:215; i
    T=table(clsa_joint.CCC_HBP_COM, clsa_joint.DEP, clsa_joint.entity_id_clsa_base,clsa_joint.COG_CONSTR_MEM_COM,clsa_joint.AGE_NMBR_COM,clsa_joint.ED_UDR11_COM,clsa_joint.SEX_ASK_COM, clsa_joint.PC1, clsa_joint.PC2, clsa_joint.PC3, clsa_joint.PC4, clsa_joint.PC5, clsa_joint.PC6, clsa_joint.PC7, clsa_joint.PC8, clsa_joint.PC9, clsa_joint.PC10, round(clsa_joint{:,i+62}),...
        'VariableNames', {'hypert', 'DEP', 'entity_id_clsa_long','COG_REYII_SCORE','AGE_BL','ED_UDR11','SEX_ASK','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'SNP'});
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP + AGE_BL + SEX_ASK +AGE_BL^2 + AGE_BL*SEX_ASK+ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_main(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'SNP'),:);
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP*hypert + AGE_BL + SEX_ASK +AGE_BL^2 + AGE_BL*SEX_ASK++ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_int(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'hypert:SNP'),:);
    mdl=fitlm(T, 'COG_REYII_SCORE ~ SNP*DEP + AGE_BL + SEX_ASK +AGE_BL^2 + AGE_BL*SEX_ASK++ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10');
    snp_stats_int_DEP(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'DEP_1:SNP'),:);
end

snp_stats_main.pfdr=mafdr(snp_stats_main.pValue, 'BHFDR', 'true')
snp_stats_main.snp_names=clsa_joint.Properties.VariableNames(63:63+214)';
snp_stats_int.snp_names=clsa_joint.Properties.VariableNames(63:63+214)';
thresh=0.01
snp_stats_int.snp_names(snp_stats_main.pValue<thresh)
snp_stats_int.snp_names(snp_stats_int.pValue<thresh)
snp_stats_int.snp_names(snp_stats_int_DEP.pValue<thresh)
