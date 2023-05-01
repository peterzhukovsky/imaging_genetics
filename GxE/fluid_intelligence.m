clear
cd D:\Canada_2020\UK_biobank\reports
minimal=readtable('minimal_ukb44677.csv');
FI=readtable('D:\Canada_2020\UK_biobank\reports\ordered\FI.csv');
FI(isnan(FI.FluidIntelligenceScore20016_0_0),:)=[];FI.eid=FI.f_eid;FI.FID=FI.f_eid;
%horzcat(FI.f_eid(~isnan(FI.FluidIntelligenceScore20016_0_0)),FI.f_eid(~isnan(FI.FluidIntelligenceScore20016_0_0)));


load('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\snpextract\fullsample_for_Gf\CT_SNPs_FI.mat');
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\snpextract\fullsample_for_Gf\ancestries.mat');
age=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\snpextract\fullsample_for_Gf\age_all.csv');age.Var1=[];
FI_table=innerjoin(FI, alldta);
FI_table=innerjoin(FI_table, minimal);
FI_table=innerjoin(FI_table, ancestries);
FI_table=innerjoin(FI_table, age);
cd D:\Canada_2020\UK_biobank\reports\AD\genetics\FI

%age, x21022; sex:  x31

modelspec=strcat('FluidIntelligenceScore20016_0_0~ x31_0_0 +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + ', FI_table.Properties.VariableNames{i+9})
(FI_table.Properties.VariableNames{i+215})
for snp=1:242
    %Outcome=cognitive_ordered.x20197_2_0; APOE=snps(:, snp);APOE=round(APOE{:,1});APOE(CT_SNPs.FID==0)=NaN;
    modelspec=strcat('FluidIntelligenceScore20016_0_0~ x31_0_0*x21022_0_0+x21022_0_0^2 +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + ', FI_table.Properties.VariableNames{snp+9})
    ix= ~strcmp(FI_table.pop,'EUR'); %  isnan(whitebritish_ordered); % %T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
    mdl = fitlm(FI_table(~ix,:),modelspec); 
    SNP_Gf_P(snp,1)=mdl.Coefficients.pValue(2);
    SNP_Gf_beta(snp,1)=mdl.Coefficients.Estimate(2);
    SNP_Gf_SE(snp,1)=mdl.Coefficients.SE(2);
end;
snps=FI_table.Properties.VariableNames(10:251)';
SNP_Gf_Pfdr=mafdr(SNP_Gf_P, 'BHFDR', true);sum(SNP_Gf_Pfdr<0.1)
summary_stats=table(snps, SNP_Gf_beta, SNP_Gf_SE, SNP_Gf_P)

tmp={'rs10052710_T';'rs10164339_G';'rs10253861_G';'rs10276148_G';'rs10498754_T';'rs10512248_T';'rs10753232_C';'rs10782438_T';'rs10821910_T';'rs10857115_G';'rs10912894_T';'rs10953403_T';'rs11063592_T';'rs11126806_C';'rs11128657_C';'rs11129967_G';'rs11197843_G';'rs112886906_A';'rs114898153_C';'rs11612673_T';'rs11692435_G';'rs11695197_G';'rs11695548_C';'rs11733915_T';'rs117421612_G';'rs1180379_G';'rs1182182_C';'rs12130849_C';'rs12187568_C';'rs12436192_G';'rs12500277_T';'rs12534154_T';'rs12603072_A';'rs1262815_C';'rs12668167_A';'rs12712508_A';'rs12886_T';'rs12904134_G';'rs12921170_G';'rs12938775_G';'rs12946564_G';'rs12979428_G';'rs12992813_T';'rs13059691_A';'rs13079368_C';'rs13082214_C';'rs13107325_C';'rs13110077_C';'rs13111112_G';'rs13128427_T';'rs13135092_A';'rs13156484_G';'rs13233398_C';'rs13424020_A';'rs13431289_G';'rs142256392_C';'rs1440802_T';'rs1441102_A';'rs1452628_A';'rs145567534_T';'rs1464047_C';'rs1466338_G';'rs147787064_C';'rs148372764_T';'rs148633354_G';'rs1511305_C';'rs1562330_C';'rs170239_C';'rs17055142_C';'rs17443495_G';'rs17496249_A';'rs17497343_T';'rs17576323_T';'rs17616633_T';'rs17785382_A';'rs1786345_A';'rs1861435_T';'rs1862664_A';'rs186347_G';'rs189807142_G';'rs190845364_T';'rs194834_G';'rs1948948_C';'rs1962848_A';'rs2002058_C';'rs2018980_A';'rs2033939_G';'rs2123161_A';'rs2241722_C';'rs224707_C';'rs2275186_C';'rs2287523_T';'rs2319627_T';'rs236607_C';'rs2415142_T';'rs245100_A';'rs2532261_A';'rs2532351_C';'rs2643222_A';'rs2677100_C';'rs2732702_T';'rs27437_A';'rs275350_C';'rs2819865_A';'rs28364628_C';'rs2844248_A';'rs28758523_C';'rs2916068_A';'rs3124591_C';'rs3200031_C';'rs333116_A';'rs341516_T';'rs34732613_T';'rs34876360_G';'rs34976331_A';'rs35021943_A';'rs35110958_C';'rs35831787_T';'rs360709_C';'rs3733574_T';'rs375326092_T';'rs3758098_C';'rs3786824_G';'rs3801650_A';'rs3816046_C';'rs3862469_C';'rs4606711_T';'rs4630220_G';'rs4630591_C';'rs4665943_C';'rs4684160_G';'rs4685022_A';'rs4773185_A';'rs4836008_T';'rs4841222_C';'rs4872084_A';'rs492623_T';'rs4937515_G';'rs4950101_T';'rs4962692_G';'rs4968317_G';'rs532_A';'rs534050929_C';'rs55823223_G';'rs56023709_A';'rs58158339_G';'rs58531798_T';'rs599823_T';'rs60371134_C';'rs61753077_G';'rs617577_T';'rs61784835_C';'rs62049363_A';'rs62051313_G';'rs62052489_T';'rs62053943_C';'rs6492273_T';'rs6492276_A';'rs6492277_G';'rs6698781_A';'rs6715800_C';'rs682571_G';'rs6877300_C';'rs6946266_G';'rs6984193_C';'rs7003744_T';'rs7094871_C';'rs71399925_G';'rs7142241_G';'rs7182018_A';'rs7218235_G';'rs729638_T';'rs73217710_G';'rs736866_C';'rs7419753_G';'rs74476491_C';'rs75112081_T';'rs7567226_C';'rs757425_A';'rs7575796_A';'rs7589627_C';'rs7592902_A';'rs7596872_C';'rs76049219_A';'rs76341705_G';'rs76928645_C';'rs77126132_G';'rs7733216_T';'rs77690628_T';'rs78580207_A';'rs7861718_A';'rs7875050_G';'rs78783493_C';'rs79104251_A';'rs79274348_C';'rs79362376_C';'rs7936928_T';'rs79852959_A';'rs8073909_T';'rs8081157_A';'rs8101493_A';'rs830179_G';'rs9322194_C';'rs9466_T';'rs9515212_A';'rs9521141_T';'rs9521760_A';'rs9559797_G';'rs9797946_T';'rs983960_C';'rs9893575_C';'rs9904766_T';'rs9926320_G';'rs9933149_T';'rs9967368_C'};
tmp=table(tmp, 'VariableNames', {'snps'});
tmp=innerjoin(tmp, summary_stats);
tmp.pfdr=mafdr(tmp.SNP_Gf_P, 'BHFDR', true);

cd D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\snpextract\fullsample_for_Gf\
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
% save alldta
ancestries=readtable('D:\Canada_2020\UK_biobank\reports\AD\all_pops_non_eur_pruned_within_pop_pc_covs.txt');
bridge=readtable('D:\Canada_2020\UK_biobank\reports\AD\ukb61530bridge31063.txt');bridge=sortrows(bridge, 2);
for i=265540:length(ancestries.s)
try; ancestries.eid(i)=bridge.Var1( bridge.Var2==ancestries.s(i) );end
end
a=ancestries;
ancestries(ancestries.eid==0,:)=[];
% save ancestries
