%% pls for chromosome 17
clear
%% load up data used in UKB_GxE.m or UKB_GxE_cross_validation.m
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\reports\crossval_inprogress.mat');
chr17snps={'rs117421612_G';'rs12603072_A';'rs12938775_G';'rs12946564_G';'rs2532261_A';'rs2532351_C';'rs2732702_T';'rs28364628_C';'rs333116_A';'rs35110958_C';'rs375326092_T';'rs4630591_C';'rs4968317_G';'rs55684829_A';'rs60371134_C';'rs62053943_C';'rs62073094_C';'rs62073102_C';'rs7206949_C';'rs7218235_G';'rs7219015_C';'rs736866_C';'rs8073909_T';'rs8081157_A';'rs9893575_C';'rs9904766_T'};
chr17snps={'rs117421612_G';'rs2532261_A';'rs2532351_C';'rs2732702_T';'rs4630591_C';'rs55684829_A';'rs62053943_C';'rs62073094_C';'rs62073102_C';'rs7206949_C';'rs736866_C';'rs9904766_T'};
chr17thicknesses={'cMFG';'cuneus';'entorhinal';'fPole';'fusiform';'iparietal';'itemporal';'isthmuscingulate';'loccipital';'lorbitofrontal';'lingual';'morbitofrontal';'mtemporal';'pocentral';'phippocampal';'popercularis';'porbitalis';'ptriangularis';'pcentral';'precentral';'precuneus';'rACC';'rMFG';'SupF';'SupPar';'SupT';'SupMarg'};%'GlobalMeanMean thickness';
chr17thicknesses={'SupF';'itemporal';'mtemporal';'fusiform';'loccipital';'rMFG';'cMFG';'iparietal';'SupMarg';'phippocampal';'porbitalis';'morbitofrontal';'rACC'};

%intersnps 
for i=1:length(chr17snps); ix=strcmp(chr17snps(i), snps.Properties.VariableNames); X(:,i)=snps{:, ix==1}; end
tmp=(thicknessFS(:,1:34)+thicknessFS(:,35:68))/2;
for i=1:length(chr17thicknesses); ix=strcmp(chr17thicknesses(i), regions_short)'; Y(:,i)=tmp(:, ix==1); end

X(CT_SNPs.FID==0,:)=NaN;
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %  isnan(whitebritish_ordered); %
X(ix,:)=[];Y(ix,:)=[];
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV,PCAs(:,1),PCAs(:,2),PCAs(:,3),PCAs(:,4),PCAs(:,5),PCAs(:,6),PCAs(:,7),PCAs(:,8),PCAs(:,9),PCAs(:,10),Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 
for i=1:12 [b,bint,r] = regress(X(:,i),[T.age,T.age.*T.sex, (T.age).^2, T.sex, T.PC1,T.PC2,T.PC3,T.PC4,T.PC5,T.PC6,T.PC7,T.PC8,T.PC9,T.PC10]); x(:,i)=r;end
for i=1:13 [b,bint,r] = regress(Y(:,i),[T.age,T.age.*T.sex, (T.age).^2, T.sex, T.TIV]); y(:,i)=r;end
X=x; Y=y;
%%% remove nans
naninx=sum(isnan(X)')'>0| sum(isnan(Y)')'>0; Y=Y(naninx==0,:);  X=X(naninx==0,:);  
Y=zscore(Y); X=zscore(X);x=X;y=Y;
ncomp=12
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%% permutation testing
permutations=5000;  
allobservations=Y;
for ncomp=12
   
    parfor n = 1:permutations
    % selecting either next combination, or random permutation
    permutation_index = randperm(length(allobservations));
    % creating random sample based on permutation index
    randomSample = allobservations(permutation_index,:);
    % running the PLS for this permutation
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,randomSample,ncomp);
    Rsq(n) = sum(PCTVAR(2,:));
    Rsq1(n) = sum(PCTVAR(1,:));  %if ncomp==4; c_perm_pls2(n,:)=corr(XS(:,2), Y);c_perm_pls1(n,:)=corr(XS(:,1), Y);c_perm_pls3(n,:)=corr(XS(:,3), Y);end
    end
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
    p(ncomp)=sum(sum(PCTVAR(2,:))<Rsq')/permutations
    p_1(ncomp)=sum(sum(PCTVAR(1,:))<Rsq1')/permutations
end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
[c,pval]=corr(XS, Y); %c(1,:)=c(1,:)*-1; c(2,:)=c(2,:)*-1;  c(3,:)=c(3,:)*-1; % reshape(mafdr(reshape(pval, [1 36]), 'BHFDR', 'True'), [4 9])
pval=reshape(pval, [1 12*13]);pval=mafdr(pval, 'BHFDR', 'true'); pval=reshape(pval, [12 13]);
figure;imagesc(abs(c)'); colormap hot; colorbar; set(gca, 'YTick', 1:27, 'YTickLabel', chr17thicknesses); set(gca, 'XTick', 1:12, 'XTickLabel', {'LV1','LV2','LV3','LV4','LV5','LV6','LV7','LV8','LV9','LV10','LV11','LV12'}, 'XTickLabelRotation',0); 
figure; imagesc(abs(XL));colormap hot; colorbar;set(gca, 'YTick', 1:26, 'YTickLabel', chr17snps); set(gca, 'XTick', 1:12, 'XTickLabel', {'LV1','LV2','LV3','LV4','LV5','LV6','LV7','LV8','LV9','LV10','LV11','LV12'}, 'XTickLabelRotation',0); 
figure; histogram(Rsq); xlim([0 0.0029]);

y_pred = [ones(size(X,1),1) X]*BETA; corr(Y, y_pred)
roi=1; figure; scatter(Y(:,roi), y_pred(:,roi),5, 'filled', 'MarkerFaceColor',[0.9 0.9 0.9]);
[ p, yhat, ci ] = polypredci( Y(:,roi), y_pred(:,roi), 1, 0.99); hold on; tmp=[Y(:,roi),yhat+ci, yhat-ci, yhat]; tmp=sortrows(tmp); plot(tmp(:,1), tmp(:,2), 'k'); plot(tmp(:,1),tmp(:,3), 'k'); plot(tmp(:,1),tmp(:,4));
corr(y_pred(:,roi), Y(:,roi))
