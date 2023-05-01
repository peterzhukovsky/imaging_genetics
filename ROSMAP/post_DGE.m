snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_CT1_META\eqtl.txt'));

snps(~strcmp(snps.uniqID, '2:27225144:A:C'),:)=[];
snps(~strcmp(snps.uniqID, '8:27098436:C:T'),:)=[];
%% global
clear
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_CT1_META\genes_eqtl.txt'));
snps(snps.eqtlMapSNPs==0,:)=[];
twas=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\TWAS_b_caudalmiddlefrontal_t_robust_normalvoom.csv');
eqtl= readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\right_oriented\FUMA_CT1_META\eqtl.txt');
 
for i=1:length(snps.HUGO)
    try
   ix=strcmp(snps.symbol{i}, twas.hugo); 
   eqtl_global(i,:)=twas(ix,:);
    end
end

for i=1:length(snps.HUGO)
    ss=eqtl(strcmp(eqtl.symbol, snps.symbol(i)) ,:);    %& strcmp(eqtl.tissue, 'Brain_Cortex')
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.8
        sign_fuma(i,1)=0;
    elseif s<0.2
        sign_fuma(i,1)=1;
    else sign_fuma(i,1)=NaN;
    end
    try; t=twas(strcmp(twas.hugo, snps.HUGO(i)),:); sign_twas(i,1)=double(t.t>0);
    catch;sign_twas(i,1)=NaN;end
end

[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)
    
%bonferroni thr=0.05/17431=2.8685e-06
logp=-log10(twas.P_Value);
figure; scatter(twas.logFC, logp, 15, 'filled','k'); hold on;
scatter(twas.logFC(logp>1.3), logp(logp>1.3), 15, [0.4 0.4 0.4],'filled'); hold on;
scatter(twas.logFC(twas.adj_P_Val<0.05), logp(twas.adj_P_Val<0.05), 15, 'filled','b'); hold on;

%scatter(twas.logFC(abs(twas.logFC)<0.5), logp(abs(twas.logFC)<0.5), 15, 'filled','k'); hold on;

ix=strcmp(twas.hugo, 'KANSL1');scatter(twas.logFC(ix), logp(ix), 25, 'filled','r'); hold on;
ix=strcmp(twas.hugo, 'ARL17A');scatter(twas.logFC(ix), logp(ix), 25, 'filled','r'); hold on;
ix=strcmp(twas.hugo, 'ARL17B');scatter(twas.logFC(ix), logp(ix), 25, 'filled','r'); hold on;
ix=strcmp(twas.hugo, 'LRRC37A');scatter(twas.logFC(ix), logp(ix), 25, 'filled','r'); hold on;
ix=strcmp(twas.hugo, 'LRRC37A2');scatter(twas.logFC(ix), logp(ix), 25, 'filled','r'); hold on;

ix=strcmp(twas.hugo, 'STMN4');scatter(twas.logFC(ix), logp(ix), 25, [1 0.62 0], 'filled'); hold on;
ix=strcmp(twas.hugo, 'DPYSL5');scatter(twas.logFC(ix), logp(ix), 25, [1 0.62 0], 'filled'); hold on;
ix=strcmp(twas.hugo, 'ARHGAP27');scatter(twas.logFC(ix), logp(ix), 25, [1 0.62 0], 'filled'); hold on;




%% caudal middle frontal
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_CT4_META\genes_eqtl.txt'));
snps(snps.eqtlMapSNPs==0,:)=[];
twas=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\TWAS_b_caudalmiddlefrontal_t_robust_normalvoom.csv');
eqtl= readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\right_oriented\FUMA_CT4_META\eqtl.txt');

clear eqtl_global
for i=1:length(snps.HUGO)
    try
   ix=strcmp(snps.symbol{i}, twas.hugo); 
   eqtl_global(i,:)=twas(ix,:);
    end
end


for i=1:length(snps.HUGO)
    ss=eqtl(strcmp(eqtl.symbol, snps.symbol(i)) ,:);    %& strcmp(eqtl.tissue, 'Brain_Cortex')
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.6
        sign_fuma(i,1)=1;
    elseif s<0.4
        sign_fuma(i,1)=0;
    else sign_fuma(i,1)=NaN;
    end
    try; t=twas(strcmp(twas.hugo, snps.HUGO(i)),:); sign_twas(i,1)=double(t.t>0);
    catch;sign_twas(i,1)=NaN;end
end
[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)

%% rostral Middle frontal
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_CT27_META\genes_eqtl.txt'));
snps(snps.eqtlMapSNPs==0,:)=[];
twas=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\TWAS_b_rostralmiddlefrontal_t_robust_normalvoom.csv');
eqtl= readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\right_oriented\FUMA_CT27_META\eqtl.txt');

clear eqtl_global
for i=1:length(snps.HUGO)
    try
   ix=strcmp(snps.symbol{i}, twas.hugo); 
   eqtl_global(i,:)=twas(ix,:);
    end
end


for i=1:length(snps.HUGO)
    ss=eqtl(strcmp(eqtl.symbol, snps.symbol(i)) ,:);    %& strcmp(eqtl.tissue, 'Brain_Cortex')
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.8
        sign_fuma(i,1)=1;
    elseif s<0.2
        sign_fuma(i,1)=0;
    else sign_fuma(i,1)=NaN;
    end
    try; t=twas(strcmp(twas.hugo, snps.HUGO(i)),:); sign_twas(i,1)=double(t.t>0);
    catch;sign_twas(i,1)=NaN;end
end
[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)

%% WMH
snps=readtable(strcat('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_WMH_META\genes.txt'));
snps(snps.eqtlMapSNPs==0,:)=[];
twas=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\WMH_TWAS_withextracolumns_t_robust_normalvoom.csv');
eqtl= readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\FUMA_WMH_META\eqtl.txt');

clear eqtl_global
for i=1:length(snps.HUGO)
    try
   ix=strcmp(snps.symbol{i}, twas.hugo); 
   eqtl_global(i,:)=twas(ix,:);
    end
end

for i=1:length(snps.HUGO)
    ss=eqtl(strcmp(eqtl.symbol, snps.symbol(i)) ,:);    %& strcmp(eqtl.tissue, 'Brain_Cortex')
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.8
        sign_fuma(i,1)=1;
    elseif s<0.2
        sign_fuma(i,1)=0;
    else sign_fuma(i,1)=NaN;
    end
    try; t=twas(strcmp(twas.hugo, snps.HUGO(i)),:); sign_twas(i,1)=double(t.t>0);
    catch;sign_twas(i,1)=NaN;end
end
[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparing direction of eqtl vs twas

clear
eq = readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\right_oriented\FUMA_CT1_META\eqtl.txt');
tt = readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\TWAS_b_caudalmiddlefrontal_t_robust_normalvoom.csv');
allsign=unique(eq.symbol);%allsign=table(allsign,'VariableNames', {'symbol'});
for i=1:length(allsign)
    ss=eq(strcmp(eq.symbol, allsign(i)),:);    
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.8
        sign_fuma(i)=0;
    elseif s<0.2
        sign_fuma(i)=1;
    else sign_fuma(i)=NaN;
    end
    try; t=tt(strcmp(tt.hugo, allsign(i)),:); sign_twas(i)=double(t.t>0);
    catch;sign_twas(i)=NaN;end
end
%imagesc([sign_fuma; sign_twas])

[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)

%%

clear
eq = readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\right_oriented\FUMA_CT4_META\eqtl.txt');
tt = readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\TWAS_b_caudalmiddlefrontal_t_robust_normalvoom.csv');
allsign=unique(eq.symbol);%allsign=table(allsign,'VariableNames', {'symbol'});
for i=1:length(allsign)
    ss=eq(strcmp(eq.symbol, allsign(i)),:);    
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.8
        sign_fuma(i)=1;
    elseif s<0.2
        sign_fuma(i)=0;
    else sign_fuma(i)=NaN;
    end
    try; t=tt(strcmp(tt.hugo, allsign(i)),:); sign_twas(i,1)=double(t.t>0);
    catch;sign_twas(i,1)=NaN;end
end
%imagesc([sign_fuma; sign_twas])

[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)
%%

clear
eq = readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\right_oriented\FUMA_CT27_META\eqtl.txt');
tt = readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\ROSMAP_DGE\TWAS_b_rostralmiddlefrontal_t_robust_normalvoom.csv');
allsign=unique(eq.symbol);%allsign=table(allsign,'VariableNames', {'symbol'});
for i=1:length(allsign)
    ss=eq(strcmp(eq.symbol, allsign(i)),:);    
    s=sum(strcmp(ss.alignedDirection, '+'))/length(ss.alignedDirection);
    if s > 0.8
        sign_fuma(i)=1;
    elseif s<0.2
        sign_fuma(i)=0;
    else sign_fuma(i)=NaN;
    end
    try; t=tt(strcmp(tt.hugo, allsign(i)),:); sign_twas(i,1)=double(t.t>0);
    catch;sign_twas(i,1)=NaN;end
end
%imagesc([sign_fuma; sign_twas])

[tbl, chi, p]=crosstab(sign_fuma, sign_twas)
[h,p,stats]=fishertest(tbl)

sum(sign_fuma==0 & sign_twas==0)
sum(sign_fuma==1 & sign_twas==1)

