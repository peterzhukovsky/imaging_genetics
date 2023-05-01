
ct_rg=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\sumstats.txt');
h2=ct_rg.h2_obs(1:34); h2_se=ct_rg.h2_obs_se(1:34);
ct_rg_square=reshape(ct_rg.rg, [34 34]); %tmp=reshape(ct_rg.p, [34 34]); ct_rg_square(tmp>0.05/(34*34))=0;
ct_rg_square(1,:)=[];ct_rg_square(:,1)=[];%ct_rg_square(eye(33)==1)=NaN;
ct_p_square=reshape(ct_rg.p, [34 34]);ct_p_square(1,:)=[];ct_p_square(:,1)=[]; ct_rg_square(ct_p_square>5e-5)=0;
figure; imagesc(ct_rg_square);tmp=regions; tmp(1)=[]; set(gca, 'XTick', 1:33, 'XTickLabel', regions_short, 'XTickLabelRotation',90); set(gca, 'YTick', 1:33, 'YTickLabel', regions_short); 
modules=community_louvain(ct_rg_square, 0.95);
[X,Y,INDSORT] = grid_communities(modules)
imagesc(ct_rg_square(INDSORT,INDSORT));hold on;plot(X,Y,'r','linewidth',2);
set(gca, 'XTick', 1:33, 'XTickLabel', regions_short(INDSORT), 'XTickLabelRotation',90);
set(gca, 'YTick', 1:33, 'YTickLabel', regions_short(INDSORT)); colormap hot
figure;imagesc(h2(INDSORT)'); colormap hot
modules_all=M;
for i=1:499
M=community_louvain(ct_rg_square, 0.95);
modules_all=[modules_all,M];
u_m_all(i)=length(unique(M));
end
modules_all(:,sum(modules_all)<50)=[];
regions_short(modules==1)
regions_short(modules==2)
regions_short(modules==3)
figure;histogram(modules)

%% genetic correlations science - reshape the uncorrected genetic correlations to compare to our data
rg=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\reports\science.aay6690_tables_regionalrg.xlsx');
rg_square=reshape(rg.rG, [33 34]);rg_square=[rg_square; ones(1,34)]; %the columns arent sorted in the same way however


%% genetic correlations

rg=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\reports\ct_rg_comprehensive.txt');
rg=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\reports\wmh_rg_comprehensive.txt');
start=4;
data = rg.rg(start:length(rg.rg));
errhigh = 1.96*rg.se(start:length(rg.rg));
errlow  = 1.96*rg.se(start:length(rg.rg));
x = 1:length(data);

figure;
bar(x,data); hold on; er = errorbar(x,data,errlow,errhigh); er.Color = [0 0 0]; 
er.LineStyle = 'none';  hold off

%% regional correlations
ct_rg=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\reports\sumstats_regional.txt');
%ct_rg=ct_rg(2:7:233,:);%ct_rg=ct_rg(7:7:238,:);
ct_rg.regions=regions(1:34);
sortrows(ct_rg,3)
