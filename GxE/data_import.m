cd D:\Canada_2020\UK_biobank\reports\AD
clear;

%% tabular data
cd D:\Canada_2020\UK_biobank\reports\ordered
%selfreport_ordered=readtable('selfreport_ordered.csv');
APOE_ordered=readtable('APOE_ordered.csv');BDNF_ordered=readtable('BDNF_ordered.csv');ADprs_ordered=readtable('ADprs_ordered.csv');
AD_ordered=readtable('AD_ordered.csv');
healthy_ordered=readtable('healthy_ordered.csv');healthydate_ordered=readtable('healthydate_ordered.csv');
minimal_ordered=readtable('minimal_ordered.csv');minimal_add_ordered=readtable('minimal_add_ordered.csv');minimal_add_ordered.Var1=[];minimal_add_ordered.eid=[];minimal_add_ordered.x20433_0_0=[];minimal_add_ordered.x20434_0_0=[];minimal_add_ordered.x20435_0_0=[];minimal_add_ordered.x20436_0_0=[];minimal_add_ordered.x20437_0_0=[];minimal_add_ordered.x20438_0_0=[];minimal_add_ordered.x20439_0_0=[];minimal_add_ordered.x20440_0_0=[];minimal_add_ordered.x20441_0_0=[];minimal_add_ordered.x20442_0_0=[];minimal_add_ordered.x20445_0_0=[];minimal_add_ordered.x20446_0_0=[];minimal_add_ordered.x20447_0_0=[];minimal_add_ordered.x20448_0_0=[];minimal_add_ordered.x20449_0_0=[];minimal_add_ordered.x20450_0_0=[];
minimal_ordered=horzcat(minimal_ordered, minimal_add_ordered);
ct_ordered=dlmread('ct_ordered.csv');
%ms_ordered=dlmread('ms_ordered.csv');
%gc_ordered_dem=dlmread('gc_ordered_dem.csv');
%gc_ordered=dlmread('gc_ordered.csv');
%yeormat_ordered=dlmread('yeormat_ordered.csv');
cognitive_ordered=readtable('cognitive_ordered.csv');
thicknessFS=readtable('thicknessFS_ordered.csv');thicknessFS=thicknessFS{:,:}; thicknessFS(:,1:2)=[];
physical_ordered=readtable('physical_ordered.csv');
%MDD_prs_ordered=dlmread('MDD_prs_ordered.csv');y=MDD_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); MDD_prs_ordered=r;
%ANX_prs_ordered=dlmread('ANX_prs_ordered.csv');y=ANX_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); ANX_prs_ordered=r;
%PTSD_prs_ordered=dlmread('PTSD_prs_ordered.csv');y=PTSD_prs_ordered;X=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10]; [b,bint,r] = regress(y,X); PTSD_prs_ordered=r;
load('whitebritish_ordered.mat');
cd D:\Canada_2020\UK_biobank\data
%ica_partial=dlmread('ica_d25_full.csv');
ica_partial=dlmread('ica_d25_par.csv');

cd D:\Canada_2020\UK_biobank\reports
datanames_rs=readtable('subjects_rs.txt');
datanames_fs=readtable('subjects_fs.txt');
datanames_all=unique(vertcat(datanames_fs.Var1, datanames_rs.Var1));
regions=readtable('thicknessFS_names.xlsx');regions=regions.names_short;regions(37)=[];regions(1)=[];%regions(72)=[];regions(69)=[];regions(37)=[];regions(36)=[];regions(33)=[]; regions(1)=[];

minimal_ordered.x2050_2_0=minimal_ordered.x2050_2_0-1; minimal_ordered.x2060_2_0=minimal_ordered.x2060_2_0-1;
minimal_ordered.x2050_2_0(minimal_ordered.x2050_2_0<0)=NaN;minimal_ordered.x2060_2_0(minimal_ordered.x2060_2_0<0)=NaN;
minimal_ordered.x2070_2_0(minimal_ordered.x2070_2_0<0)=NaN; minimal_ordered.x2080_2_0(minimal_ordered.x2080_2_0<0)=NaN; 
llabel_names=readtable('D:\Canada_2020\UK_biobank\data\ROI_names\lh.rsfc_HCP.txt');llabel_names=llabel_names.Var5;
rlabel_names=readtable('D:\Canada_2020\UK_biobank\data\ROI_names\rh.rsfc_HCP.txt');rlabel_names=rlabel_names.Var5;
label_names_all=vertcat(llabel_names,rlabel_names);  label_names_HOA={'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Postcentral Gyrus';'Precentral/postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Middle Frontal Gyrus';'Occipital Cortex';'Posterior Cingulate Cortex';'Precuneous ';'Occipital Cortex';'Lateral Occipital Cortex ';'Temporal/Occipital Fisiform Cortex';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Heschl''s Gyrus';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Superior Frontal Gyrus';'Precuneous ';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Precuneous ';'Precuneous ';'Precuneous ';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Posterior Cingulate Cortex';'Superior Parietal Lobule';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Superior Parietal Lobule';'Supplementary Motor Cortex';'Superior Frontal Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Superior Parietal Lobule ';'Superior Parietal Lobule ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Postcentral Gyrus';'Postcentral Gyrus';'Middle Postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (ventral, anterior)';'Paracingulate Gyrus';'Medial Superior Frontal Gyrus';'Medial Superior Frontal Gyrus';'Paracingulate Gyrus (anterior)';'Frontal Orbital Cortex';'Middle Frontal Gyrus';'Superior/Middle Frontal Gyrus';'Superior Frontal Gyrus/Frontal Pole';'Frontal Pole';'Frontal Pole (anterior)';'Frontal Pole (anterior)';'Middle Frontal Gyrus/ Inferior Frontal Gyrus';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Frontal Orbital Cortex';'Frontal Orbital Cortex/Frontal Pole (lateral)';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus/Inferior Frontal Gyrus';'Inferior Frontal Gyrus/Middle Frontal Gyrus';'Inferior Frontal Gyrus, pars triangularis ';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Frontal Pole';'Frontal Pole (dorsal anterior)';'Frontal Medial Cortex ';'Frontal Pole (ventral)';'Frontal Pole (ventral)';'Frontal Orbital Cortex/Frontal Pole (medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Superior Parietal Lobule';'Precentral Gyrus (medial)';'Middle Frontal Gyrus ';'Middle Frontal Gyrus ';'Precentral Gyrus (lateral ventral)';'Central Opercular Cortex';'Parietal Opercular Cortex';'Insular/Parietal Opercular Cortex';'Heschl''s Gyrus';'Parietal Opercular Cortex';'Parietal Opercular Cortex';'Insular Cortex (posterior)';'Insular Cortex (ventral)';'Frontal Operculum Cortex ';'Insular Cortex (dorsal)';'Insular Cortex (ventral)';'Insular Cortex (anterior)/Frontal Orbital Cortex';'Insular Cortex (anterior)';'Central Opercular Cortex ';'Central Opercular Cortex ';'Insular Cortex (dorsal)';'Supramarginal Gyrus (anterior)';'Superior Parietal Lobule';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Posterior Cingulate Cortex (ventral)';'Parahippocampal Gyrus';'Temporal Pole (aTL)';'Planum Temporale (STG)';'Superior Temporal Gyrus (anterior)';'Parahippocampal/Lingual Gyrus';'Temporal Fusiform Cortex';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus (posterior)';'Middle Temporal Gyrus (posterior)';'Temporal Pole';'Middle Temporal Gyrus (anterior)';'Inferior Temporal Gyrus, temporooccipital ';'Inferior Temporal Gyrus, posterior division ';'Temporal Fusiform Cortex, posterior division ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, superior division ';'Precuneous';'Lateral Occipital Cortex, superior division ';'Supramarginal Gyrus, posterior division ';'Lateral Occipital Cortex, superior division ';'Lateral Occipital Cortex, superior division ';'Postcentral Gyrus (ventral)';'Supramarginal Gyrus, anterior division ';'Angular Gyrus/Supramarginal Gyrus';'Angular Gyrus/Lateral Occipital Cortex';'Angular Gyrus/Lateral Occipital Cortex';'Lateral Occipital Cortex';'Lingual Gyrus';'Temporal Occipital Fusiform Cortex';'Temporal Fusiform Cortex, posterior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lingual Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Temporal Occipital Fusiform Cortex ';'Subcallosal Cortex (Subgenual ACC)';'Paracingulate Gyrus (anterior, ventral)';'Subcallosal Cortex (Subgenual ACC)';'Insular Cortex';'Insular Cortex/Parietal Operculum';'Frontal Operculum Cortex ';'Frontal Pole (medial)';'Frontal Pole (lateral)';'Inferior Temporal Gyrus, anterior division ';'Heschl''s Gyrus ';'Planum Temporale (STG)';'Planum Temporale (STG)';'Middle Temporal Gyrus, posterior division';'Middle Temporal Gyrus, posterior division ';'Planum Polare';'Anterior Cingulate';'Anterior Cingulate';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Postcentral Gyrus';'Precentral/postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Middle Frontal Gyrus';'Occipital Cortex';'Posterior Cingulate Cortex';'Precuneous ';'Occipital Cortex';'Lateral Occipital Cortex ';'Temporal/Occipital Fisiform Cortex';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Heschl''s Gyrus';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Superior Frontal Gyrus';'Precuneous ';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Precuneous ';'Precuneous ';'Precuneous ';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Posterior Cingulate Cortex';'Superior Parietal Lobule';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Superior Parietal Lobule';'Supplementary Motor Cortex';'Superior Frontal Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Superior Parietal Lobule ';'Superior Parietal Lobule ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Postcentral Gyrus';'Postcentral Gyrus';'Middle Postcentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Precentral Gyrus';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate Cortex (ventral, anterior)';'Paracingulate Gyrus';'Medial Superior Frontal Gyrus';'Medial Superior Frontal Gyrus';'Paracingulate Gyrus (anterior)';'Frontal Orbital Cortex';'Middle Frontal Gyrus';'Superior/Middle Frontal Gyrus';'Superior Frontal Gyrus/Frontal Pole';'Frontal Pole';'Frontal Pole (anterior)';'Frontal Pole (anterior)';'Middle Frontal Gyrus/ Inferior Frontal Gyrus';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Frontal Orbital Cortex';'Frontal Orbital Cortex/Frontal Pole (lateral)';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus/Inferior Frontal Gyrus';'Inferior Frontal Gyrus/Middle Frontal Gyrus';'Inferior Frontal Gyrus, pars triangularis ';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Middle Frontal Gyrus';'Frontal Pole';'Frontal Pole (dorsal anterior)';'Frontal Medial Cortex ';'Frontal Pole (ventral)';'Frontal Pole (ventral)';'Frontal Orbital Cortex/Frontal Pole (medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Frontal Orbital Cortex (Medial)';'Superior Parietal Lobule';'Precentral Gyrus (medial)';'Middle Frontal Gyrus ';'Middle Frontal Gyrus ';'Precentral Gyrus (lateral ventral)';'Central Opercular Cortex';'Parietal Opercular Cortex';'Insular/Parietal Opercular Cortex';'Heschl''s Gyrus';'Parietal Opercular Cortex';'Parietal Opercular Cortex';'Insular Cortex (posterior)';'Insular Cortex (ventral)';'Frontal Operculum Cortex ';'Insular Cortex (dorsal)';'Insular Cortex (ventral)';'Insular Cortex (anterior)/Frontal Orbital Cortex';'Insular Cortex (anterior)';'Central Opercular Cortex ';'Central Opercular Cortex ';'Insular Cortex (dorsal)';'Supramarginal Gyrus (anterior)';'Superior Parietal Lobule';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Parahippocampal Gyrus';'Posterior Cingulate Cortex (ventral)';'Parahippocampal Gyrus';'Temporal Pole (aTL)';'Planum Temporale (STG)';'Superior Temporal Gyrus (anterior)';'Parahippocampal/Lingual Gyrus';'Temporal Fusiform Cortex';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus (posterior)';'Middle Temporal Gyrus (posterior)';'Temporal Pole';'Middle Temporal Gyrus (anterior)';'Inferior Temporal Gyrus, temporooccipital ';'Inferior Temporal Gyrus, posterior division ';'Temporal Fusiform Cortex, posterior division ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, superior division ';'Precuneous';'Lateral Occipital Cortex, superior division ';'Supramarginal Gyrus, posterior division ';'Lateral Occipital Cortex, superior division ';'Lateral Occipital Cortex, superior division ';'Postcentral Gyrus (ventral)';'Supramarginal Gyrus, anterior division ';'Angular Gyrus/Supramarginal Gyrus';'Angular Gyrus/Lateral Occipital Cortex';'Angular Gyrus/Lateral Occipital Cortex';'Lateral Occipital Cortex';'Lingual Gyrus';'Temporal Occipital Fusiform Cortex';'Temporal Fusiform Cortex, posterior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, inferior division ';'Lingual Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Temporal Occipital Fusiform Cortex ';'Subcallosal Cortex (Subgenual ACC)';'Paracingulate Gyrus (anterior, ventral)';'Subcallosal Cortex (Subgenual ACC)';'Insular Cortex';'Insular Cortex/Parietal Operculum';'Frontal Operculum Cortex ';'Frontal Pole (medial)';'Frontal Pole (lateral)';'Inferior Temporal Gyrus, anterior division ';'Heschl''s Gyrus ';'Planum Temporale (STG)';'Planum Temporale (STG)';'Middle Temporal Gyrus, posterior division';'Middle Temporal Gyrus, posterior division ';'Planum Polare';'Anterior Cingulate';'Anterior Cingulate'};
ica2yeo7=readtable('ica2yeo7.csv');

minimal_ordered.x21000_0_0; tmp=repmat( {'Missing'} , length(datanames_all),1 );
tmp(minimal_ordered.x21000_0_0==1|minimal_ordered.x21000_0_0==1001|minimal_ordered.x21000_0_0==1002|minimal_ordered.x21000_0_0==1003)={'White'};
tmp(~(minimal_ordered.x21000_0_0==1|minimal_ordered.x21000_0_0==1001|minimal_ordered.x21000_0_0==1002|minimal_ordered.x21000_0_0==1003))={'Non-White'};
minimal_ordered.x21000_0_0=tmp;
tmp=repmat( {'Missing'} , length(datanames_all),1 ); 
tmp(minimal_ordered.x54_2_0==11025)={'Cheadle'}; tmp(minimal_ordered.x54_2_0==11026)={'Reading'}; tmp(minimal_ordered.x54_2_0==11027)={'Newcastle'}; 
minimal_ordered.x54_2_0=tmp; clear tmp;


%% %% %% %%%% %%%% %%%% %%%% %%
%% depression vs anxiety vs stress -  rsfmri ICA; setting up the clinical variables
%% %% %% %%%% %%%% %%%% %%%% %%
id_depressed_clin=( (minimal_ordered.x130895_0_0>=30 & minimal_ordered.x130895_0_0<=51) | (healthy_ordered.x130897_0_0 >=30 & healthy_ordered.x130897_0_0<=51) ); sum(id_depressed_clin>0)
id_anxious_clin=(healthy_ordered.x130907_0_0>=30 & healthy_ordered.x130907_0_0<=51); sum(id_anxious_clin>0)
id_stress_clin=(healthy_ordered.x130911_0_0>=30 & healthy_ordered.x130911_0_0<=51); sum(id_stress_clin)
id_stress=(id_stress_clin>0 & id_anxious_clin==0 & id_depressed_clin==0); 
id_anxious=(id_stress_clin==0 & id_anxious_clin>0 & id_depressed_clin==0); 
id_dep=(id_stress_clin==0 & id_anxious_clin==0 & id_depressed_clin>0); 
id_depanx=(id_stress_clin==0 & id_anxious_clin>0 & id_depressed_clin>0); 
id_healthy=sum(~isnan(healthy_ordered{:,8:145})')';
%for i=1:length(id_healthy);
%    if minimal_ordered.x21003_2_0(i)>60 & id_healthy(i)==0
%        id_healthy(i)=round(rand);
%end; end
sum(id_healthy==0)

clinical={'a'}; %clinical(age<60)={[]};
clinical(id_healthy==0)={'ahc'};
clinical(id_dep==1)={'dep'};
clinical(id_depanx==1)={'depanx'};
clinical(id_anxious==1)={'anx'};
clinical(id_stress==1)={'str'};clinical(40699)=clinical(40697);
age=minimal_ordered.x21003_2_0;sex=minimal_ordered.x31_0_0; 

ancestries_ordered=readtable( 'D:\Canada_2020\UK_biobank\reports\ordered\ancestries_ordered.csv');ancestries_ordered.pop(ancestries_ordered.eid==0)={'NaN'};      
               
               
               
               
               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
%histogram(APOE_ordered.rs429358_T); histogram(APOE_ordered.rs7412_C)
TIV=AD_ordered.x26521_2_0; %eTIV
AD_ordered.x25781_2_0(AD_ordered.x25781_2_0>15000)=NaN; histogram(AD_ordered.x25781_2_0); title('WMH');% WHM
AD_ordered.x22506_0_0=AD_ordered.x22506_0_0-110;AD_ordered.x22506_0_0(AD_ordered.x22506_0_0 <0)=NaN;%tobacco smoking
physical_ordered.x6150_0_0(physical_ordered.x6150_0_0<0)=0; %cardio
physical_ordered.x6150_0_0(physical_ordered.x6150_0_0>0)=1; %physical_ordered.x6150_0_0(physical_ordered.x6150_0_0>0 & physical_ordered.x6150_0_0<4)=1;
regions={'GlobalMeanMean thickness (left hemisphere)';'bankssts (left hemisphere)';'caudalanteriorcingulate (left hemisphere)';'caudalmiddlefrontal (left hemisphere)';'cuneus (left hemisphere)';'entorhinal (left hemisphere)';'fusiform (left hemisphere)';'inferiorparietal (left hemisphere)';'inferiortemporal (left hemisphere)';'isthmuscingulate (left hemisphere)';'lateraloccipital (left hemisphere)';'lateralorbitofrontal (left hemisphere)';'lingual (left hemisphere)';'medialorbitofrontal (left hemisphere)';'middletemporal (left hemisphere)';'parahippocampal (left hemisphere)';'paracentral (left hemisphere)';'parsopercularis (left hemisphere)';'parsorbitalis (left hemisphere)';'parstriangularis (left hemisphere)';'pericalcarine (left hemisphere)';'postcentral (left hemisphere)';'posteriorcingulate (left hemisphere)';'precentral (left hemisphere)';'precuneus (left hemisphere)';'rostralanteriorcingulate (left hemisphere)';'rostralmiddlefrontal (left hemisphere)';'superiorfrontal (left hemisphere)';'superiorparietal (left hemisphere)';'superiortemporal (left hemisphere)';'supramarginal (left hemisphere)';'frontalpole (left hemisphere)';'transversetemporal (left hemisphere)';'insula (left hemisphere)';'GlobalMeanMean thickness (right hemisphere)';'bankssts (right hemisphere)';'caudalanteriorcingulate (right hemisphere)';'caudalmiddlefrontal (right hemisphere)';'cuneus (right hemisphere)';'entorhinal (right hemisphere)';'fusiform (right hemisphere)';'inferiorparietal (right hemisphere)';'inferiortemporal (right hemisphere)';'isthmuscingulate (right hemisphere)';'lateraloccipital (right hemisphere)';'lateralorbitofrontal (right hemisphere)';'lingual (right hemisphere)';'medialorbitofrontal (right hemisphere)';'middletemporal (right hemisphere)';'parahippocampal (right hemisphere)';'paracentral (right hemisphere)';'parsopercularis (right hemisphere)';'parsorbitalis (right hemisphere)';'parstriangularis (right hemisphere)';'pericalcarine (right hemisphere)';'postcentral (right hemisphere)';'posteriorcingulate (right hemisphere)';'precentral (right hemisphere)';'precuneus (right hemisphere)';'rostralanteriorcingulate (right hemisphere)';'rostralmiddlefrontal (right hemisphere)';'superiorfrontal (right hemisphere)';'superiorparietal (right hemisphere)';'superiortemporal (right hemisphere)';'supramarginal (right hemisphere)';'frontalpole (right hemisphere)';'transversetemporal (right hemisphere)';'insula (right hemisphere)'};
physical_ordered.x21001_2_0(physical_ordered.x21001_2_0>45)=NaN; %BMI bigger 45

AD_ordered.x25019_2_0=clean(AD_ordered.x25019_2_0, 4); % Hippo left
AD_ordered.x25020_2_0=clean(AD_ordered.x25020_2_0, 4);% Hippo right
AD_ordered.x26793_2_0 =clean(AD_ordered.x26793_2_0 , 4);% Ento left
AD_ordered.x26894_2_0 =clean(AD_ordered.x26894_2_0 , 4);% Ento right
AD_ordered.x6138_0_0(AD_ordered.x6138_0_0<0)=NaN; %education remove not indicated

Hearing=repmat(NaN, 40699,1);Hearing(AD_ordered.x2247_0_0==1|AD_ordered.x2247_2_0==1|AD_ordered.x2257_0_0==1|AD_ordered.x2257_2_0==1)=1;
Hearing(AD_ordered.x2247_0_0==0|AD_ordered.x2247_2_0==0|AD_ordered.x2257_0_0==0|AD_ordered.x2257_2_0==0)=0;
HRT=repmat(NaN, 40699,1); HRT(AD_ordered.x2814_0_0==1|AD_ordered.x2814_2_0==1)=1; HRT(AD_ordered.x2814_0_0==0|AD_ordered.x2814_2_0==0)=0; HRT(sex==1)=NaN;
AD_ordered.x3581_0_0(AD_ordered.x3581_0_0<1)=NaN; %age at period
AD_ordered.x20414_0_0(AD_ordered.x20414_0_0<0)=NaN; %freq alc
AD_ordered.x20416_0_0(AD_ordered.x20416_0_0<1)=NaN; %freq alc 6 units

cognitive_ordered.x6350_2_0(cognitive_ordered.x6350_2_0> (nanmean(cognitive_ordered.x6350_2_0)+4*nanstd(cognitive_ordered.x6350_2_0)) | cognitive_ordered.x6350_2_0<100)=NaN; %Duration to complete alphanumeric path 
cognitive_ordered.x6348_2_0(cognitive_ordered.x6348_2_0<100 | cognitive_ordered.x6348_2_0> (nanmean(cognitive_ordered.x6350_2_0)+3*nanstd(cognitive_ordered.x6348_2_0)) )=NaN; %Duration to complete numeric/easy path 
cognitive_ordered.tmt_cor=(cognitive_ordered.x6350_2_0 + 5*cognitive_ordered.x6351_2_0);% - (cognitive_ordered.x6348_2_0 +5*cognitive_ordered.x6349_2_0) ;

% depression PHQ2
PHQ2=(minimal_ordered.x2050_2_0 + minimal_ordered.x2060_2_0)>=2 ; sum(PHQ2>0)
% prs for brain variables
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\hippovol\hippoprs_FUMA_indSNP.mat');
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\WMH\wmhprs_FUMA_genes.mat');
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\thickness\thickness_prs.mat');
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\MD\mdprs.mat');
load('D:\Canada_2020\UK_biobank\reports\AD\genetics\CT\1\ct1prs.mat');

PCAs=[healthydate_ordered.x22009_0_1,healthydate_ordered.x22009_0_2, healthydate_ordered.x22009_0_3, healthydate_ordered.x22009_0_4, healthydate_ordered.x22009_0_5,healthydate_ordered.x22009_0_6,healthydate_ordered.x22009_0_7,healthydate_ordered.x22009_0_8,healthydate_ordered.x22009_0_9,healthydate_ordered.x22009_0_10];
[b,bint,r] = regress(HippoPRS',PCAs); HippoPRS=r';HippoPRS=clean(HippoPRS, 5);
[b,bint,r] = regress(WMH_PRS',PCAs); WMH_PRS=r';WMH_PRS=clean(WMH_PRS, 5);
[b,bint,r] = regress(thickness_PRS',PCAs); thickness_PRS=r';thickness_PRS=clean(thickness_PRS, 5);
[b,bint,r] = regress(MD_PRS',PCAs); MD_PRS=r';MD_PRS=clean(MD_PRS, 6);
[b,bint,r] = regress(CT1_PRS',PCAs); CT1_PRS=r';CT1_PRS=clean(CT1_PRS, 6);

%load('D:\Canada_2020\UK_biobank\reports\AD\genetics\CT\15\ctprs.mat');
%[b,bint,r] = regress(CT_PRS',PCAs); CT_PRS=r';CT15_PRS=clean(CT_PRS, 6);
ancestries_ordered=readtable( 'D:\Canada_2020\UK_biobank\reports\ordered\ancestries_ordered.csv');ancestries_ordered.pop(ancestries_ordered.eid==0)={'NaN'};


               
