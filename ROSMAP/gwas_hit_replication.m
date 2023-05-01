module load MATLAB/2020b
matlab -nodesktop -nodisplay

cp /external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_* ./

map_thick_l=readtable('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_l_map_share.csv','TreatAsMissing',{'.','NULL'});
map_thick_r=readtable('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_r_map_share.csv','TreatAsMissing',{'.','NULL'});
ros_thick_l=readtable('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_l_ros_share.csv','TreatAsMissing',{'.','NULL'});
ros_thick_r=readtable('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_r_ros_share.csv','TreatAsMissing',{'.','NULL'});
ros_etiv=readtable('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_subcortical_v6_ros_share.csv','TreatAsMissing',{'.','NULL'});
map_etiv=readtable('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_subcortical_v6_map_share.csv','TreatAsMissing',{'.','NULL'});

rosmap_thick_l=vertcat(ros_thick_l,map_thick_l);
rosmap_thick_r=vertcat(ros_thick_r,map_thick_r);
rosmap_thick_r.etiv=vertcat(ros_etiv.estimatedtotalintracranialvol , map_etiv.estimatedtotalintracranialvol);
rosmap_thick_l.etiv=vertcat(ros_etiv.estimatedtotalintracranialvol , map_etiv.estimatedtotalintracranialvol);
pheno=readtable('/external/rprshnas01/tigrlab/scratch/pzhukovsky/biobank/genetics/ROSMAP_rep/ROSMAP_master.csv','TreatAsMissing',{'.','NA'});

pcas=readtable('/KIMEL/tigrlab/scratch/pzhukovsky/biobank/genetics/ROSMAP_rep/rosmap_pca10.eigenvec.txt');

for roi=3:37; rosmap_thick_l{:,roi}=(rosmap_thick_l{:,roi}+rosmap_thick_r{:,roi})/2; end

rosmap_thick_l.Properties.VariableNames(37)={'GlobalMeanThickness'};
rosmap_thick_l{:,37}=mean(rosmap_thick_l{:,3:36}')';
rosmap_thick_r.Properties.VariableNames(37)={'GlobalMeanThickness'};
rosmap_thick_r{:,37}=mean(rosmap_thick_r{:,3:36}')';

rois={'caudalmiddlefrontal';'inferiorparietal';'lateraloccipital';'lingual';'postcentral';'precentral';'lingual';'lateralorbitofrontal';'middletemporal';'parstriangularis';'superiortemporal';'pericalcarine';'caudalanteriorcingulate';'lateralorbitofrontal';'medialorbitofrontal';'inferiorparietal';'postcentral';'precentral';'inferiorparietal';'lingual';'postcentral';'caudalmiddlefrontal';'inferiorparietal';'lateraloccipital';'lingual';'postcentral';'isthmuscingulate';'isthmuscingulate';'lingual';'rostralanteriorcingulate';'rostralanteriorcingulate';'rostralanteriorcingulate';'rostralanteriorcingulate';'lingual';'inferiortemporal';'superiortemporal';'parsorbitalis';'frontalpole';'superiortemporal';'superiortemporal';'superiortemporal';'superiortemporal';'insula';'paracentral';'postcentral';'parsorbitalis';'GlobalMeanThickness';'caudalmiddlefrontal';'inferiorparietal';'middletemporal';'superiorfrontal';'supramarginal';'parahippocampal';'lateralorbitofrontal';'precuneus'}

                          
                          

%################################# MAIN EFFECTS MAIN EFFECTS MAIN EFFECTS  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################%################################# LEFT SIDE  #######################

%################################# LEFT SIDE  #######################

dta=readtable('for_rosmap_rep.csv'); rois=dta.roi; snps_OI=dta.rsid;

rm_thick_r_pheno=innerjoin(rosmap_thick_l,pheno);
%writetable(rm_thick_r_pheno,rm_thick_r_pheno.csv); 
rm_thick_r_pheno.FID=[];
snps=readtable('ROSMAP_rep.txt');snps=innerjoin(snps, pcas, 'Keys',{'IID','IID'} );
rm_thick_r_pheno_snps=innerjoin(rm_thick_r_pheno,snps, 'Keys',{'IID','IID'} );
%writetable(rm_thick_r_pheno_snps,'rm_thick_r_pheno_snps.csv');
tmp=innerjoin(rm_thick_r_pheno_snps, rosmap_thick_l, 'Keys',{'projid','projid'});
roilabels=rm_thick_r_pheno_snps.Properties.VariableNames(3:36);

for j=1:length(rois) 
try
ix=contains(rm_thick_r_pheno_snps.Properties.VariableNames, snps_OI{j});APOE=rm_thick_r_pheno_snps{:,ix}; rsID(j)=rm_thick_r_pheno_snps.Properties.VariableNames(ix);
ix=contains(rm_thick_r_pheno_snps.Properties.VariableNames, rois(j)); Outcome=rm_thick_r_pheno_snps{:,ix};

cardio=rm_thick_r_pheno_snps.hypertension_bl; %	cardio=rm_thick_r_pheno_snps.hypertension_ever;
cardio=rm_thick_r_pheno_snps.stroke_bl==1 | rm_thick_r_pheno_snps.claudication_bl==1 | rm_thick_r_pheno_snps.cvda>0 | rm_thick_r_pheno_snps.hypertension_bl;
PC1=rm_thick_r_pheno_snps.PC1;PC2=rm_thick_r_pheno_snps.PC2;PC3=rm_thick_r_pheno_snps.PC3;PC4=rm_thick_r_pheno_snps.PC4;PC5=rm_thick_r_pheno_snps.PC5;PC6=rm_thick_r_pheno_snps.PC6;PC7=rm_thick_r_pheno_snps.PC7;PC8=rm_thick_r_pheno_snps.PC8;PC9=rm_thick_r_pheno_snps.PC9;PC10=rm_thick_r_pheno_snps.PC10;
etiv=rm_thick_r_pheno_snps.etiv;
T=table(rm_thick_r_pheno_snps.projid, APOE, cardio, rm_thick_r_pheno_snps.age_bl, rm_thick_r_pheno_snps.msex, rm_thick_r_pheno_snps.study, rm_thick_r_pheno_snps.visit, etiv, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10, Outcome, 'VariableNames', {'ID','APOE','cardio', 'age', 'sex','site','visit', 'etiv', 'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'});  

mdl = fitlme(T,'Outcome~APOE+sex+age+site+visit+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+etiv+(1|ID)'); %+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10

SNP_Cardio_CT_P(j)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Name,'APOE'));
SNP_Cardio_CT_beta(j)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'APOE'));
SNP_Cardio_CT_SE(j)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Name,'APOE'));
end
end
t=table(rsID',  SNP_Cardio_CT_beta', SNP_Cardio_CT_SE', SNP_Cardio_CT_P');writetable(t, 'sumstats.csv')


%################################# RIGHT SIDE  #######################
dta=readtable('for_rosmap_rep.csv'); rois=dta.roi; snps_OI=dta.rsid;

rm_thick_r_pheno=innerjoin(rosmap_thick_r,pheno);
%writetable(rm_thick_r_pheno,rm_thick_r_pheno.csv); 
rm_thick_r_pheno.FID=[];
snps=readtable('ROSMAP_rep.txt');snps=innerjoin(snps, pcas, 'Keys',{'IID','IID'} );
rm_thick_r_pheno_snps=innerjoin(rm_thick_r_pheno,snps, 'Keys',{'IID','IID'} );
%writetable(rm_thick_r_pheno_snps,'rm_thick_r_pheno_snps.csv');
tmp=innerjoin(rm_thick_r_pheno_snps, rosmap_thick_l, 'Keys',{'projid','projid'});
roilabels=rm_thick_r_pheno_snps.Properties.VariableNames(3:36);


for j=1:length(rois) 
try
ix=contains(rm_thick_r_pheno_snps.Properties.VariableNames, snps_OI{j});APOE=rm_thick_r_pheno_snps{:,ix}; rsID(j)=rm_thick_r_pheno_snps.Properties.VariableNames(ix);
ix=contains(rm_thick_r_pheno_snps.Properties.VariableNames, rois(j)); Outcome=rm_thick_r_pheno_snps{:,ix};

cardio=rm_thick_r_pheno_snps.hypertension_bl; %	cardio=rm_thick_r_pheno_snps.hypertension_ever;
cardio=rm_thick_r_pheno_snps.stroke_bl==1 | rm_thick_r_pheno_snps.claudication_bl==1 | rm_thick_r_pheno_snps.cvda>0 | rm_thick_r_pheno_snps.hypertension_bl;
PC1=rm_thick_r_pheno_snps.PC1;PC2=rm_thick_r_pheno_snps.PC2;PC3=rm_thick_r_pheno_snps.PC3;PC4=rm_thick_r_pheno_snps.PC4;PC5=rm_thick_r_pheno_snps.PC5;PC6=rm_thick_r_pheno_snps.PC6;PC7=rm_thick_r_pheno_snps.PC7;PC8=rm_thick_r_pheno_snps.PC8;PC9=rm_thick_r_pheno_snps.PC9;PC10=rm_thick_r_pheno_snps.PC10;
etiv=rm_thick_r_pheno_snps.etiv;
T=table(rm_thick_r_pheno_snps.projid, APOE, cardio, rm_thick_r_pheno_snps.age_bl, rm_thick_r_pheno_snps.msex, rm_thick_r_pheno_snps.study, rm_thick_r_pheno_snps.visit, etiv, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10, Outcome, 'VariableNames', {'ID','APOE','cardio', 'age', 'sex','site','visit', 'etiv', 'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Outcome'});  

mdl = fitlme(T,'Outcome~APOE+sex+age+visit+site+age*sex+age^2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+etiv+(1|ID)'); %+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10

SNP_Cardio_CT_P(j)=mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Name,'APOE'));
SNP_Cardio_CT_beta(j)=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,'APOE'));
SNP_Cardio_CT_SE(j)=mdl.Coefficients.SE(strcmp(mdl.Coefficients.Name,'APOE'));
end
end
t=table(rsID',  SNP_Cardio_CT_beta', SNP_Cardio_CT_SE', SNP_Cardio_CT_P');writetable(t, 'sumstats.csv')

%% demographics
ix=contains(rm_thick_r_pheno_snps.Properties.VariableNames, snps_OI{j});APOE=rm_thick_r_pheno_snps{:,ix}; rsID(j)=rm_thick_r_pheno_snps.Properties.VariableNames(ix);
ix=contains(rm_thick_r_pheno_snps.Properties.VariableNames, rois(j)); Outcome=rm_thick_r_pheno_snps{:,ix};

cardio=rm_thick_r_pheno_snps.hypertension_bl; %	cardio=rm_thick_r_pheno_snps.hypertension_ever;
cardio=rm_thick_r_pheno_snps.stroke_bl==1 | rm_thick_r_pheno_snps.claudication_bl==1 | rm_thick_r_pheno_snps.cvda>0 | rm_thick_r_pheno_snps.hypertension_bl;
PC1=rm_thick_r_pheno_snps.PC1;PC2=rm_thick_r_pheno_snps.PC2;PC3=rm_thick_r_pheno_snps.PC3;PC4=rm_thick_r_pheno_snps.PC4;PC5=rm_thick_r_pheno_snps.PC5;PC6=rm_thick_r_pheno_snps.PC6;PC7=rm_thick_r_pheno_snps.PC7;PC8=rm_thick_r_pheno_snps.PC8;PC9=rm_thick_r_pheno_snps.PC9;PC10=rm_thick_r_pheno_snps.PC10;
etiv=rm_thick_r_pheno_snps.etiv;
T=table(rm_thick_r_pheno_snps.projid, APOE, cardio, rm_thick_r_pheno_snps.age_bl, rm_thick_r_pheno_snps.msex, rm_thick_r_pheno_snps.study, rm_thick_r_pheno_snps.visit, etiv, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10, Outcome, 

rm_thick_r_pheno_snps.cardio=rm_thick_r_pheno_snps.hypertension_bl; %	cardio=rm_thick_r_pheno_snps.hypertension_ever;
rm_thick_r_pheno_snps.cardio=rm_thick_r_pheno_snps.stroke_bl==1 | rm_thick_r_pheno_snps.claudication_bl==1 | rm_thick_r_pheno_snps.cvda>0 | rm_thick_r_pheno_snps.hypertension_bl;

IDs=unique(rm_thick_r_pheno_snps.projid)
for i=1:length(IDs)
tmp=rm_thick_r_pheno_snps(rm_thick_r_pheno_snps.projid==IDs(i),:);
tmp=tmp(1,:);
Tcross(i,:)=tmp;
end

sum(Tcross.cardio)
std(Tcross.age_bl)
mean(Tcross.age_bl)
                          
