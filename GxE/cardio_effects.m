for i=1:68
%Outcome=(thicknessFS(:,15) + thicknessFS(:,49))/2; % 6 40 entorhinal %15 49 hippo %28 62 sup front %        Outcome=AD_ordered.x4282_2_0;
Outcome=thicknessFS(:,i);Outcome=clean(Outcome, 4); 
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %isnan(whitebritish_ordered)|
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV, Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 

mdl = fitlm(T,'Outcome~cardio+sex+age+HMotion+site+age*sex+age^2'); %plotInteraction(mdl,'sex','cardio','predictions'); ylabel('Cortical Thickness'); 
Cardio(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'cardio'),:);

mdl = fitlm(T,'Outcome~Dep+sex+age+HMotion+site+age*sex+age^2'); 
Dep(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Dep_1'),:);

end

for i=1:360;i
Outcome=ct_ordered(:,i);Outcome=clean(Outcome, 4); APOE=thickness_PRS'; APOE=CT1_PRS';
ix= AD_ordered.eid==0 |  ~strcmp(ancestries_ordered.pop,'EUR'); %isnan(whitebritish_ordered)|
T=table(APOE, PHQ2, age, sex, Hearing, HRT, AD_ordered.x6138_0_0, minimal_ordered.x25741_2_0, AD_ordered.x3581_0_0,AD_ordered.x20416_0_0, physical_ordered.x21001_2_0 ,physical_ordered.x6150_0_0,AD_ordered.x22506_0_0,minimal_ordered.x54_2_0, TIV, Outcome, 'VariableNames', {'APOE','Dep' ,'age', 'sex', 'Hearing', 'HRT', 'Education', 'HMotion','AgeLastPeriod', 'AlcFreq', 'BMI', 'cardio','tobacco','site', 'TIV', 'Outcome'}); T( ix,:)=[];%dep: strcmp(clinical, 'dep')' 

mdl = fitlm(T,'Outcome~cardio+sex+age+HMotion+site+age*sex+age^2'); %plotInteraction(mdl,'sex','cardio','predictions'); ylabel('Cortical Thickness'); 
Cardio(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'cardio'),:);

mdl = fitlm(T,'Outcome~Dep+sex+age+HMotion+site+age*sex+age^2'); 
Dep(i,:)=mdl.Coefficients(strcmp(mdl.Coefficients.Properties.RowNames,'Dep_1'),:);

end


ix=(Cardio.pValue<0.05/360); tmp=ct_ordered(:,ix); tmp=nanmean(tmp(physical_ordered.x6150_0_0==0,:)')';tmp(isnan(tmp))=[];
figure;plot_histogram_shaded(tmp,'Alpha',0.3, 'Normalization', 'pdf');%xlim([-1.5 3]);
ix=(Cardio.pValue<0.05/360); tmp=ct_ordered(:,ix); tmp=nanmean(tmp(physical_ordered.x6150_0_0==1,:)')';tmp(isnan(tmp))=[];
hold on;plot_histogram_shaded(tmp,'Alpha',0.3, 'Normalization', 'pdf');xlim([2.2 3.4]);
