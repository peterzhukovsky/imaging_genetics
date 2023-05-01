map=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\chroms_viz\cytoBand.txt');
pos=readtable('D:\Canada_2020\UK_biobank\reports\AD\genetics\reports\summary_stats_ALL_cut.xlsx');
for i=1:length(pos.pos)
    ix= (map.Var2<pos.pos(i)) & (map.Var3>pos.pos(i) & map.Var6==pos.chr(i));
    towrite(i,1)=strcat(num2str(pos.chr(i) ), map.Var4(ix));
end
cd D:\Canada_2020\UK_biobank\reports\AD\genetics\EURnoWB\METAL\chroms_viz
writetable(table(towrite), 'CT_all_loci_cyto.txt')
