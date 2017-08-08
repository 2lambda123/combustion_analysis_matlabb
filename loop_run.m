%%loop run
jjjj=0;
for jj=25:25:75
   for  jjj=0:5:15
       jjjj = jjjj + 1; 
       f_in = strcat('L', num2str(jj),'H',num2str(jjj),'.csv');
       f_out_1_1 = strcat('L', num2str(jj),'H',num2str(jjj),'_values.xlsx');
       f_out_1_2 = strcat('L', num2str(jj),'H',num2str(jjj),'_values.txt');
       f_out_2_1 = strcat('L', num2str(jj),'H',num2str(jjj),'_curves.xlsx');
       f_out_2_2 = strcat('L', num2str(jj),'H',num2str(jjj),'_curves.txt');
       run('convert.m');
       mmm(jjjj,:)=mm_v;
       clearvars -except jj jjj jjjj mmm
   end
end
fin=table(mmm);
writetable(fin,'fin_values.xlsx','Sheet',1)
