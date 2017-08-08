%%%This program reads a CSV File contaning Crank Angle (Index), In-Cylinder Pressuer and Fuel Line Pressure%%%

%%A=input ('Enter A'); % note enter "A" - case sensitive
% filename = 'L50H0.csv';
filename = f_in;
A=csvread(filename);
a=A(:,2); % value of all pressures 
b=A(:,3); % value of all fuel injection 
c=A(:,1);% value of all index 


x1=1:5000;
x=x1;
p1=a; 
% plot(x,p1)% TDC is visually aligned. 
p=p1;
%%% sub-routine to ensure that p always has positive value%%% 
[m,n]=size(p); 
for i=1:m;
%     mini=min(p(i,:)); 
    mini=min(p); 
    p(i,:)=p(i,:)-1*mini+.001; 
end 

%%% sub routine ends%%% 
pmean=mean(p); % Calculation of pmean 
x=1:5000; 
figure(4)
plot(x*.1,pmean) 
pr=p';
[pmax,theta]=max(pr); % Pmax and it's position

%%% sub-routine to ensure that injp always has positive value%%% 
[m,n]=size(b); 
for i=1:m;
%     mini=min(p(i,:)); 
    mini=min(b); 
    b(i,:)=b(i,:)-1*mini+.001; 
end 
inj =b;
cn = c+360;
dinj=gradient(inj,cn);
p=smooth(cn,p,10,'moving');
inj=smooth(cn,inj,10,'moving');
run('heatrelease_att.m');

figure(4)
yyaxis right
plot(cn,p);
yyaxis left
plot(cn,inj);
axis([cn(SOI_ind) eoc_idx*0.144 -inf inf])
% 
% %%Write Values%%
mm_n = {'cn_SOI', 'cn_EOI', 'Inj_dur', 'igdelay', 'cn_SOC', 'eoc', 'comb_dur', 'work', 'mep', 'gwork', 'gimep', 'pmep'};
mm_v =  [cn_SOI cn_EOI Inj_dur igdelay cn_SOC eoc comb_dur work mep gwork gimep pmep];
values = table (cn_SOI, cn_EOI, Inj_dur, igdelay, cn_SOC, eoc, comb_dur, work, mep, gwork, gimep, pmep);
% writetable(values,'L25H20.xlsx','Sheet',1)
% writetable(values,'L25H20.txt','Delimiter',' ')  
writetable(values,f_out_1_1,'Sheet',1)
writetable(values,f_out_1_2,'Delimiter',' ') 
%%Write values end%%

%%Write Cuves%%
cn_tdc_c=cn-360;
mm_calcu = [cn_tdc_c'; dpf'; avg_dpf'; heatrelt'; avg_heatrel'; ahr'; chr' ];
mm_calcu = mm_calcu.';
mm_calcu_name = {'cn_tdc_c'; 'dpf'; 'avg_dpf'; 'heatrelt'; 'avg_heatrel'; 'ahr'; 'chr' };
curves = table (cn_tdc_c, dpf, avg_dpf, heatrelt, avg_heatrel, ahr, chr);
% writetable(curves,'L25H20_curves.xlsx','Sheet',1)
% writetable(curves,'L25H20_curves.txt','Delimiter',' ')  
writetable(curves,f_out_2_1,'Sheet',1)
writetable(curves,f_out_2_2,'Delimiter',' ')  
% type('L25H0.txt')
