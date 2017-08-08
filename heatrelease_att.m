c_length=0.234; %%meter
r_throw=0.055;
r_stroke=0.110;
r_ratio=c_length/r_throw;
CR=17.5;
pis_area2=(pi*0.0875^2)/4;
vd2=(pi*0.0875^2*r_stroke)/4; 
vc2=(pis_area2*2*r_throw)/CR;
%%vc2=4.01e-05;
vol2=0;
cnr=(cn*pi)/180;

%%volume calculations%%
for i=1:1:5000
%     z1(i)=sqrt(c_length^2-((r_throw^2)*(sin(cnr(i,:))^2)));
%     z2(i)=c_length+r_throw-r_throw*cos(cnr(i))-z1(i);
%     v2(i,:)=vc2 + pis_area2*z2(i);
    ax(i,:)=r_ratio+(1-cos(cnr(i,:)))-sqrt(r_ratio^2-(sin(cnr(i,:))^2));
    v_theta(i,:)=1+((ax(i,:)*pis_area2*r_throw)/vc2);
    v_theta(i,:)=v_theta(i,:)*vc2;
    %%vol2=cat(2,vol2,v);
end
v2=v_theta;
% % [vmax,theta]=max(v2)
% % [vmin,theta]=min(v2)
% vol2=vol2(:,2:end); 
%%v2=transpose(v2); 

%%% Calculation Ends%%% 

dvol2=gradient(v2,cn);  % Volumechange for every crankangle
dpf=gradient(p,cn); % pressurechange for every crankangle
dcrankd=gradient(cn);

% averge to smooth PPR calculations
n = 10; % average every n values
% avg_dpf = arrayfun(@(i) mean(dpf(i:i+n-1)),1:n:length(dpf)-n+1)';% the averaged PPR
% avg_cn = arrayfun(@(i) mean(cn(i:i+n-1)),1:n:length(cn)-n+1)';
avg_dpf=smooth(cn,dpf,n,'moving');
clear n;
% calcuation end
%clear avg_cn

%%polytropic calculations
gamma_e=1.25;
poly_e=(1/(gamma_e-1));

gamma_c=1.15;
poly_c=(1/(gamma_c-1));


%Nmep=trapz(a,v2)

for ii=1:2500
    y1(ii)=(poly_c*gamma_c*p(ii,:)*dvol2(ii,:));
    y2(ii)=(poly_c*v2(ii,:)*dpf(ii,:));
%     heatrel(i,:)=y1(i+1)+y2(i+1);
end
heatrel_c=y1+y2;
clear y1 y2 ii;
for ii=2500:5000
    y1(ii)=(poly_e*gamma_e*p(ii,:)*dvol2(ii,:));
    y2(ii)=(poly_e*v2(ii,:)*dpf(ii,:));
%     heatrel(i,:)=y1(i+1)+y2(i+1);
    heatrel_c(ii)=0;
end
heatrel_e=y1+y2;
clear y1 y2 ii;
heatrel=heatrel_c+heatrel_e;
heatrel=heatrel*1000;

% averge to smooth heat release calculations
n = 10; % average every n values
% avg_heatrel = arrayfun(@(i) mean(heatrel(i:i+n-1)),1:n:length(heatrel)-n+1)';% the averaged vector
% avg_cn = arrayfun(@(i) mean(cn(i:i+n-1)),1:n:length(cn)-n+1)';
avg_heatrel=smooth(cn,heatrel,n,'moving');
% calcuation end
heatrel=avg_heatrel;
heatrelt=avg_heatrel;
hrr_mat=[cn, heatrelt];

% clear avg_cn;

%%SOI, SOC and EOC calculations start%%
%%SOI calcualtions
n = 10;
m = 4;
c = conv(double(inj >= m), ones(1, n));
%alter%   %c = conv((dat >= m)*1, ones(1, n))%
SOI_ind=min(find(c == n)) - n + 1;
clear n m

heatrelt_temp=heatrelt(SOI_ind:end);
cn_SOI=cn(SOI_ind)-360;

%%EOI calcualtions
[eoi_ind_min, EOI_ind]=min(dinj);
Inj_dur=EOI_ind-SOI_ind;
Inj_dur=Inj_dur*0.144;
 
clear n m heatrelt_temp

heatrelt_temp=heatrelt(SOI_ind:end);
cn_EOI=cn(SOI_ind)-360;

%%SOC calculations
n = 10;
m = 0;
c = conv(double(heatrelt_temp >= m), ones(1, n));
soc_ind=min(find(c == n)) - n + 1;
igdelay = soc_ind*0.144;
soc_ind=soc_ind+SOI_ind-1;
clear n m heatrelt_temp;
heatrelt_temp=heatrelt(soc_ind:end);
cn_SOC=cn(soc_ind)-360;



%%EOC calcualtions
z_v2=v2.^gamma_e;
z_eoc = p.*z_v2;
[eoc_val, eoc_idx]=max(z_eoc);
eoc=eoc_idx*0.144-360;
comb_dur=eoc-cn_SOC;

%%calcuations end%%

%%engien timing%%
ivc = -144.5;
evo = 144.5;
ivc_idx = round((ivc+360)/0.144);
evo_idx = round((evo+360)/0.144);
%%timing end%%

%Appearnet heat release calculations%
% heatrelt_ahr=heatrelt(soc_ind:eoc_idx);
% cn_ahr=cn(soc_ind:eoc_idx);

% heatrelt_ahr=heatrelt(EOI_ind:eoc_idx);
% cn_ahr=cn(EOI_ind:eoc_idx);
% ahr = cumtrapz(heatrelt_ahr,cn_ahr);
% chr= cumsum(heatrelt_ahr);

ahr = cumtrapz(heatrelt,cn);
chr = cumsum(heatrelt);
%Ahr calcuation end%

%Plots%
figure('Name','cummulative HR');
yyaxis left;
plot(cn,ahr);
yyaxis right;
plot(cn,chr);
axis([cn(SOI_ind), eoc_idx*0.144, -inf, inf])

figure('Name','Pressure and PRR');
yyaxis right
plot(cn,avg_dpf);
yyaxis left
plot(cn,p);
axis([cn(SOI_ind) eoc_idx*0.144 -inf inf])

figure('Name','heat release')
yyaxis left
plot (cn,heatrel);
yyaxis right;
% plot (avg_cn,avg_heatrel)
plot (cn,avg_heatrel)
axis([cn(SOI_ind) eoc_idx*0.144 -inf inf])
% heatrelt=transpose(heatrel);

%%workdone and mep calculations%%
% figure(10);
% plot(v2,p);
work=trapz(v2,p);%kN
mep=work/vd2;

gwork=trapz(v2(ivc_idx:evo_idx),p(ivc_idx:evo_idx)); %kN
gimep=gwork/vd2;
pmep=gimep-mep;

% %%Write Values%%
% mm_n = {'cn_SOI', 'cn_EOI', 'Inj_dur', 'igdelay', 'cn_SOC', 'eoc', 'comb_dur', 'work', 'mep', 'gwork', 'gimep', 'pmep'};
% mm_v =  [cn_SOI cn_EOI Inj_dur igdelay cn_SOC eoc comb_dur work mep gwork gimep pmep];
% values = table (cn_SOI, cn_EOI, Inj_dur, igdelay, cn_SOC, eoc, comb_dur, work, mep, gwork, gimep, pmep);
% writetable(values,'L25H0.xlsx','Sheet',1)
% writetable(values,'L25H0.txt','Delimiter',' ')  
% %%Write values end%%
% 
% %%Write Cuves%%
% cn_tdc_c=cn-360;
% mm_calcu = [cn_tdc_c'; dpf'; avg_dpf'; heatrelt'; avg_heatrel'; ahr'; chr' ];
% mm_calcu = mm_calcu.';
% mm_calcu_name = {'cn_tdc_c'; 'dpf'; 'avg_dpf'; 'heatrelt'; 'avg_heatrel'; 'ahr'; 'chr' };
% curves = table (cn_tdc_c, dpf, avg_dpf, heatrelt, avg_heatrel, ahr, chr);
% writetable(curves,'L25H0_curves.xlsx','Sheet',1)
% writetable(curves,'L25H0_curves.txt','Delimiter',' ')  
% % type('L25H0.txt')
