%% this programtakes Pressure matrix formed and averaged smoothed pressure data is returned

p=input('value of p'); % enter "p" where is "p" is the 
% pressure matrix formed by th.m
pavg=mean(ptemp); % pressure data averaged 
psmo=p;
[m,n]=size(p); 
i=0; 
for j=0:10:n-10 
n1=psmo(j+1:j+10); 
i=[i,mean(n1)]; 
end 
pdata=i(2:end); 
dp=diff(pdata); 
[i,j]=size(dp);
x=1:j; 
plot(x,dp) 
%%% Subroutine to ensure thatintake pressure is 1 Atm%%% 
if pdata(1)<1 
add=1-pdata(1); 
pdata=pdata+add; 
end 
% %% Subroutine ends%%% 
