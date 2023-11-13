% tic;
% clear;
% M=1;%膜面积单位是um^2
% N=1:170;  %channel个数
% g_ch=35*10^(-12); %g_ch单位是ps
% 
% V_l=-90*10^(-3); %V_l单位是mV
% V_ch=0*10^(-3); %V_ch单位是mV
% 
% A=20;       %浓度，单位是nM
% k11=2*0.1*A;   %k1=8.15/80=0.1 ms^(-1)\muM^(-1)
% k22=0.1*A; 
% k_1=8.150;     %%单位是ms^(-1)
% k_2=2*8.150;    %%单位是ms^(-1)
% alpha=0.714;     %%单位是ms^(-1)
% beta=30.600;     %%单位是ms^(-1)
% 
% k1=10^(3)*(k11*k22*beta*alpha)/(alpha*(k_1*k_2+k11*k_2+k11*k22));
% k2=10^(3)*alpha;%%单位是ms^(-1)
% p=k1/(k1+k2);
% 
% K_b=1.380649*10^-23;%玻尔兹曼常数,单位为J/K，J为焦耳，K是热力学温度
% T=294.15;%绝对温度，单位K
% 
% C_m=0.01*10^(-12)*M;  %C_m单位是pF/um^2*um^2
% g_l=3*10^(-12)*M;%g_l单位是ps/um^2*um^2 
% 
% EU=zeros(1,length(N));
% EU2=zeros(1,length(N));
% DU=zeros(1,length(N));
% EU2N=zeros(1,length(N));
% DUN=zeros(1,length(N));
% for j=1:length(N)
%   
%   sigma=zeros(N(j)+1);
%   for i=1:N(j)+1
%          sigma(i,i)= sqrt(2*K_b*T*(g_l+(i-1)*g_ch));
%   end
%   sigma_m = sigma^2/(C_m^2);
%   N_1=zeros(N(j)+1);
%   N_1(1,1)=-N(j)*k1;N_1(1,2)=k2;
%  for i=2:N(j)
%          N_1(i,i)= -((N(j)-i+1)*k1+(i-1)*k2);
%          N_1(i,i-1)=(N(j)-i+2)*k1;
%          N_1(i,i+1)=i*k2; 
%  end
%  N_1(N(j)+1,N(j))=k1;N_1(N(j)+1,N(j)+1)=-N(j)*k2;
%  N_1;
% 
%   D_g=zeros(N(j)+1);
%   for i=1:N(j)+1
%          D_g(i,i)= -(g_l+(i-1)*g_ch)/C_m;
%   end
%   D_g;
% 
%   D_V=zeros(N(j)+1);
%   for i=1:N(j)+1
%         D_V(i,i)= (g_l*V_l+(i-1)*g_ch*V_ch)/(g_l+(i-1)*g_ch);
%   end
%   D_V;
% 
%   D_gV=D_g*D_V;
%   A_1=N_1+D_g;
% 
%   P = zeros(N(j)+1,1);
%   for i=1:N(j)+1
%         P(i)= (factorial(N(j)) /(factorial(i-1)*factorial(N(j)-i+1)))*(p^(i-1))*((1-p)^(N(j)-i+1));
%   end
%   P;
%   A_2= N_1+2*D_g;
%   EU_m = pinv(A_1)*D_gV*P;
%   EU(j) = sum(EU_m);%一阶矩
%   EU2(j) =sum( pinv(A_2)*(2*D_gV*EU_m));%二阶矩
%   EU2N(j)=sum(pinv(A_2)*(-sigma_m*P));
%   DU(j) = EU2(j)-EU(j)^2;%方差
%   DUN(j)= EU2N(j);%噪声
%   
%   
%   
% end
% save varcn_achN_201 N DU



%%%%%%%%%%%%%%%%%% M固定，随N变化的D_cn
load('varcn_achN_1020.mat')
plot(N,DU,'b');hold on 
load('varcn_achN_10100.mat')
plot(N,DU,'b');hold on 
load('varcn_achN_2050.mat')
plot(N,DU,'r');hold on 
load('varcn_achN_1050.mat')
plot(N,DU,'r');hold on 
load('varcn_achN_201.mat')
plot(N,DU,'k');hold on 
x=140:170;
plot(x,1./(x.^3),'k--')
% 
% 
% 
