% tic;
% clear;
% M=400:5:3000;  %膜面积单位是um^2
% N=10;  %channel个数
% g_ch=35*10^(-12); %g_ch单位是ps
% 
% V_l=-90*10^(-3); %V_l单位是mV
% V_ch=0*10^(-3); %V_ch单位是mV
% 
% A=10;       %浓度，单位是nM
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
% 
% EU=zeros(1,length(M));
% EU2=zeros(1,length(M));
% DU=zeros(1,length(M));
% EU2N=zeros(1,length(M));
% DUN=zeros(1,length(M));
% for j=1:length(M)
%     C_m=0.01*10^(-12)*M(j);  %C_m单位是pF/um^2*um^2
%     g_l=3*10^(-12)*M(j);%g_l单位是ps/um^2*um^2 
%     sigma=zeros(N+1);
%     for i=1:N+1
%          sigma(i,i)= sqrt(2*K_b*T*(g_l+(i-1)*g_ch));
%     end
%     sigma_m = sigma^2/(C_m^2);
%   N_1=zeros(N+1);
%   N_1(1,1)=-N*k1;N_1(1,2)=k2;
%  for i=2:N
%          N_1(i,i)= -((N-i+1)*k1+(i-1)*k2);
%          N_1(i,i-1)=(N-i+2)*k1;
%          N_1(i,i+1)=i*k2; 
%  end
%  N_1(N+1,N)=k1;N_1(N+1,N+1)=-N*k2;
%  N_1;
% 
%   D_g=zeros(N+1);
%   for i=1:N+1
%          D_g(i,i)= -(g_l+(i-1)*g_ch)/C_m;
%   end
%   D_g;
% 
%   D_V=zeros(N+1);
%   for i=1:N+1
%         D_V(i,i)= (g_l*V_l+(i-1)*g_ch*V_ch)/(g_l+(i-1)*g_ch);
%   end
%   D_V;
% 
%   D_gV=D_g*D_V;
%   A_1=N_1+D_g;
% 
%   P = zeros(N+1,1);
%   for i=1:N+1
%         P(i)= (factorial(N) /(factorial(i-1)*factorial(N-i+1)))*(p^(i-1))*((1-p)^(N-i+1));
%   end
%   P;
%   A_2= N_1+2*D_g;
%   EU_m = pinv(A_1)*D_gV*P;
%   EU(j) = sum(EU_m);%一阶矩
%   EU2(j)=sum( pinv(A_2)*(2*D_gV*EU_m));
%   EU2N(j)=sum(pinv(A_2)*(-sigma_m*P));
%   DU(j) = EU2(j)-EU(j)^2;%电导方差
%   DUN(j)= EU2N(j);%噪声方差
%   DUS(j)= DU(j)+DUN(j);
%   
%   
% end
% save varwn_exachM_1010 M DUN

%%%%%%%%%%%%N固定，随M的变化DN_wn
load('varwn_exachM_150.mat')
 plot(M,DUN,'b');hold on 
load('varwn_exachM_1050.mat')
plot(M,DUN,'r');hold on 
load('varwn_exachM_2050.mat')
plot(M,DUN,'b');hold on 
load('varwn_exachM_1010.mat')
plot(M,DUN,'r');hold on 
% %load('varcn_ach_10100.mat')
% %plot(M,DU,'r');hold on 
x=2000:3000;
plot(x,1e-6*1./(x),'k--')
% 
% % 
