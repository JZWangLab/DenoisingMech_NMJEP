tic;
clear;
M=2;%膜面积单位是um^2
N=30;  %channel个数
g_ch=35*10^(-12); %g_ch单位是ps

V_l=-90*10^(-3); %V_l单位是mV
V_ch=0*10^(-3); %V_ch单位是mV

A=20;       %浓度，单位是nM
k11=2*0.1*A;   %k1=8.15/80=0.1 ms^(-1)\muM^(-1)
k22=0.1*A; 
k_1=8.150;     %%单位是ms^(-1)
k_2=2*8.150;    %%单位是ms^(-1)
alpha=0.714;     %%单位是ms^(-1)
beta=30.600;     %%单位是ms^(-1)

k1=10^(3)*(k11*k22*beta*alpha)/(alpha*(k_1*k_2+k11*k_2+k11*k22));
k2=10^(3)*alpha;%%单位是ms^(-1)
p=k1/(k1+k2);

K_b=1.380649*10^-23;%玻尔兹曼常数,单位为J/K，J为焦耳，K是热力学温度
T=294.15;%绝对温度，单位K

C_m=0.01*10^(-12)*M;  %C_m单位是pF/um^2*um^2
g_l=3*10^(-12)*M;%g_l单位是ps/um^2*um^2 

sigma=zeros(N+1);
 for i=1:N+1
         sigma(i,i)= sqrt(2*K_b*T*(g_l+(i-1)*g_ch));
 end
sigma_m = sigma^2/(C_m^2);

% EU=zeros(1,length());
% EU2=zeros(1,length());
% EU2N=zeros(1,length());

N_1=zeros(N+1);
N_1(1,1)=-N*k1;N_1(1,2)=k2;
 for i=2:N
         N_1(i,i)= -((N-i+1)*k1+(i-1)*k2);
         N_1(i,i-1)=(N-i+2)*k1;
         N_1(i,i+1)=i*k2; 
 end
 N_1(N+1,N)=k1;N_1(N+1,N+1)=-N*k2;
 N_1;
 
 D_g=zeros(N+1);
 for i=1:N+1
         D_g(i,i)= -(g_l+(i-1)*g_ch)/C_m;
 end
 D_g;
 
 AA1=N_1+D_g;
 AA2= N_1+2*D_g;
 

 D_V=zeros(N+1);
 for i=1:N+1
       D_V(i,i)= (g_l*V_l+(i-1)*g_ch*V_ch)/(g_l+(i-1)*g_ch);
 end
 D_V;

 D_gV=D_g*D_V;
 

 P = zeros(N+1,1);
 for i=1:N+1
        P(i)= (factorial(N) /(factorial(i-1)*factorial(N-i+1)))*(p^(i-1))*((1-p)^(N-i+1));
 end
 P;

 
 
 EU1 = pinv(AA1)*D_gV*P;
 EU2_noise = pinv(AA2)*-sigma_m*P;
 EU2 = pinv(AA2)*2*D_gV*EU1;
 

A1=zeros(N+1);
B1=zeros(N+1);
A2=zeros(N+1);
B2=zeros(N+1);
freq1=0.001:001:20;
freq2=20:500:1000000;
freq=[freq1,freq2];
%  Omega1=10;
%  Omega2=1000;
% freq1=0.001:0.001:Omega1;
% freq2=Omega1:1:Omega2;
% freq3=Omega2:200:10000;
% freq=[freq1,freq2,freq3];


PSSUM1=zeros(1,length(freq));
PSSUM2=zeros(1,length(freq));
PSSUM3=zeros(1,length(freq));

for  ii_freq=1:length(freq)
    
    A1=1i*freq(ii_freq)*eye(N+1)-N_1;  %A(iw)  
    B1=1i*freq(ii_freq)*eye(N+1)-AA1;  %B(iw)
    A2=-1i*freq(ii_freq)*eye(N+1)-N_1; %A(-iw)
    B2=-1i*freq(ii_freq)*eye(N+1)-AA1; %B(-iw)  
    
   PSSUM1(ii_freq)=ones(1,N+1)*(pinv(B1)+pinv(B2))*EU2;
   PSSUM2(ii_freq)=ones(1,N+1)*(pinv(B1)*pinv(A1)+pinv(B2)*pinv(A2))*D_gV*EU1;
   PSSUM3(ii_freq)=PSSUM1(ii_freq)-PSSUM2(ii_freq);%不加噪声部分的功率谱
   
   PSSUM4(ii_freq)=ones(1,N+1)*(pinv(B1)+pinv(B2))*EU2_noise;%噪声项的功率谱
   PSSUM(ii_freq)=PSSUM3(ii_freq)+PSSUM4(ii_freq);
   
   
end
% %save spectrumach_10-1-40 freq PSSUM3 PSSUM4 PSSUM 
% 
% % subplot(1,2,1)
% %load('spectrum3-6.mat')
plot(freq,PSSUM3,'-.r');hold on
plot(freq,PSSUM4,'-.b');hold on
plot(freq,PSSUM,'--k');hold on
d=PSSUM3-PSSUM4;
ix1 = find(d>0,1,'last'); 
x11=freq(ix1);
y11=PSSUM3(ix1);
plot(x11,y11,'gp');hold on
x=100000:200000;
plot(x,1./(x.^4));hold on
plot(x,1./(x.^2))
% title('Power Spectrum')
% subplot(1,2,2)
% load('spectrum3-60.mat')
% plot(freq,PSSUM3,'-.r');hold on
% plot(freq,PSSUM4,'--b');hold on
% plot(freq,PSSUM,'-k');hold on
% d=PSSUM3-PSSUM4;
% ix1 = find(d>0,1,'last'); 
% x11=freq(ix1);
% y11=PSSUM3(ix1);
% plot(x11,y11,'gp');hold on
% x=100:1000;
% plot(x,1./(x.^4));hold on
% plot(x,1./(x.^2))
% subplot(2,2,3)
% load('spectrum30-6.mat')
% plot(freq,PSSUM3,'-.r');hold on
% plot(freq,PSSUM4,'--b');hold on
% plot(freq,PSSUM,'-k');hold on
% d=PSSUM3-PSSUM4;
% ix1 = find(d>0,1,'last'); 
% x11=freq(ix1);
% y11=PSSUM3(ix1);
% plot(x11,y11,'gp');hold on
% x=10:100;
% plot(x,1./(x.^4));hold on
% plot(x,1./(x.^2))
% % title('Power Spectrum')
% subplot(2,2,4)
% load('spectrum30-60.mat')
% plot(freq,PSSUM3,'-.r');hold on
% plot(freq,PSSUM4,'--b');hold on
% plot(freq,PSSUM,'-k');hold on
% d=PSSUM3-PSSUM4;
% ix1 = find(d>0,1,'last'); 
% x11=freq(ix1);
% y11=PSSUM3(ix1);
% plot(x11,y11,'gp');hold on
% x=100:1000;
% plot(x,1./(x.^4));hold on
% plot(x,1./(x.^2))
% % title('Power Spectrum')
