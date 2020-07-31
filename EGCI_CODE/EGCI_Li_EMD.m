clc;
close all;
tic;
clear;

addpath utilities/;
addpath GC_GC_clean-mainline/GCcal/;
addpath GC_clean-mainline/tools_and_utilities/;

%addpath('~\GC_clean-mainline\GCcal');
%addpath('~\GC_clean-mainline\tools_and_utilities');
settspath;
Choose_Ex=input('Which Example Will Choose =  ');
global Ref_Number;
global test_cond_num;
global Couple_ie;
global Reference_Points_Choose;

test_cond_num=0;
Ref_Number=100;
Reference_Points_Choose=0;
Couple_ie=1;

if Choose_Ex==1
    str_name=(['EGCI_Chen_Results/gc_c2r1_e3_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end
elseif Choose_Ex==2
  str_name=(['EGCI_Chen_Results/gc_c2r3_e1_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end    
elseif Choose_Ex==3
    str_name=(['EGCI_Chen_Results/gc_rev_v3_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end 
elseif Choose_Ex==4
    str_name=(['EGCI_Chen_Results/gc_c1r1_e1_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end    
end
Resultsmat=[pwd,'/',str_name];

%% Generate the data
len = 1e5;
if Choose_Ex==1
x = zeros(len, 1);
y = zeros(len, 1);
%err = 1e-6 * rand(n, 2);

err0 = 0.3;
x(1) = err0 * rand(1);
x(2) = err0 * rand(1);
y(1) = err0 * rand(1);
y(2) = err0 * rand(1);

c = 0.5;

for k = 3:len
  x(k) = 3.4 * x(k-1) * (1 - x(k-1)^2) * exp(-x(k-1)^2) + 0.8 * x(k-2);
  y(k) = 3.4 * y(k-1) * (1 - y(k-1)^2) * exp(-y(k-1)^2) + 0.5 * y(k-2) + c * x(k-2);
end
%% gc_c2r3_e1
elseif Choose_Ex==2
c = 0.5;

odefun = @(t, y) [
  -0.25*y(1) + y(2) - y(2)^3;
  y(1) - y(2) - y(1)*y(2);
  -0.25*y(3) + y(4) - y(4)^3 + c * y(1)^2;
  y(3) - y(4) - y(3)*y(4)];
y0 = [1.6130    0.6268    4.3849    0.8178];
%[tt, yy] = ode45(odefun, linspace(0, 5*1000, 100*1000), y0);
yy = zeros(4, len);
yy(:,1) = y0;
dt = 0.05;
for j=2:len
  yy(:,j) = yy(:,j-1) + dt * odefun(j, yy(:,j-1)) + sqrt(dt)*0.01*rand(4,1);
end

x = yy(1, :)';
y = yy(3, :)';
%% gc_rev_v3
elseif Choose_Ex==3
c = 0.5;

odefunx = @(x1,y1,z1,x2,y2,z2) [
-(y1 + z1);
x1 + 0.2*y1;
0.2 + z1 * (x1 - 4.7);
-(y2 + z2) + c * x1;
x2 + 0.2 * y2;
0.2 + z2 * (x2 - 4.7) ];

odefun = @(y) odefunx(y(1),y(2),y(3),y(4),y(5),y(6));

y0 = randn(1, 6);

len = 10000;
yy = zeros(length(y0), len);

dt = 0.02;
yy(:, 1) = y0;
for j = 2:len
  yy(:, j) = yy(:, j-1) + dt * odefun(yy(:, j-1)) + sqrt(dt)*1e-1*rand(6,1);
end

x = yy(1, :);
y = yy(2, :);
z = yy(3, :);

%% gc_c1r1_e3
elseif Choose_Ex==4
len = 1e5;
x = zeros(1, len);
y = zeros(1, len);
z = zeros(1, len);

err = randn(3, len);

for k = 2 : len
  x(k) = 0.2 * x(k-1) + err(1, k);
  y(k) = 0.5 * y(k-1) + 0.5 * x(k-1) + err(2, k);
  z(k) = 0.4 * z(k-1) + 0.3 * y(k-1) + c * x(k-1) + err(3, k);
end
end

x=(x-min(x))./(max(x)-min(x));
y=(y-min(y))./(max(y)-min(y));

%% Add the Noise?
  Noise_Add=input('The Noise is Add? 1:Add  2:Non_Add ');% 
if Noise_Add==1
    SNR_Determined=input('The SNR is Determined :  ');
if SNR_Determined==0
sigma_white=input('The variance of white noise is :  ');
end

max_k=length(x);
if  SNR_Determined==1
 %sigma_white=sqrt(var(y))/10^(3)
sigma_white_x=sqrt(var(x))/10^(3.2965)
sigma_white_y=sqrt(var(y))/10^(3.2965)
sigma_white=(sigma_white_x+sigma_white_y)/2;
end
White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white.*randn(1,max_k);

x=White_Noise+x;
y=White_Noise1+y;
end
X = [x;y];
%% Extended GC
% Time Delay and Embeding Dimension
Auto_Tau=input('Auto_Tau is = ');
 if Auto_Tau==1
     taux=0;
     tauy=0;
s1=signal(x(1:50000,1));
a1=amutual2(s1,10);
a1_1=data(a1);
for i=1:length(a1_1)-2
if (a1_1(i+1)<a1_1(i)) && (a1_1(i+1)<a1_1(i+2))
taux=i+1
break;
end
end
if taux==0
a1=amutual2(s1,100);
a1_1=data(a1);
for i=1:length(a1_1)-2
if (a1_1(i+1)<a1_1(i)) && (a1_1(i+1)<a1_1(i+2))
taux=i+1
break;
end
end 
end

s2=signal(y(1:50000,1));
a2=amutual2(s2,40);
a2_1=data(a2);
for i=1:length(a2_1)-2
if (a2_1(i+1)<a2_1(i)) && (a2_1(i+1)<a2_1(i+2))
tauy=i+1
break;
end
end
if tauy==0
a2=amutual2(s2,1000);
a2_1=data(a2);
for i=1:length(a2_1)-2
   if (a2_1(i+1)<a2_1(i)) && (a2_1(i+1)<a2_1(i+2))
     tauy=i+1
     break;
   end
end
    
end
c1 = cao(s1,30,taux,3,1000);c2 = cao(s2,18,tauy,3,1000);
figure;view(c1);figure;view(c2)
clear s1;
clear s2;

c11=data(c1);
c21=data(c2);

for i=1:length(c11)-1
EE(i,1)=var(c11(i:end));
end
for i=1:length(c21)-1
EE2(i,1)=var(c21(i:end));
end

HH=diff(EE);
HH2=diff(EE2);

for i=1:length(HH)
   if abs(HH(i))<10^(-4)
       mx=i+1
       break;
   end
    
end

for i=1:length(HH2)
   if abs(HH2(i))<10^(-4)
       my=i+1
       break;
   end
    
end

clear EE;
clear EE2;

 elseif Auto_Tau==2
    mx=input('The mx is = ');% 
    my=input('The my is = ');%
    Tau_xy=input('The Tau_xy is = ');%
    taux=Tau_xy
end
 
mx
m=mx+my;