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
    str_name=(['GC_DMI_Results/gc_c2r1_e3_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end
elseif Choose_Ex==2
  str_name=(['GC_DMI_Results/gc_c2r3_e1_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end    
elseif Choose_Ex==3
    str_name=(['GC_DMI_Results/gc_rev_v3_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end 
elseif Choose_Ex==4
    str_name=(['GC_DMI_Results/gc_c1r1_e1_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end 
elseif Choose_Ex==5
    str_name=(['GC_DMI_Results/gc_eq_noneq_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end 
    
end
Resultsmat=[pwd,'/',str_name];

%% Generate the data
len = 1e7;
if Choose_Ex==1
%% gc_c2r1_e3
A = [0.0, 0.0;
     0.1, 0.0];
De = [1.0, 0.0; 0.0, 1.0];
X0 = gendata_linear(A, De, len);

x = X0(1,:).^2;
y = X0(2,:);

%% gc_c2r3_e1
elseif Choose_Ex==2
A = [ 0.3, 0.0;
      0.9,-0.3];

De = [1.0, 0.0; 0.0, 1.0];
X0 = gendata_linear(A, De, len);

x = X0(1,:);
y = X0(2,:);

x = ((x+abs(x))/2).^5;
%% gc_rev_v3
elseif Choose_Ex==3
A = zeros(2, 2*9);
[a, de] = gen_hfreq_coef(0.8, 0.1, 8);
A(1, 1:2:2*8) = a;
A(2, 1:2:2*9) = 100*[1 a];
A(2, 2:2:2*8) = a;
De = [de, 0.0; 0.0, 1.0];
X0 = gendata_linear(A, De, len+1e4);
X0 = X0(:, 1e4+1:end);

x = X0(1,:);
y = X0(2,:);

x = x .^ 5;
%% gc_c1r1_e3
elseif Choose_Ex==4
A = [0.0, 0.1;
     0.1, 0.0];
De = [1.0, 0.0; 0.0, 1.0];
X0 = gendata_linear(A, De, len);

x = X0(1,:).^2;
y = X0(2,:);

%% gc_eq_noneq
elseif Choose_Ex==5
A = zeros(2, 2*9);
[a, de] = gen_hfreq_coef(0.8, 0.1, 8);
A(1, 1:2:2*8) = a;
A(2, 2:2:2*8) = a;
A(2, 1:2:2*9) = 0.5*[1 a];
A(1, 2:2:2*9) = 0.5*[1 a];
De = [de, 0.0; 0.0, de];
X0 = gendata_linear(A, De, len+1e4);
X0 = X0(:, 1e4+1:end);

x = X0(1,:);
y = X0(2,:);

x = tanh(10*x);
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
s1=signal(x(1,1:5000000)');
a1=amutual2(s1,40);
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

s2=signal(y(1,1:5000000)');
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