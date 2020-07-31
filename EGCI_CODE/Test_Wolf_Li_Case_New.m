
%% Songting Li,.et.al 2018 PRE
% gc_c2r3_e1.m
clc;
close all;
tic;
clear;

addpath utilities/;
addpath('~\GC_clean-mainline\GCcal');
addpath('~\GC_clean-mainline\tools_and_utilities');
settspath;
Choose_Ex=input('Which Example Will Choose =  ');
 if Choose_Ex==7
     Choose_Wolf=input('Which Lya Method Will Choose: 1. Wolf   2.lyaprosen x  3.lyaprosen y  4.lyaprosen z ');
 else
Choose_Wolf=input('Which Lya Method Will Choose: 1. Wolf   2.lyaprosen x  3.lyaprosen y  ');
 end
if Choose_Wolf==1
    len = 1e5;
    Choose_Wolf_x=input('Which Time Series Will Choose: 1. X Time Series; 2. Y Time Series ');
    if Choose_Ex==7
       Choose_Wolf_x=input('Which Time Series Will Choose: 1. X Time Series; 2. Y Time Series  3. Z Time Series');     
    end
m=input('The Embedding Dimension =  ');
tau=input('The Time Delay = ');
else
    len=1e4;
end

global Ref_Number;
global test_cond_num;
global Couple_ie;
global Reference_Points_Choose;

test_cond_num=0;
Ref_Number=100;
Reference_Points_Choose=0;
Couple_ie=1;
GC_Matrix=zeros(1,2);
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
elseif Choose_Ex==6
    str_name=(['GC_DMI_Results/ex4_logistic_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end  
elseif Choose_Ex==7
  str_name=(['GC_DMI_Results/ex7_Conditional_case_',num2str(1),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end     
end
Resultsmat=[pwd,'/',str_name];

%% Generate the data
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
elseif Choose_Ex==6
A = zeros(2, 2*9);
[a, de] = gen_hfreq_coef(0.8, 0.1, 8);
A(1, 1:2:2*8) = a;
A(2, 1:2:2*9) = 100*[1 a];
A(2, 2:2:2*8) = a;
De = [de, 0.0; 0.0, 1.0];
X0 = gendata_linear(A, De, len+1e4);
X0 = X0(:, 1e4+1:end);

X0(1,:)=(X0(1,:)-min(X0(1,:)))./(max(X0(1,:))-min(X0(1,:)));
X0(2,:)=(X0(2,:)-min(X0(2,:)))./(max(X0(2,:))-min(X0(2,:)));

x = 3.98.*X0(1,:).*(1-X0(1,:));
y = X0(2,:);
elseif Choose_Ex==7 % Three Time Series Case
    c=0.24;
    Transfer=0;
    seed_experiment=load('seed.mat');
seed=seed_experiment.seed;
    nr=len;
rng(seed(1,1),'twister'); 
 e1=randn(1,nr);
 rng(seed(1,2),'twister'); 
 e2=randn(1,nr);
 rng(seed(1,3),'twister'); 
 e3=randn(1,nr);
  
  eps=zeros(1,nr);
  eps2=zeros(1,nr);
  eps3=zeros(1,nr);
  zz=zeros(1,nr);
  xx=zeros(1,nr);
  yy=zeros(1,nr);
 
     for i=2:nr
    eps2(1,i)=e2(i);
    yy(1,i)=0.9*yy(1,i-1)+eps2(1,i);
    end
    
    for i=2:nr
    eps3(1,i)=sqrt(0.2)*e3(1,i);
    zz(1,i)=0.5*zz(1,i-1)+0.5*yy(1,i-1)^2+eps3(1,i);
    end
   
   for i=2:nr
    eps(1,i)=sqrt(0.3)*e1(1,i);
    xx(1,i)=0.8*xx(1,i-1)+0.4*zz(1,i-1)+c*yy(1,i-1)^2+eps(1,i);
   end
if Transfer==1
x = 3.98.*xx(1,:).*(1-xx(1,:));
y = yy(1,:);
z=tanh(10*zz(1,:));
else
x = xx(1,:);
y = yy(1,:);
z=zz(1,:);   
end
end

x=(x-min(x))./(max(x)-min(x));
y=(y-min(y))./(max(y)-min(y));
if Choose_Ex==7
 z=(z-min(z))./(max(z)-min(z));   
end

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
if Choose_Ex==7
sigma_white_z=sqrt(var(z))/10^(3.2965)
sigma_white=(sigma_white_x+sigma_white_y+sigma_white_z)/2;

White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white.*randn(1,max_k);
White_Noise2=sigma_white.*randn(1,max_k);
else
sigma_white=(sigma_white_x+sigma_white_y)/2;
White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white.*randn(1,max_k);
end

end

if Choose_Ex==7
  x=White_Noise+x;
y=White_Noise1+y;
z=White_Noise2+z;
else
x=White_Noise+x;
y=White_Noise1+y;
end
end

if Choose_Ex==7
 if Choose_Wolf==1
    k1 = 30000;             % 前面的迭代点数
    k2 = 41000;              % 后面的迭代点数
    if Choose_Wolf_x==1
Lyapunov=lyapunov_wolf(x,k2,m,tau,19) % lyapunov_wolf(data,N,m,tau,P)
    elseif Choose_Wolf_x==2
  Lyapunov=lyapunov_wolf(y,k2,m,tau,19) % lyapunov_wolf(data,N,m,tau,P)      
    end
elseif Choose_Wolf==2
[LLE_x,Lambda_x,tau_x,m_x]=lyaprosen_New(x,0,0)
elseif Choose_Wolf==3
[LLE_y,Lambda_y,tau_y,m_y]=lyaprosen_New(y,0,0)
 elseif Choose_Wolf==4
[LLE_z,Lambda_z,tau_z,m_z]=lyaprosen_New(z,0,0)     
end   
else

if Choose_Wolf==1
    k1 = 30000;             % 前面的迭代点数
    k2 = 41000;              % 后面的迭代点数
    if Choose_Wolf_x==1
Lyapunov=lyapunov_wolf(x,k2,m,tau,19) % lyapunov_wolf(data,N,m,tau,P)
    elseif Choose_Wolf_x==2
  Lyapunov=lyapunov_wolf(y,k2,m,tau,19) % lyapunov_wolf(data,N,m,tau,P)      
    end
elseif Choose_Wolf==2
[LLE_x,Lambda_x,tau_x,m_x]=lyaprosen_New(x,0,0)
elseif Choose_Wolf==3
[LLE_y,Lambda_y,tau_y,m_y]=lyaprosen_New(y,0,0)
end
end
%Lyapunov=largest_lyapunov_exponent(x,k2,5,14,3)