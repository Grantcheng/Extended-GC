% EGCI for f1=0.28 f2=0.282 S12=0.16 S21=0 Amp1=Amp2=10
clc
clear
% close all
tic;
display('The HH_28401 All Couple Program ......   ');
global Threshold_Spike_x;
global Threshold_Spike_y;
global Ref_Number;
global test_cond_num;
global Couple_ie;
global Reference_Points_Choose;
global Ration_S_NS;
global global_prb_gc;
global global_GC_Threshold;
addpath utilities/;
settspath;
% N = input('input N  ');
% DetalN=input('input DetalN  ');
epstot=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2]

time=0;
GC_Matrix=zeros(1,2);
data_chose=28401;
Choose_Wolf_No_Sampling=0;
Chose_One_Couple=input('Whether use all ie  Yes =  1,  No = 0  ');
Realization=input('Input the realization .... ');
test_cond_num=input('The test_cond_num is : Yes = 1, No = 0 ');
Reference_Points_Choose=input('Reference_Points_Choose is : 0 - Random;  1 - Threshold_Spike_x;  2 - Step_Fix;  3 - Threshold_Spike_x and Threshold_Spike_y  ');
Noise_Add=input('Whether Noise_Add: Yes =1 , No = 0 ');

if Noise_Add==1
SNR_Determined=input('The SNR is Determined :  ');
if SNR_Determined==0
sigma_white=input('The variance of white noise is :  ');
end
end

if (Reference_Points_Choose==1 || Reference_Points_Choose==3)
Ration_S_NS=input('Ration_S_NS =  ');
end

if Chose_One_Couple==0
ie_3=input('The Couple Index is = ');
Leng_IE=ie_3;
elseif Chose_One_Couple==1
Leng_IE=length(epstot);
    ie_3=input('The Start of Couple Index is = ');
end
if Chose_One_Couple==0
All_Neighbor=input('All_Neighbor  Yes = 1, No = 0 ');
if All_Neighbor==1
DetalN=10;
M=0.1:(1/DetalN):1;
else
DetalN=1;
M=0.1; 
end
else
DetalN=1;
M=0.1;    
end

Var=zeros(length(M),Realization);RSS=zeros(length(M),Realization);Var1=zeros(length(M),Realization);RSS1=zeros(length(M),Realization);
SVar=zeros(length(M),1);SRSS=zeros(length(M),1);SVar1=zeros(length(M),1);SRSS1=zeros(length(M),1);
Preprocessing=0;

ie_1=input('The mx and my reset = ');
if ie_1==0
mx_vect=9;
my_vect=7;
Tau_xy=29;
elseif ie_1==1
        mx_vect=input('The mx is = ');
        my_vect=input('The my is = ');
        Tau_xy=input('Tau_xy is = ');
end
Auto_Tau=input('Auto_Tau is = ');
Ref_Number=input('Ref_Number is = ');
display('The New Code In EGCI ');
display(['The Chose_One_Couple is: ' num2str(Chose_One_Couple)]);

display(['The data_chose is: ' num2str(data_chose)]);
for ie=ie_3:Leng_IE
    mx=0;
    my=0;
    taux=0;
    tauy=0;
    C=epstot(ie);
    Couple_ie=C;
display(['The Strength is : ' num2str(C)]);

if Choose_Wolf_No_Sampling==1    
      display('Sampling is 0.05 ');
V1=load('Couple=0.17/HH_Solution_w1_0.28_0.txt');
V2=load('Couple=0.17/HH_Solution_w1_0.28_1.txt');

Sampling=0.05;
stv=0.01;% the sample interval
dt=0.01;% the data interval
stv_L=stv/dt;
X=[V1(1,1:stv_L:length_N);V2(1,1:stv_L:length_N)];
else
     display('Sampling is 0.01 ');
     if data_chose==3766
V1=load('HH_3766/HH_Solution_w1_0.28_0.txt');
V2=load('HH_3766/HH_Solution_w1_0.28_1.txt'); 
     elseif data_chose==28401
if ie==1
    V1=load('HH_28401/Couple=0/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0/HH_Solution_w1_0.28_1.txt');
   elseif ie==2
    V1=load('HH_28401/Couple=0.02/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.02/HH_Solution_w1_0.28_1.txt');   
    
    elseif ie==3
    V1=load('HH_28401/Couple=0.04/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.04/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==4
    V1=load('HH_28401/Couple=0.06/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.06/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==5
    V1=load('HH_28401/Couple=0.08/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.08/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==6
    V1=load('HH_28401/Couple=0.1/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.1/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==7
     V1=load('HH_28401/Couple=0.11/HH_Solution_w1_0.28_0.txt');
     V2=load('HH_28401/Couple=0.11/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==8
    V1=load('HH_28401/Couple=0.12/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.12/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==9
    V1=load('HH_28401/Couple=0.13/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.13/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==10
    V1=load('HH_28401/Couple=0.14/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.14/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==11
    V1=load('HH_28401/Couple=0.15/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.15/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==12
    V1=load('HH_28401/Couple=0.16/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.16/HH_Solution_w1_0.28_1.txt');
       
    elseif ie==13
    V1=load('HH_28401/Couple=0.17/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.17/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==14
    V1=load('HH_28401/Couple=0.18/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.18/HH_Solution_w1_0.28_1.txt');
            
    elseif ie==15
    V1=load('HH_28401/Couple=0.19/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.19/HH_Solution_w1_0.28_1.txt');
     
    elseif  ie==16
    V1=load('HH_28401/Couple=0.2/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_28401/Couple=0.2/HH_Solution_w1_0.28_1.txt'); 
end
     elseif data_chose==6
V1=load('Amp2=6/Couple=0.17/HH_Solution_w1_0.28_0.txt');
V2=load('Amp2=6/Couple=0.17/HH_Solution_w1_0.28_1.txt');
V1=V1(1,1000000:6000001);
V2=V2(1,1000000:6000001);
     end
Sampling=0.01;
stv=0.01;% the sample interval
dt=0.01;% the data interval
stv_L=stv/dt;
X=[V1;V2];


end
N=length(X(1,:));
display(['The N is: ' num2str(N)]);
clear V1;
clear V2;
x = X(1,:);
y = X(2,:);
Threshold_Spike_x=(-55-min(x))./(max(x)-min(x));
Threshold_Spike_y=(-55-min(y))./(max(y)-min(y));

x=(x-min(x))./(max(x)-min(x));
y=(y-min(y))./(max(y)-min(y));
X=[x;y];
display(['The Auto_Tau is: ' num2str(Auto_Tau)]);
for mx_i=1:length(mx_vect)
    for my_i=1:length(my_vect)
for r_i=1:Realization
 x=X(1,:);
 y=X(2,:);

if Noise_Add==1
max_k=length(x);
if  SNR_Determined==1
 %sigma_white=sqrt(var(y))/10^(3)
sigma_white=sqrt(var(y))/10^(3.2965)
end
White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white.*randn(1,max_k);

x=White_Noise+x;
y=White_Noise1+y;
end


%% ------------------len=length(X(1,:));
mm=0;
    %% -------Set the embedding dimension and time delay
 if Auto_Tau==1
     taux=0;
     tauy=0;
s1=signal(x(1,1:1000000)');
a1=amutual2(s1,250);
a1_1=data(a1);
for i=1:length(a1_1)-2
if (a1_1(i+1)<a1_1(i)) && (a1_1(i+1)<a1_1(i+2))
taux=i
break;
end
end
if taux==0
a1=amutual2(s1,400);
a1_1=data(a1);
for i=1:length(a1_1)-2
if (a1_1(i+1)<a1_1(i)) && (a1_1(i+1)<a1_1(i+2))
taux=i
break;
end
end 
end

s2=signal(y(1,1:1000000)');
a2=amutual2(s2,250);
a2_1=data(a2);
for i=1:length(a2_1)-2
if (a2_1(i+1)<a2_1(i)) && (a2_1(i+1)<a2_1(i+2))
tauy=i
break;
end
end
if tauy==0
a2=amutual2(s2,350);
a2_1=data(a2);
for i=1:length(a2_1)-2
if (a2_1(i+1)<a2_1(i)) && (a2_1(i+1)<a2_1(i+2))
tauy=i
break;
end
end
    
end
c1 = cao(s1,18,taux,3,1000);c2 = cao(s2,18,taux,3,1000);
figure;view(c1);figure;view(c2)
clear s1;
clear s2;

c11=data(c1);
c21=data(c2);
% for i=1:length(c11)-1
% EE(i,1)=c11(i+1)/c11(i);
% end
% for i=1:length(c21)-1
% EE2(i,1)=c21(i+1)/c21(i);
% end
% Temp=abs(EE-ones(length(EE),1));
% Temp2=abs(EE2-ones(length(EE2),1));

% [min_mx,index_mx]=min(Temp);
% [min_my,index_my]=min(Temp2);
% mx = index_mx+1    %  Embedding Dimension
% my = index_my+1



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
    taux=Tau_xy
    mx=mx_vect(mx_i)
    my=my_vect(my_i)
end
 
tauy=taux % 41

if mx<8
   mx=8; 
end

mx
m=mx+my;
%--------------------------------------------------
% (phase space reconstruction)

[xn,dn1,L1] = PhaSpaRecon(x,taux,mx);
[yn,dn2,L2] = PhaSpaRecon(y,tauy,my);
%%
global Resultsmat;
if data_chose==3766 
    str_name=(['HH_3766/mx_',num2str(mx),' my_',num2str(my),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end
Resultsmat=[pwd,'/',str_name];

 elseif data_chose==28401
     if Noise_Add==1
         str_name=(['HH_28401/F_Test/Noise_Results_All/mx_',num2str(mx),' my_',num2str(my),'Sigma_',num2str(sigma_white),' Results mat']);    
     else
 str_name=(['HH_28401/F_Test/mx_',num2str(mx),' my_',num2str(my),' Results mat']);
     end
if ~exist(str_name)
    mkdir(str_name)
end
Resultsmat=[pwd,'/',str_name];

elseif data_chose==6
    str_name=(['Amp2=6/Couple=0.17/mx_',num2str(mx),' my_',num2str(my),' Results mat']);
if ~exist(str_name)
    mkdir(str_name)
end
Resultsmat=[pwd,'/',str_name];

end

%% Reconstruct phase attractor
% xn1=xn(1,:);
% yn1=yn(1,:);

L=min(L1,L2);
% Attractor(xn1,yn1,yyn1,L,m)

%% ------------------------------------------ 
ZX=zeros(mx+my,L);ZY=zeros(my+mx,L);

for i=1:L
   ZX(:,i)=[xn(:,i);yn(:,i)];%��X��Y ����
   ZY(:,i)=[yn(:,i);xn(:,i)];%��Y��X ����
end
clear xn;
clear yn;
%-------------------------------------------------------
%-------------------------------------------------
for u=1:length(M)
    display(['The Detal is: ' num2str(M(u))]);
% u=15;
% figure(1);
%-------------------------------------------------------
Choose_Modify_GC_New=1;
if Choose_Modify_GC_New==0
tic;
% [Varxy,Lxy]=Granger_Causality(ZX,M(u),L,mx,my,tau,1);%
[Varxy,Lxy]=Modify_Granger_Causality(ZX,M(u),L,mx,my,tauy,1);%

time=toc+time;
toc;
  
tic;
% [Varyx,Lyx]=Granger_Causality(ZY,M(u),L,my,mx,tau,2);%

[Varyx,Lyx]=Modify_Granger_Causality(ZY,M(u),L,my,mx,taux,2);%
toc;
time=toc+time;

time=time+toc;
elseif Choose_Modify_GC_New==1
    tic;

    [Varxy,Varyx,Lxy]=Modify_Granger_Causality_New_2(ZX,M(u),L,mx,my,taux,tauy,1);%
    Lyx=Lxy;
    toc;
    time=toc+time;
end
clear ZX;
%% -----------------------------------------
Var2=Varxy(1:Lxy,1);  %
RSS2=Varyx(1:Lyx,1);%

Var(u,r_i)=mean(Var2);%
RSS(u,r_i)=mean(RSS2);%

GC_Matrix(global_prb_gc<global_GC_Threshold) = 1;

Var(u,r_i)=GC_Matrix(1,1)*Var(u,r_i);%
RSS(u,r_i)=GC_Matrix(1,2)*RSS(u,r_i);%

MX_Temp(1,r_i)=mx;
MY_Temp(1,r_i)=my;

%% ----------------------------------------
display(['GC_xy = ',num2str(Var(u,r_i))]);
display(['GC_yx = ',num2str(RSS(u,r_i))]);

if Noise_Add==1
display(['SNR = ',num2str(log10(sqrt(var(y))/sigma_white))]);
end

mm=mm+1
end
time=toc+time;
display('The total time cost is :');
time
%-------------------------------------------------
%  fid=fopen('HH_28401/HH_28401_EGCI_Results.txt','a+');
if data_chose==3766 
      fid=fopen('HH_3766/HH_3766_EGCI_Noise_Results.txt','a+');

 elseif data_chose==28401
     fid=fopen('HH_28401/HH_28401_EGCI_Noise_Results.txt','a+');
     
elseif data_chose==6
     fid=fopen('Amp2=6/Couple=0.17/HH_EGCI_Noise_Results.txt','a+');
end

 fprintf(fid,'%s \t','mx');
 fprintf(fid,'%s \t','my');
 fprintf(fid,'%s \t','T_x');
 fprintf(fid,'%s \n','T_y');
 fprintf(fid,'%d \t',mx);
 fprintf(fid,'%d \t',my);
 fprintf(fid,'%d \t',taux);
 fprintf(fid,'%d \n',tauy);
 
 fprintf(fid,'%s \t','Var');
 for Var_ii=1:length(Var)
 fprintf(fid,'%d \t',Var(Var_ii));
     if Var_ii==length(Var)
          fprintf(fid,'%d \n',Var(Var_ii));
     end
 end
 fprintf(fid,'%s \t','RSS');
  for Var_ii=1:length(RSS)
 fprintf(fid,'%d \t',RSS(Var_ii));
     if Var_ii==length(RSS)
          fprintf(fid,'%d \n',RSS(Var_ii));
     end
  end
  fprintf(fid,'%s \n','-------------------------------------------------------------');
 fclose(fid);
 
if Realization==1
figure('visible','off');
plot(M,Var(:,r_i),'r+',M,RSS(:,r_i),'b*');hold on;
plot(M,Var(:,r_i),'r',M,RSS(:,r_i),'b');
axis([0 1 0 1]);
title(['H-H Model GC N= ',num2str(L) ,' S= ',num2str(C),' m_x= ',num2str(mx),' m_y= ',num2str(my),' taux= ',num2str(taux),' tauy= ',num2str(tauy), ' Sampling= ',num2str(Sampling)]); 
xlabel('\delta ');           
ylabel('Extended Granger Causality Index');          
legend('x->y','y->x');
end
    
end
    end
end

if Chose_One_Couple==1
Var_xy(ie,1:Realization)=Var;%
Var_yx(ie,1:Realization)=RSS;%
MX(ie,1:Realization)=MX_Temp;
MY(ie,1:Realization)=MY_Temp;
if Noise_Add==1
SNR(ie,1:Realization)=log10(sqrt(var(X(2,:)))/sigma_white);
end
Tau(ie,1:Realization)=taux;  
else
Var_xy(1,1:Realization)=Var;%
Var_yx(1,1:Realization)=RSS;%
MX=MX_Temp;
MY=MY_Temp;
if Noise_Add==1
SNR(1,1:Realization)=log10(sqrt(var(X(2,:)))/sigma_white);
end
Tau(1,1:Realization)=taux; 
    
end

end

display(['Noise_Add: ' num2str(Noise_Add)]);

