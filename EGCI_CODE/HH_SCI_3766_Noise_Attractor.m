% EGCI for f1=0.28 f2=0.282 S12=0.16 S21=0 Amp1=Amp2=10
clc
clear
% close all
tic;
display('The HH_3766 All Couple Program ......   ');

global Couple_ie;

epstot=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2]

time=0;
GC_Matrix=zeros(1,2);
data_chose=3766;
Choose_Wolf_No_Sampling=0;
Chose_One_Couple=input('Whether use all ie  Yes =  1,  No = 0  ');
Noise_Add=input('Whether Noise_Add: Yes =1 , No = 0 ');

if Noise_Add==1
sigma_white=input('The variance of white noise is :  ');
end

if Chose_One_Couple==0
ie_3=input('The Couple Index is = ');
Leng_IE=ie_3;
elseif Chose_One_Couple==1
Leng_IE=length(epstot);
    ie_3=input('The Start of Couple Index is = ');
end
        mx=input('The mx is = ');
        my=input('The my is = ');
        taux=input('Tau_xy is = ');
        tauy=taux;
display(['The Chose_One_Couple is: ' num2str(Chose_One_Couple)]);

display(['The data_chose is: ' num2str(data_chose)]);
for ie=ie_3:Leng_IE
    
    C=epstot(ie);
    Couple_ie=C;
     display(['The Strength is : ' num2str(C)]);
     display('Sampling is 0.01 ');
if ie==1
    V1=load('HH_3766/Couple=0/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0/HH_Solution_w1_0.28_1.txt');
   elseif ie==2
    V1=load('HH_3766/Couple=0.02/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.02/HH_Solution_w1_0.28_1.txt');   
    
    elseif ie==3
    V1=load('HH_3766/Couple=0.04/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.04/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==4
    V1=load('HH_3766/Couple=0.06/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.06/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==5
    V1=load('HH_3766/Couple=0.08/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.08/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==6
    V1=load('HH_3766/Couple=0.1/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.1/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==7
     V1=load('HH_3766/Couple=0.11/HH_Solution_w1_0.28_0.txt');
     V2=load('HH_3766/Couple=0.11/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==8
    V1=load('HH_3766/Couple=0.12/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.12/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==9
    V1=load('HH_3766/Couple=0.13/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.13/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==10
    V1=load('HH_3766/Couple=0.14/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.14/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==11
    V1=load('HH_3766/Couple=0.15/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.15/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==12
    V1=load('HH_3766/Couple=0.16/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.16/HH_Solution_w1_0.28_1.txt');
       
    elseif ie==13
    V1=load('HH_3766/Couple=0.17/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.17/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==14
    V1=load('HH_3766/Couple=0.18/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.18/HH_Solution_w1_0.28_1.txt');
            
    elseif ie==15
    V1=load('HH_3766/Couple=0.19/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.19/HH_Solution_w1_0.28_1.txt');
     
    elseif  ie==16
    V1=load('HH_3766/Couple=0.2/HH_Solution_w1_0.28_0.txt');
    V2=load('HH_3766/Couple=0.2/HH_Solution_w1_0.28_1.txt'); 
end
Sampling=0.01;
stv=0.01;% the sample interval
dt=0.01;% the data interval
stv_L=stv/dt;
X=[V1;V2];
N=length(X(1,:));
display(['The N is: ' num2str(N)]);
clear V1;
clear V2;
x = X(1,:);
y = X(2,:);

x=(x-min(x))./(max(x)-min(x));
y=(y-min(y))./(max(y)-min(y));
X=[x;y];

 x=X(1,:);
 y=X(2,:);

if Noise_Add==1
max_k=length(x);

White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white.*randn(1,max_k);

x=White_Noise+x;
y=White_Noise1+y;
end
%--------------------------------------------------
% (phase space reconstruction)
[xn,dn1,L1] = PhaSpaRecon(x(1,1:10000000),taux,mx);
[yn,dn2,L2] = PhaSpaRecon(y(1,1:10000000),tauy,my);
%% Reconstruct phase attractor

L=min(L1,L2);

figure;
plot(xn(1,1:10000),xn(2,1:10000),'k.','MarkerSize',3);
title(['Driving attractor Example-A  m= ',num2str(mx)]); 
xlabel('V_{1}(n)'); 
ylabel('V_{1}(n-1)');  

figure;
plot(yn(1,1:10000),yn(2,1:10000),'k.','MarkerSize',3);
title(['Driven attractor Example-A  m= ',num2str(my)]); 
xlabel('V_{2}(n)'); 
ylabel('V_{2}(n-1)');

end