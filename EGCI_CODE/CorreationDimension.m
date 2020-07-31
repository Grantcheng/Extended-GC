clc
clear
tic;
addpath utilities/;
settspath;

N_0=1000000;
time=0;
length_N=1000000;
Only_One_Couple=0;

epstot=0.17:0.001:0.2;
D2=zeros(length(epstot),1);

for ie=1:length(epstot)
 display(['The Couple is: ' num2str(epstot(ie))]);
    
    if Only_One_Couple==1
%       V1=load('S=0.1/ML_Solution_w1_0.029_0.txt');
%       V2=load('S=0.1/ML_Solution_w1_0.029_1.txt');
%     
     V1=load('Amp2=6/Couple=0.17/HH_Solution_w1_0.28_0.txt');
     V2=load('Amp2=6/Couple=0.17/HH_Solution_w1_0.28_1.txt');
       
    else
    if ie==1
    V1=load('Amp2=6/Couple=0.17/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.17/HH_Solution_w1_0.28_1.txt');
      
    elseif ie==2
    V1=load('Amp2=6/Couple=0.171/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.171/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==3
    V1=load('Amp2=6/Couple=0.172/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.172/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==4
    V1=load('Amp2=6/Couple=0.173/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.173/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==5
    V1=load('Amp2=6/Couple=0.174/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.174/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==6
    V1=load('Amp2=6/Couple=0.175/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.175/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==7
    V1=load('Amp2=6/Couple=0.176/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.176/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==8
    V1=load('Amp2=6/Couple=0.177/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.177/HH_Solution_w1_0.28_1.txt');
       
    elseif ie==9
    V1=load('Amp2=6/Couple=0.178/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.178/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==10
    V1=load('Amp2=6/Couple=0.179/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.179/HH_Solution_w1_0.28_1.txt');
            
    elseif ie==11
    V1=load('Amp2=6/Couple=0.18/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.18/HH_Solution_w1_0.28_1.txt');
     
    elseif  ie==12
    V1=load('Amp2=6/Couple=0.181/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.181/HH_Solution_w1_0.28_1.txt');
          
    elseif ie==13
    V1=load('Amp2=6/Couple=0.182/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.182/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==14
    V1=load('Amp2=6/Couple=0.183/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.183/HH_Solution_w1_0.28_1.txt');
   
     elseif ie==15
    V1=load('Amp2=6/Couple=0.184/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.184/HH_Solution_w1_0.28_1.txt');
     
     elseif ie==16
    V1=load('Amp2=6/Couple=0.185/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.185/HH_Solution_w1_0.28_1.txt');
      
     elseif ie==17
    V1=load('Amp2=6/Couple=0.186/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.186/HH_Solution_w1_0.28_1.txt');
    
     elseif ie==18
    V1=load('Amp2=6/Couple=0.187/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.187/HH_Solution_w1_0.28_1.txt');
    
     elseif ie==19
    V1=load('Amp2=6/Couple=0.188/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.188/HH_Solution_w1_0.28_1.txt');

    elseif ie==20
    V1=load('Amp2=6/Couple=0.189/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.189/HH_Solution_w1_0.28_1.txt');
       
     elseif ie==21
    V1=load('Amp2=6/Couple=0.19/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.19/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==22
    V1=load('Amp2=6/Couple=0.191/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.191/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==23
    V1=load('Amp2=6/Couple=0.192/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.192/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==24
    V1=load('Amp2=6/Couple=0.193/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.193/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==25
    V1=load('Amp2=6/Couple=0.194/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.194/HH_Solution_w1_0.28_1.txt');
         
    elseif ie==26
    V1=load('Amp2=6/Couple=0.195/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.195/HH_Solution_w1_0.28_1.txt');
       
    elseif ie==27
    V1=load('Amp2=6/Couple=0.196/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.196/HH_Solution_w1_0.28_1.txt');

    elseif ie==28
    V1=load('Amp2=6/Couple=0.197/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.197/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==29
    V1=load('Amp2=6/Couple=0.198/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.198/HH_Solution_w1_0.28_1.txt');
   
    elseif ie==30
    V1=load('Amp2=6/Couple=0.199/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.199/HH_Solution_w1_0.28_1.txt');  
    
    elseif ie==31
    V1=load('Amp2=6/Couple=0.2/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=6/Couple=0.2/HH_Solution_w1_0.28_1.txt');
    end
    end

stv=0.01;% the sample interval
dt=0.01;% the data interval
stv_L=stv/dt;

X=[V1(1,(N_0+1):stv_L:length_N+N_0);V2(1,(N_0+1):stv_L:length_N+N_0)];
N=length(X(1,:));

clear V1;
clear V2;

X(1,:)=(X(1,:)-min(X(1,:)))./(max(X(1,:))-min(X(1,:)));
X(2,:)=(X(2,:)-min(X(2,:)))./(max(X(2,:))-min(X(2,:)));

s=signal(X');
clear X;
display(['The C_D is: ' num2str(ie)]);
% D2 = takens_estimator2(s, n, range, past);
D2(ie,1) = takens_estimator(s, -1, 0.1, 1)
clear s;
end

