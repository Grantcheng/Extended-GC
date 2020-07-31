clc
clear
close all
tic;
addpath utilities/;
settspath;
% N = input('input N  ');
% DetalN=input('input DetalN  ');
N=50000;
PVAL=0.01;
time=0;
epstot=0.1:0.1:1;
Realization=length(epstot);

Var=zeros(1,Realization);RSS=zeros(1,Realization);Var1=zeros(1,Realization);RSS1=zeros(1,Realization);

Preprocessing=0;

C=1;


for r_i=1:Realization
tic;
    display(['The Couple is: ' num2str(r_i)]);
    if r_i==1
    V1=load('Amp2=0/Couple=0.1/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.1/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];      
    elseif r_i==2
    V1=load('Amp2=0/Couple=0.2/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.2/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];
    elseif r_i==3
    V1=load('Amp2=0/Couple=0.3/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.3/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];
    elseif r_i==4
    V1=load('Amp2=0/Couple=0.4/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.4/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];    
    elseif r_i==5
    V1=load('Amp2=0/Couple=0.5/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.5/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];
    elseif r_i==6
    V1=load('Amp2=0/Couple=0.6/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.6/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];
    elseif  r_i==7
    V1=load('Amp2=0/Couple=0.7/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.7/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];  
    elseif r_i==8
    V1=load('Amp2=0/Couple=0.8/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.8/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)]; 
    elseif r_i==9
    V1=load('Amp2=0/Couple=0.9/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=0.9/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];   
    elseif r_i==10
    V1=load('Amp2=0/Couple=1/HH_Solution_w1_0.28_0.txt');
    V2=load('Amp2=0/Couple=1/HH_Solution_w1_0.28_1.txt');
    X=[V1(1,end-N:end);V2(1,end-N:end)];  
    end
toc

clear V1;
clear V2;
%% ------------------Data Preprocessing
if Preprocessing==1
% detrend and demean data
disp('detrending and demeaning data');
X = cca_detrend(X);
X = cca_rm_temporalmean(X);

% check covariance stationarity
disp('checking for covariance stationarity ...');
uroot = cca_check_cov_stat(X,10);
inx = find(uroot);
if sum(uroot) == 0,
    disp('OK, data is covariance stationary by ADF');
else
    disp('WARNING, data is NOT covariance stationary by ADF');
    disp(['unit roots found in variables: ',num2str(inx(1,1))]);
end

% check covariance stationarity again using KPSS test
[kh,kpss] = cca_kpss(X);
inx = find(kh==0);
if isempty(inx),
    disp('OK, data is covariance stationary by KPSS');
else
    disp('WARNING, data is NOT covariance stationary by KPSS');
    disp(['unit roots found in variables: ',num2str(inx(1))]);
end
end
%% ------------------
len=length(X(1,:));
x = X(1,len-N+1:end);
y=X(2,len-N+1:end);
clear X;
LEN=length(x);
mm=0;
 %% -------Set the embedding dimension and time delay
m = 8       %  Embedding Dimension
R = 6    % Combinations of order R of different monomial terms.

%--------------------------------------------------
% Polynomial embedding vertor (PEV)
Degree_M=0;
for i=1:R
Degree_M=factorial(m+i-1)/(factorial(i)*factorial(m-1))+Degree_M;
end
xn=zeros(LEN-m,Degree_M);
yn=zeros(LEN-m,Degree_M);

Coe_Num=0;
tic;
for i=0:R 
    for j=0:R
        for k=0:R
            for l=0:R
                for mm=0:R
                    for n=0:R
                    for i_i=0:R
                        for j_j=0:R
            if((i+j+k+l+mm+n+i_i+j_j)>=1 && (i+j+k+l+mm+n+i_i+j_j)<=R)
            Coe_Num=Coe_Num+1;
            Coeiff(Coe_Num,:)=[i j k l mm n i_i j_j];
            end
                        end
                    end
                    end
                end
            end
        end
    end  
end
toc
%% Reconstruct phase attractor
x=fliplr(x);
y=fliplr(y);

tic;
for i=2:LEN-m
    
   Temp_X=(repmat(x(1,i:i+m-1),Coe_Num,1).^Coeiff);
   Temp_Y=(repmat(y(1,i:i+m-1),Coe_Num,1).^Coeiff);
   
   Temp_XX=Temp_X(:,1).*Temp_X(:,2);
   Temp_YY=Temp_Y(:,1).*Temp_Y(:,2);
   
   for j=3:size(Temp_X,2)
    Temp_XX=Temp_XX.*Temp_X(:,j);
    Temp_YY=Temp_YY.*Temp_Y(:,j);
   end
   
   xn(i,:)=Temp_XX';
   yn(i,:)=Temp_YY';

end
toc
clear Coeiff;
L1=size(xn,1);
L2=size(yn,1);
xn=xn';
yn=yn';
%% ------------------------------------------ 
ZX=zeros(2*Degree_M,L1);ZY=zeros(2*Degree_M,L2);

for i=1:L1
   ZX(:,i)=[xn(:,i);yn(:,i)];%��X��Y ����
   ZY(:,i)=[yn(:,i);xn(:,i)];%��Y��X ����
end
clear xn;
clear yn;
%-------------------------------------------------------
LY1=y(1,1:LEN-m)';
LY2=x(1,1:LEN-m)';

clear x;
clear y;

LX1=ZX';
LX2=ZY';

clear ZX;
clear ZY;

O_LX1=LX1(:,Degree_M+1:end);
O_LX2=LX2(:,Degree_M+1:end);


nobs_1 = size(LX1',2);
LX1_M = mean(LX1);
mall = repmat(LX1_M',1,nobs_1);
LX1 = LX1-mall';
LY1_M = mean(LY1);
mall = repmat(LY1_M',1,nobs_1);
LY1 = LY1-mall';

nobs_2 = size(LX2',2);
LX2_M = mean(LX2);
mall = repmat(LX2_M',1,nobs_2);
LX2 = LX2-mall';
LY2_M = mean(LY2);
mall = repmat(LY2_M',1,nobs_2);
LY2 = LY2-mall';

nobs_3 = size(O_LX1',2);
LX1_M = mean(O_LX1);
mall = repmat(LX1_M',1,nobs_3);
O_LX1 = O_LX1-mall';


nobs_4 = size(O_LX2',2);
LX2_M = mean(O_LX2);
mall = repmat(LX2_M',1,nobs_4);
O_LX2 = O_LX2-mall';

%% -----------------------------------------------------------------
[Varxy]=LSE(LX1,LY1,L1,2*Degree_M,1,1,0); 
[Vary]=LSE(O_LX1,LY1,L1,Degree_M,1,1+4,0); 
nlags=Degree_M;
n2 =L1;
ftest = ((Vary-Varxy)/nlags)/(Varxy/n2);    % causality x->y
ret.prb= 1 - cca_cdff(ftest,nlags,n2);
GC=1-(Varxy/Vary);
[PR,q] = cca_findsignificance(ret,PVAL,3);
Var(1,r_i) = GC.*PR;

[Varxy]=LSE(LX2,LY2,L1,2*Degree_M,1,1,0); 
[Vary]=LSE(O_LX2,LY2,L1,Degree_M,1,1+4,0); 
nlags=Degree_M;
n2 =L1;
ftest = ((Vary-Varxy)/nlags)/(Varxy/n2);    % causality y->x
ret.prb= 1 - cca_cdff(ftest,nlags,n2);
GC=1-(Varxy/Vary);
[PR,q] = cca_findsignificance(ret,PVAL,3);
RSS(1,r_i) = GC.*PR;


%------------------------------------------------
% display('The total time cost is :');
% 
% %-------------------------------------------------
% if Realization==1
% figure;
% plot(M,Var(:,r_i),'r+',M,RSS(:,r_i),'b*');hold on;
% plot(M,Var(:,r_i),'r',M,RSS(:,r_i),'b');
% axis([0 1 0 1]);
% title(['Nonlinear driving case Granger causality N= ',num2str(L) 'DN= ',num2str(DetalN),' m= ',num2str(m)]); 
% xlabel('\delta ');           
% ylabel('Extended Granger Causality Index');          
% legend('x->y','y->x');
% end
% clear X;
    
end

figure;
plot(epstot,Var(:,r_i),'r+',epstot,RSS(:,r_i),'b*');hold on;
plot(epstot,Var(:,r_i),'r',epstot,RSS(:,r_i),'b');
axis([0 1 0 1]);
title(['Nonlinear driving case Granger causality R= ',num2str(R),' m= ',num2str(m)]); 
xlabel('c');           
ylabel('causality x->y');          
legend('x->y','y->x');
clear X;
% 
% if Realization>1
% %% -----Errorbar plot
% varyx2y=mean(Var');
% varyy2x=mean(RSS');
%  
% figure;
% plot(M,varyx2y,'r+',M,varyy2x,'b*');
% hold on;
% plot(M,varyx2y,'r',M,varyy2x,'b');
% ex2y = std(Var');
% ey2x = std(RSS');
% errorbar( M, varyx2y, ex2y,'r');
% errorbar(M,varyy2x,ey2x,'b');
% axis([0 1 0 1]);
% title(['Nonlinear driving case Granger causality N= ',num2str(L) 'DN= ',num2str(DetalN),' m= ',num2str(m)]);%��ͼ�α��� 
% xlabel('\delta ');           
% ylabel('Extended Granger Causality Index');          
% legend('x->y','y->x');
% end

