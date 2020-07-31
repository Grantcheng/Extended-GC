clc
clear
close all
tic;
addpath utilities/;
settspath;
% N = input('input N  ');
% DetalN=input('input DetalN  ');
N=500000;
DetalN=30;
time=0;
Realization=1;
M=0.1:(1/DetalN):1;
Var=zeros(length(M),Realization);RSS=zeros(length(M),Realization);Var1=zeros(length(M),Realization);RSS1=zeros(length(M),Realization);
SVar=zeros(length(M),1);SRSS=zeros(length(M),1);SVar1=zeros(length(M),1);SRSS1=zeros(length(M),1);
Preprocessing=1;
 
for r_i=1:Realization
X=example1(N);
% load('dataxy1.mat')            % ��ɢ��� x,y,yy�� 2N����
%% ------------------Data Preprocessing
if Preprocessing==1
display('Data Preprocessing ');
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
x = X(1,:);
y=X(2,:);
yy=X(3,:);
mm=0;
 %% -------Set the embedding dimension and time delay
tau = 1      % Delay Time   
m =2       %  Embedding Dimension

%--------------------------------------------------
% �������е���ռ��ع� (phase space reconstruction)

[xn,dn1,L1] = PhaSpaRecon(x,tau,m);
[yn,dn2,L2] = PhaSpaRecon(y,tau,m);
[yyn,dn2,L3] = PhaSpaRecon(yy,tau,m);

%% Reconstruct phase attractor
xn1=xn(1,:);
yn1=yn(1,:);
yyn1=yyn(1,:);
L=length(xn1);
% Attractor(xn1,yn1,yyn1,L,m)

%% ------------------------------------------ 
ZX=zeros(2*m,L1);ZZX=zeros(m,L1);ZZY=zeros(m,L1);ZY=zeros(2*m,L1);ZYY=zeros(2*m,L1);
ZY1=zeros(2*m,L1);
for i=1:L1
   ZX(:,i)=[xn(:,i)',yn(:,i)'];%��X��Y ����
   ZY(:,i)=[yn(:,i)',xn(:,i)'];%��Y��X ����
   ZY1(:,i)=[yyn(:,i)',xn(:,i)'];%��Y��X ������
   ZYY(:,i)=[xn(:,i)',yyn(:,i)'];%��X��Y ������
end
%-------------------------------------------------------
%-------------------------------------------------
for u=1:length(M)
% u=15;
% figure(1);
%-------------------------------------------------------
tic;
[Varxy,Lxy]=Granger_Causality(ZX,M(u),L,m,m,tau,1);%��X��Y ����
time=toc+time;
toc;
  
tic;
[Varyx,Lyx]=Granger_Causality(ZY,M(u),L,m,m,tau,2);%��Y��X ����
toc;
time=toc+time;

tic;
 [Varxy1,Lxy1]=Granger_Causality(ZYY,M(u),L,m,m,tau,3);%��X��Y ������
toc;
time=time+toc;

tic;
 [Varyx1,Lyx1]=Granger_Causality(ZY1,M(u),L,m,m,tau,4);%��Y��X ������
toc;
time=time+toc;
%% -----------------------------------------
Var2=Varxy(1:Lxy,1);  %��X��Y ����
Var3=Varxy1(1:Lxy1,1);%��X��Y ������
RSS2=Varyx(1:Lyx,1);%��Y��X ����
RSSyx=Varyx1(1:Lyx1,1);%��Y��X ������

Var(u,r_i)=mean(Var2);%��X��Y ����
Var1(u,r_i)=mean(Var3);%��X��Y ������
RSS(u,r_i)=mean(RSS2);%��Y��X ����
RSS1(u,r_i)=mean(RSSyx);%��Y��X ������
%% ----------------------------------------

mm=mm+1
end
time=toc+time;
display('The total time cost is :');
time
%-------------------------------------------------
if Realization==1
figure;
plot(M,Var(:,r_i),'r+',M,RSS(:,r_i),'b*');hold on;
plot(M,Var(:,r_i),'r',M,RSS(:,r_i),'b');
axis([0 1 0 1]);
title(['Linear driving case Granger causality N= ',num2str(L) 'DN= ',num2str(DetalN) ' m= ',num2str(m)]);
xlabel('\delta ');           
ylabel('Extended Granger Causality Index');          
legend('x->y','y->x');

figure;
plot(M,Var1(:,r_i),'r+',M,RSS1(:,r_i),'b*');hold on;
plot(M,Var1(:,r_i),'r',M,RSS1(:,r_i),'b');
axis([0 1 0 1]);
title(['Nonlinear driving case Granger causality N= ',num2str(L) 'DN= ',num2str(DetalN),' m= ',num2str(m)]);%��ͼ�α��� 
xlabel('\delta ');           
ylabel('Extended Granger Causality Index');           
legend('x->y','y->x');
end
clear X;
 
end
 

if Realization>1
%% -----Errorbar plot
 varyx2y=mean(Var');
 varyx2yn=mean(Var1');
 varyy2x=mean(RSS');
 varyy2xn=mean(RSS1');
 
 figure;
plot(M,varyx2y,'r+',M,varyy2x,'b*');
hold on;
plot(M,varyx2y,'r',M,varyy2x,'b');
ex2y = std(Var');
ey2x = std(RSS');
errorbar( M, varyx2y, ex2y,'r');
errorbar(M,varyy2x,ey2x,'b');
axis([0 1 0 1]);
title(['Linear driving case Granger causality N= ',num2str(L) 'DN= ',num2str(DetalN) ' m= ',num2str(m)]);
xlabel('\delta ');           
ylabel('Extended Granger Causality Index');         
legend('x->y','y->x');

figure;
plot(M,varyx2yn,'r+',M,varyy2xn,'b*');
hold on;
plot(M,varyx2yn,'r',M,varyy2xn,'b');
ex2yn = std(Var1');
ey2xn = std(RSS1');
errorbar( M, varyx2yn, ex2yn,'r');
errorbar(M,varyy2xn,ey2xn,'b');
axis([0 1 0 1]);
title(['Nonlinear driving case Granger causality N= ',num2str(L) 'DN= ',num2str(DetalN),' m= ',num2str(m)]);%��ͼ�α��� 
xlabel('\delta ');           
ylabel('Extended Granger Causality Index');          
legend('x->y','y->x');
end

