function [Var,pa,Flag，RSS] = LSE(X,Y,N,m,number,option,off)
%  fai=((X'*X)\(X'))*Y;
% e=ones(m,1);I=eye(m,m);ee=ones(m,1);
%N is the length of X ;
if off==1
switch option
    case 1
        mingzi='X��Y ����';
    case 2
        mingzi='Y��X ����';
    case 3
        mingzi='X��Y ������';
    case 4
        mingzi='X��Y ������';
        case 5
        mingzi='Y ����';
    case 6
        mingzi='X ����';
    case 7
        mingzi='Y ������';
    case 8
        mingzi='Y ������';
end
end
y1=zeros(N,1);
A=X'*X;
b=X'*Y;
[ra,pa]=chol(A);
display(['The pa is: ' num2str(pa)]);
[eigvec,eigval]=eig(A);
nn=size(A,1);
%% --- Chose Ling Hui Gui or Zhu Cheng Feng
chose_L_Z=2;

k_dcgxx=eigval(nn,nn)/eigval(1,1);

% if  abs(k_dcgxx)>=100
if eigval(1,1)<1e-5
    disp(' Multicollinearity ');  
    Flag=1;
if chose_L_Z==1;
%% ------------- Ling Hui Gui
    %-------Determin the k_data
    X01=zeros(m,1);wucha1=10^(-8);
    %fai1=CG(A,b,X01,2,wucha1) ;
    fai1=A\b;
    fai_a=eigvec'*fai1;
    e_a=zeros(N,1);
     
     e_a(1:N,1)=Y(1:N,1)-X(1:N,:)*fai1 ;
     
     RSS_a=e_a'*e_a;
     Var_a=covariance(e_a,e_a,N);
     k_data=Var_a/max(fai_a);
    %---------------------------    
    X0=zeros(m,1);wucha=10^(-8);
    
    % fai=CG(A+k_data.*eye(nn,nn),b,X0,2,wucha) ;
    
    B=A+k_data.*eye(nn,nn);
    fai=B\b;
    else
%% ------------- Zhu Cheng Feng Fen Xi
Z=X*eigvec;
eigval_ZCF=eig(A);
Jieshu=size(eigval,2);
eigval_ZCF_0=eigval_ZCF<1e-5;
Delete_GS=0;
for i=1:Jieshu-1
if eigval_ZCF_0(i)==1 && eigval_ZCF_0(i+1)==0
    Delete_GS=i;
end 
end
Z1=Z(:,Delete_GS+1:end);
Alpha=(Z1'*Z1)\(Z1'*Y);
fai=eigvec(:,Delete_GS+1:end)*Alpha;

end
else
    Flag=0;
    disp(' OK  '); 
    X0=zeros(m,1);wucha=10^(-8);
    % fai=CG(A,b,X0,2,wucha) ;
    fai=A\b;

end
%% ----------------------------------
%fai=A\b;
%fai=GMRES(A,b,X0,m,wucha);
 e=zeros(N,1);
 
     e(1:N,1)=Y(1:N,1)-X(1:N,:)*fai ;
     y1(1:N,1)=X(1:N,:)*fai;
  if off==1
  figure(6);
  clf  ;
  gcf=figure(6) ;
  plot(y1,e,'k+');%ɢ��ͼ�����Ƿ��������Իع�ĸ�˹-���Ʒ����
  saveas(gcf,[mingzi ,num2str(number) ],'jpg');
 end
  RSS=e'*e;
  % Var=RSS/(N-m);
  Var=covariance(e,e,N);
end

