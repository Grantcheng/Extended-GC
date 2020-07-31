function [Var_1,hh]= Granger_Causality(pointset,detal,N,m1,m2,Tau,Option)
%����
%   Detailed explanation goes here
% s Ϊʱ������
% detal Ϊ�������
% Distance Ϊ����֮��ľ���
% Z0 Ϊ�̶��ĵ�

PVAL=0.05;
m=m1+m2;
count=0;
Var=zeros(N,1);%VarΪ���յķ���
lm=m;       % The least points in sphere of radius

% if Option==1
%  figure(1) ;
% elseif Option==2 
%      figure(2) ;
% elseif Option==3 
%      figure(3) ;
% elseif Option==4 
%      figure(4) ;
% end
%  clf  ;

%% -------The number of reference points around the attractor is 100
K=randperm(length(pointset')-1);
ref=300;
pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');
referenceset=K(1:ref);
[count_r, neighbors] = range_search(pointset', atria, referenceset, detal, 0);
%% ---------------------------
for i=1:length(count_r)
    h_0=count_r(i);
if h_0>lm
    
    Next_Z0_i=referenceset(i)+Tau;
    
    if(Next_Z0_i<N+1)
    
    Next_Y_i=neighbors{i,1}+Tau.*ones(1,h_0);  
    index=find(Next_Y_i<N+1); 
    Next_Y_ii=Next_Y_i(index);
    Z0=pointset(:,referenceset(i));
    DATA=pointset(:,neighbors{i,1}); % neighbors{i,1} is the nearest points, neighbors{i,2}is the distance bewteen nearest points and reference points
    D0=pointset(:,Next_Z0_i);
    D=pointset(:,Next_Y_ii);
    else
        continue;
    end
    h=length(index);
    count=count+1;
    
    DA=[Z0,DATA];
    DD=[D0,D];
      
LX1=zeros(h+1,m1); 
LX2=zeros(h+1,m2); 
 
DATA1=DA(:,1:h+1);
D1=DD(:,1:h+1);
  
plot(DATA1(m,1:h),DATA1(m-1,1:h),'r.');  hold on;
plot(Z0(m,1),Z0(m-1,1),'b*');hold on;
title(['The neighborhood \delta  =' num2str(detal)]);
DLT1=D1';
%---------------------------��ȡY=XB�е�Y
   LY1=DLT1(1:h+1,m);
   LY2=DLT1(1:h+1,m);
%---------------------------��ȡY=XB�е�X
for j=1:h+1
    for k=1:m
    LX1(j,k)=DATA1(k,j);
    end
end
for j=1:h+1
    for k=m1+1:m
    LX2(j,k-m1)=DATA1(k,j);
    end
end
number=count;
%---------------------------
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

%% -----------------------------------------------------------------
[Varxy]=LSE(LX1,LY1,h+1,m,number,Option,0); 
[Vary]=LSE(LX2,LY2,h+1,m2,number,Option+4,0); 
nlags=m;
n2 = h+1;
ftest = ((Vary-Varxy)/nlags)/(Varxy/n2);    % causality x->y
ret.prb= 1 - cca_cdff(ftest,nlags,n2);
% if Vary<1e-30 && Varxy<1e-30
% GC=1;
% else
% GC=1-(Varxy/Vary);
% end
if Vary==0
GC=1;
else
GC=1-(Varxy/Vary);
end
[PR,q] = cca_findsignificance(ret,PVAL,3);
Var(count,1) = GC.*PR;
clear DATA DA DD LX2 LX1 D1 DATA1
end
end
% INDEX=find(Var<=1&Var>0);
% Var_1=Var(INDEX);

Var_0=Var(1:count,1);
INDEX1=Var_0>0;
INDEX2=Var_0<=1;
INDEX=INDEX1.*INDEX2;
Var_1=Var_0.*INDEX;

%  Var_1=Var(1:count,1);

hh=length(Var_1);%��ȡ��ĸ���
% saveas(gcf,['Overlap the attractor, neighborhood is  ' , num2str(detal) ],'png');
 
end