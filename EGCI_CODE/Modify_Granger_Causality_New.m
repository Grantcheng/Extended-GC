function [Var_xy,Var_yx,hh]= Modify_Granger_Causality_New(pointset,detal,N,m1,m2,Taux,Tauy,Option)
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
Var_yx0=zeros(N,1);
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
K=randperm(length(pointset'));
ref=100;
pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');
referenceset=K(1:ref);
[count_r, neighbors] = range_search(pointset', atria, referenceset, detal, 0);
%% ---------------------------
for i=1:length(count_r)
    h_0=count_r(i);
if h_0>10*lm
     
    Next_Y_ix=neighbors{i,1}+Tauy.*ones(1,h_0);  
    Next_Y_iy=neighbors{i,1}+Taux.*ones(1,h_0);  
    
    index=find(Next_Y_ix<N+1); 
    Next_Y_ii=Next_Y_ix(index);
    
    index_y=find(Next_Y_iy<N+1); 
    Next_Y_ii_y=Next_Y_iy(index_y);
    
    Z0=pointset(:,referenceset(i));
    DATA=pointset(:,neighbors{i,1}); % neighbors{i,1} is the nearest points, neighbors{i,2}is the distance bewteen nearest points and reference points
    
    DD=pointset(:,Next_Y_ii);
    DD_y=pointset(:,Next_Y_ii_y);
    
    h=min(length(index),length(index_y));
    count=count+1;
    
    DA=DATA;
      
LX1=zeros(h,m); 
LX2=zeros(h,m2); 
LX2_1=zeros(h,m1); 
 
DATA1=DA(:,1:h);
DLT1=DD(:,1:h)';
DLT1_y=DD_y(:,1:h)';  
%---------------------------X->Y
   LY1=DLT1(1:h,m);
   LY2=DLT1(1:h,m);
%---------------------------Y->X   
   LY_1=DLT1_y(1:h,m1);
   LY_2=DLT1_y(1:h,m1);
%---------------------------
for j=1:h
    for k=1:m
    LX1(j,k)=DATA1(k,j);
    end
end
for j=1:h
    for k=m1+1:m
    LX2(j,k-m1)=DATA1(k,j);
    end
end

for j=1:h
    for k=1:m1
    LX2_1(j,k)=DATA1(k,j);
    end
end
number=count;
%---------------------------
Norma_Choose=1;

nobs_1 = size(LX1',2);
LX1_M = mean(LX1);
mall = repmat(LX1_M',1,nobs_1);
LX1 = LX1-mall';
LY1_M = mean(LY1);
mall = repmat(LY1_M',1,nobs_1);
if Norma_Choose==1
    Norma_02=LY1-mall';
    LY1 = (LY1-mall')/sqrt(Norma_02'*Norma_02);
elseif Norma_Choose==0
LY1 = LY1-mall';
end

nobs_2 = size(LX2',2);
LX2_M = mean(LX2);
mall = repmat(LX2_M',1,nobs_2);
LX2 = LX2-mall';
LY2_M = mean(LY2);
mall = repmat(LY2_M',1,nobs_2);
if Norma_Choose==1
Norma_2=LY2-mall';
LY2 = (LY2-mall')/sqrt(Norma_2'*Norma_2);
elseif Norma_Choose==0
    LY2=LY2-mall';
end

nobs_3 = size(LX2_1',2);
LX2_M = mean(LX2_1);
mall = repmat(LX2_M',1,nobs_3);
LX2_1 = LX2_1-mall';
LY2_M = mean(LY_2);
mall = repmat(LY2_M',1,nobs_3);
if Norma_Choose==1
    Norma_01=LY_2-mall';
    LY_2 = (LY_2-mall')/sqrt(Norma_01'*Norma_01);
elseif Norma_Choose==0
LY_2 = LY_2-mall';
end

LY2_M = mean(LY_1);
mall = repmat(LY2_M',1,nobs_3);
if Norma_Choose==1
Norma_1=LY_1-mall';
LY_1 = (LY_1-mall')/sqrt(Norma_1'*Norma_1);
elseif Norma_Choose==0
    LY_1 = (LY_1-mall');
end
nlags=m;
n2 = h-m;
%% -----------------------------------------------------------------
% LSE Old
% [Varxy]=LSE(LX1,LY1,h,m,number,Option,0); 
% [Vary]=LSE(LX2,LY2,h,m2,number,Option+4,0); 
% 
% [Varyx]=LSE(LX1,LY_1,h,m,number,Option,0); 
% [Varx]=LSE(LX2_1,LY_2,h,m2,number,Option+4,0); 

% LSE New

% display('Now is The GC_xy: ')
[Varxy,pxy,Flagxy,RSS0]=LSE_New(LX1,LY1,h,m,number,Option,0); 
[Vary,py,Flagy,RSS1]=LSE_New(LX2,LY2,h,m2,number,Option+4,0); 

ftest = ((RSS0-RSS1)/nlags)/(RSS1/n2);    % causality x->y
ret.prb= 1 - cca_cdff(ftest,nlags,n2);
if Vary==0
GC=1;
else
GC=1-(Varxy/Vary);
end
[PR,q] = cca_findsignificance(ret,PVAL,3);
Var(count,1) = GC.*PR;
if (Var(count,1))<0
    LX1;
    LX2;
    %display(['The GC_xy is small Zero: ' num2str(1)]); 
end
% display(' ')
% display('Now is The GC_yx: ')
[Varyx,pyx,Flagyx,RSS0_0]=LSE_New(LX1,LY_1,h,m,number,Option,0); 
[Varx,px,Flagx,RSS1_1]=LSE_New(LX2_1,LY_2,h,m1,number,Option+4,0); 

ftest_1 = ((RSS0_0-RSS1_1)/nlags)/(RSS1_1/n2);    % causality y->x
ret_1.prb= 1 - cca_cdff(ftest_1,nlags,n2);
if Varx==0
GC_1=1;
else
GC_1=1-(Varyx/Varx);
end
[PR_1,q_1] = cca_findsignificance(ret_1,PVAL,3);
Var_yx0(count,1) = GC_1.*PR_1;
if Var_yx0(count,1)<0
    LX1;
    LX2_1;
   % display(['The GC_yx is small Zero: ' num2str(1)]); 
end
% if Vary<1e-30 && Varxy<1e-30
% GC=1;
% else
% GC=1-(Varxy/Vary);
% end
Zhengding_xy(count,:)=[pxy py] ;
Zhengding_yx(count,:)=[pyx px];
MLL_xy(count,:)=[Flagxy Flagy];
MLL_yx(count,:)=[Flagyx Flagx];

clear DATA DA DD LX2 LX1 D1 DATA1
end
end
% INDEX=find(Var<=1&Var>0);
% Var_1=Var(INDEX);

Var_0=Var(1:count,1);
INDEX1=Var_0>=0;
INDEX2=Var_0<=1;
INDEX=INDEX1.*INDEX2;
Var_xy=Var_0.*INDEX

Var_00=Var_yx0(1:count,1);
INDEX10=Var_00>=0;
INDEX20=Var_00<=1;
INDEX0=INDEX10.*INDEX20;
Var_yx=Var_00.*INDEX0

%  Var_1=Var(1:count,1);

hh=length(Var_xy);
% saveas(gcf,['Overlap the attractor, neighborhood is  ' , num2str(detal) ],'png');
 
end
