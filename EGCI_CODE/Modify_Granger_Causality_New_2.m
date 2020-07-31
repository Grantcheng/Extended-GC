function [Var_xy,Var_yx,hh]= Modify_Granger_Causality_New_2(pointset,detal,N,m1,m2,Taux,Tauy,Option)
global Resultsmat;
global Threshold_Spike_x;
global Ref_Number
global test_cond_num;
global Couple_ie;
global Reference_Points_Choose; % The Reference points are chose at spike region
global Ration_S_NS;
global global_GC_Threshold;
global global_prb_gc;
ref=Ref_Number;
 disp(['The Ref_Number is: ',num2str(ref)])
PVAL=0.01;
%test_cond_num=0;
GC_xy_Threshold=ones(1,ref).*NaN;
GC_yx_Threshold=ones(1,ref).*NaN;

S_NS=Ration_S_NS*ref
m=m1+m2;
count=0;
Var=zeros(N,1);%VarΪ���յķ���
Var_yx0=zeros(N,1);
GC_Method=Option; % 1-- F-Test; 2-- GC_Distribution
trails=50;
lm=m;       % The least points in sphere of radius
%% -------The number of reference points around the attractor is ref
if Reference_Points_Choose==0
K=randperm(length(pointset'));
pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');
referenceset=K(1:ref);
elseif Reference_Points_Choose==1
HH=pointset(1,:)>Threshold_Spike_x;

K=randperm(length(pointset'));

pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');

K2=(find(HH~=0));
K2_NS=(find(HH==0));
KK=randperm(length(K2'));
KK_NS=randperm(length(K2_NS'));
if S_NS==ref
    referenceset=K2(KK(1:S_NS));
else
referenceset0=K2(KK(1:S_NS));
referenceset1=K2_NS(KK_NS(1:(ref-S_NS)));
referenceset=[referenceset0 referenceset1];
end
elseif Reference_Points_Choose==3
HH=(pointset(1,:)>Threshold_Spike_x)&(pointset(2,:)>0.8);

K=randperm(length(pointset'));

pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');

K2=(find(HH~=0));
K2_NS=(find(HH==0));
KK=randperm(length(K2'));
KK_NS=randperm(length(K2_NS'));
if S_NS==ref
    referenceset=K2(KK(1:S_NS));
else
referenceset0=K2(KK(1:S_NS));
referenceset1=K2_NS(KK_NS(1:(ref-S_NS)));
referenceset=[referenceset0 referenceset1];
end
elseif Reference_Points_Choose==2
    Step=round(length(pointset')/ref);
    referenceset=1:Step:length(pointset');
    pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');
end
[count_r, neighbors] = range_search(pointset', atria, referenceset, detal, 0);
%% ---------------------------
for i=1:length(count_r)
    h_0=count_r(i);
        disp(['The Progress: ',num2str(i)])
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
    
    Z_0_Reserve(1:length(Z0),i)=Z0;
    Z_0_Reserve(1+length(Z0),i)=i;
    Z_0_Reserve(2+length(Z0),i)=referenceset(i);
    
    DA=DATA;
clear DATA;

LX1=zeros(h,m); 
LX2=zeros(h,m2); 
LX2_1=zeros(h,m1); 
 
DATA1=DA(:,1:h);
DLT1=DD(:,1:h)';
DLT1_y=DD_y(:,1:h)';  

clear DA; 
clear DD;
clear DD_y;
%---------------------------X->Y
   LY1=DLT1(1:h,m);
   LY2=DLT1(1:h,m);
   clear DLT1;
%---------------------------Y->X   
   LY_1=DLT1_y(1:h,m1);
   LY_2=DLT1_y(1:h,m1);
   clear DLT1_y;
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
% First Method to obtain GC
if GC_Method==1
% LSE New

% display('Now is The GC_xy: ')
beta_xy=LX1\LY1;
xpred_xy = LX1*beta_xy;
u_xy = LY1-xpred_xy;
Varxy=var(u_xy,1);
RSS0=u_xy'*u_xy;

beta_y=LX2\LY2;
xpred_y = LX2*beta_y;
u_y = LY2-xpred_y;
Vary=var(u_y,1);
RSS1=u_y'*u_y;
if Vary==0
GC=1;
else
GC=1-(Varxy/Vary);
end

if test_cond_num==1
[Varxy_test,pxy,Flagxy,RSS0_test]=LSE_New(LX1,LY1,h,m,number,1,0); % non-restricted regressions
[Vary_test,py,Flagy,RSS1_test]=LSE_New(LX2,LY2,h,m2,number,1,0); % restricted regressions

if Vary_test==0
GC_test=1;
else
GC_test=1-(Varxy_test/Vary_test);
end
end

% display('Now is The GC_yx: ')

beta_yx=LX1\LY_1;
xpred_yx = LX1*beta_yx;
u_yx = LY_1-xpred_yx;
Varyx=var(u_yx,1);
RSS0_0=u_yx'*u_yx;

beta_x=LX2_1\LY_2;
xpred_x = LX2_1*beta_x;
u_x = LY_2-xpred_x;
Varx=var(u_x,1);
RSS1_1=u_x'*u_x;
if Varx==0
GC_1=1;
else
GC_1=1-(Varyx/Varx);
end
if test_cond_num==1
[Varyx_test,pyx,Flagyx,RSS0_0_test]=LSE_New(LX1,LY_1,h,m,number,1,0); 
[Varx_test,px,Flagx,RSS1_1_test]=LSE_New(LX2_1,LY_2,h,m1,number,1,0); 
if Varx_test==0
GC_1_test=1;
else
GC_1_test=1-(Varyx_test/Varx_test);
end
end

%% GC_xy <0 or GC_yx<0
if GC<0
    disp('GC<0')
end

if GC_1<0
    disp('GC<0')
end

% Test Suppose
prb = ones(2).*NaN;
 
ftest = ((RSS1-RSS0)/m1)/(RSS0/n2);    % causality x->y
prb_xy= 1 - cca_cdff(ftest,m1,n2);

ftest_1 = ((RSS1_1-RSS0_0)/m2)/(RSS0_0/n2);    % causality y->x
prb_yx= 1 - cca_cdff(ftest_1,m2,n2);

%   do r-squared and check whiteness, consistency
% nobs-nlags=h;
% nvar*nlags=m;

    df_error = h-m;
    df_total = h;
   
        xvec = LY_1';
        rss2_x = xvec*xvec';
        rss(1) = 1 - (RSS0_0./ rss2_x);
        rss_adj(1) = 1 - ((RSS0_0/df_error) / (rss2_x/df_total) );
        waut(1,1) = sum(diff(u_x).^2)/sum(u_x.^2);
       
   
        yvec= LY1';
        rss2_y = yvec*yvec';
        rss(2) = 1 - (RSS0./ rss2_y);
        rss_adj(2) = 1 - ((RSS0/df_error) / (rss2_y/df_total) );
        waut(1,2) = sum(diff(u_y).^2)/sum(u_y.^2);
        

prb(2,1)=prb_xy;
prb(1,2)=prb_yx;
ret.prb=prb;

[PR,q] = cca_findsignificance(ret,PVAL,3);

Var(count,1) = GC.*PR(2,1);% G_xy_Ftest_ori ??GC ? q?????? GC?????????????? ???????
Var_yx0(count,1) = GC_1.*PR(1,2);

Temp_prb(i,1)=prb(2,1);
Temp_prb(i,2)=prb(1,2);


G_xy_Ftest_Ori(i,1)=GC;
G_xy_Ftest_Ori(i,2)=PR(2,1);
G_xy_Ftest_Ori(i,3)=q;
G_xy_Ftest_Ori(i,4)=prb(2,1);
G_xy_Ftest_Ori(i,5)=cond(LX1);
G_xy_Ftest_Ori(i,6)=cond(LX2);
if test_cond_num==1
    G_xy_Ftest_Ori(i,7)=GC_test;
    G_xy_Ftest_Ori(i,8)=Varxy;
    G_xy_Ftest_Ori(i,9)=Vary;
    G_xy_Ftest_Ori(i,10)=Varxy_test;
    G_xy_Ftest_Ori(i,11)=Vary_test;
end
    

G_yx_Ftest_Ori(i,1)=GC_1;
G_yx_Ftest_Ori(i,2)=PR(1,2);
G_yx_Ftest_Ori(i,3)=q;
G_yx_Ftest_Ori(i,4)=prb(1,2);
G_yx_Ftest_Ori(i,5)=cond(LX1);
G_yx_Ftest_Ori(i,6)=cond(LX2_1);
if test_cond_num==1
    G_yx_Ftest_Ori(i,7)=GC_1_test;
    G_yx_Ftest_Ori(i,8)=Varyx;
    G_yx_Ftest_Ori(i,9)=Varx;
    G_yx_Ftest_Ori(i,10)=Varyx_test;
    G_yx_Ftest_Ori(i,11)=Varx_test;
end


if (Var(count,1))<0
    LX1;
    LX2;
    %display(['The GC_xy is small Zero: ' num2str(1)]); 
end

if Var_yx0(count,1)<0
    LX1;
    LX2_1;
   % display(['The GC_yx is small Zero: ' num2str(1)]); 
end

clear LX2 LX1 D1 DATA1 LX2_1

else
    
% display('Now is The GC_xy: ')
beta_xy=LX1\LY1;
xpred_xy = LX1*beta_xy;
u_xy = LY1-xpred_xy;
Varxy=var(u_xy,1);
RSS0=u_xy'*u_xy;

beta_y=LX2\LY2;
xpred_y = LX2*beta_y;
u_y = LY2-xpred_y;
Vary=var(u_y,1);
RSS1=u_y'*u_y;
if Vary==0
GC=1;
else
GC=1-(Varxy/Vary);
end

if test_cond_num==1
[Varxy_test,pxy,Flagxy,RSS0_test]=LSE_New(LX1,LY1,h,m,number,1,0); % non-restricted regressions
[Vary_test,py,Flagy,RSS1_test]=LSE_New(LX2,LY2,h,m2,number,1,0); % restricted regressions

if Vary_test==0
GC_test=1;
else
GC_test=1-(Varxy_test/Vary_test);
end
end

% display('Now is The GC_yx: ')

beta_yx=LX1\LY_1;
xpred_yx = LX1*beta_yx;
u_yx = LY_1-xpred_yx;
Varyx=var(u_yx,1);
RSS0_0=u_yx'*u_yx;

beta_x=LX2_1\LY_2;
xpred_x = LX2_1*beta_x;
u_x = LY_2-xpred_x;
Varx=var(u_x,1);
RSS1_1=u_x'*u_x;
if Varx==0
GC_1=1;
else
GC_1=1-(Varyx/Varx);
end
if test_cond_num==1
[Varyx_test,pyx,Flagyx,RSS0_0_test]=LSE_New(LX1,LY_1,h,m,number,1,0); 
[Varx_test,px,Flagx,RSS1_1_test]=LSE_New(LX2_1,LY_2,h,m1,number,1,0); 
if Varx_test==0
GC_1_test=1;
else
GC_1_test=1-(Varyx_test/Varx_test);
end
end

chose_Two_Test_Method=2;
if chose_Two_Test_Method==2
rand_gc_point=zeros(h,m+1);
end
% Test :
    for trails_i=1:trails
     if chose_Two_Test_Method==1
       for num_point_region=1:h
        rand_gc_point=randperm(m+1) ; 
        LX1_Temp=[LX1(num_point_region,1:m) LY1(num_point_region)];
        
        LX11(num_point_region,1:m)=LX1_Temp(rand_gc_point(1:m));
        LY11(num_point_region,1)=LX1_Temp(rand_gc_point(end));
        LY21(num_point_region,1)=LX1_Temp(rand_gc_point(end));
        LX21(num_point_region,1:m2)=LX1_Temp(rand_gc_point(m1+1:m));
        
        LX2_Temp=[LX1(num_point_region,1:m) LY_1(num_point_region)];
        LX11_2(num_point_region,1:m)=LX2_Temp(rand_gc_point(1:m));
        LY_12(num_point_region,1)=LX2_Temp(rand_gc_point(end));
        LY_22(num_point_region,1)=LX2_Temp(rand_gc_point(end));
        LX2_12(num_point_region,1:m1)=LX2_Temp(rand_gc_point(1:m1));
        
       end
       
        elseif chose_Two_Test_Method==2
            
        for num_point_region=1:h 
        rand_gc_point(num_point_region,:)=(num_point_region-1)*(m+1)+randperm(m+1) ; 
        end
        LX1_Temp0=[LX1(:,1:m) LY1(:)]';
        LX1_Temp=LX1_Temp0(rand_gc_point);
        LX11=LX1_Temp(:,1:m);
        LY11=LX1_Temp(:,end);
        LY21=LX1_Temp(:,end);
        LX21=LX1_Temp(:,m1+1:m);
        
        LX2_Temp0=[LX1(:,1:m) LY_1(:)]';
        LX2_Temp=LX2_Temp0(rand_gc_point);
        LX11_2=LX2_Temp(:,1:m);
        LY_12=LX2_Temp(:,end);
        LY_22=LX2_Temp(:,end);
        LX2_12=LX2_Temp(:,1:m1);
     end
        
        % display('Now is The GC_xy: ')
        [Varxy,pxy,Flagxy,RSS0]=LSE_New(LX11,LY11,h,m,number,0,0); 
        [Vary,py,Flagy,RSS1]=LSE_New(LX21,LY21,h,m2,number,0,0); 


        if Vary==0
            GC_temp=1;
        else
            GC_temp=1-(Varxy/Vary);
        end

        % display('Now is The GC_yx: ')
            [Varyx,pyx,Flagyx,RSS0_0]=LSE_New(LX11_2,LY_12,h,m,number,0,0); 
            [Varx,px,Flagx,RSS1_1]=LSE_New(LX2_12,LY_22,h,m1,number,0,0); 

            if Varx==0
                GC_1_temp=1;
            else
                GC_1_temp=1-(Varyx/Varx);
            end   
            
            GC_xy_Distri(1,trails_i)=GC_temp;
            GC_yx_Distri(1,trails_i)=GC_1_temp;
        
    end
    GC_xy_Distri_Sort=sort(GC_xy_Distri,2);
    GC_yx_Distri_Sort=sort(GC_yx_Distri,2);
    
    
 fid=fopen('HH_28401/XYY_Threshold/HH_28401_GC_Threshold.txt','a+');
 fprintf(fid,'%s \t','Cycle');
 fprintf(fid,'%s \n','Trails');
 fprintf(fid,'%d \t',i);
 fprintf(fid,'%d \n',trails);
 
 fprintf(fid,'%s \t','GC_xy_Distri_Sort');
 for GC_xy_Distri_Sort_i=1:length(GC_xy_Distri_Sort)
 fprintf(fid,'%d \t',GC_xy_Distri_Sort(GC_xy_Distri_Sort_i));
     if GC_xy_Distri_Sort_i==length(GC_xy_Distri_Sort)
          fprintf(fid,'%d \n',GC_xy_Distri_Sort(GC_xy_Distri_Sort_i));
     end
 end
 fprintf(fid,'%s \t','GC_yx_Distri_Sort');
  for GC_yx_Distri_Sort_i=1:length(GC_yx_Distri_Sort)
 fprintf(fid,'%d \t',GC_yx_Distri_Sort(GC_yx_Distri_Sort_i));
     if GC_yx_Distri_Sort_i==length(GC_yx_Distri_Sort)
          fprintf(fid,'%d \n',GC_yx_Distri_Sort(GC_yx_Distri_Sort_i));
     end
  end
  fprintf(fid,'%s \n','-------------------------------------------------------------');
 fclose(fid);
    
    
    Length_GC=size(GC_xy_Distri,1);
    for Length_GC_i=1:Length_GC

        GC_xy_Threshold(Length_GC_i,i)=GC_xy_Distri_Sort(round(0.99*trails));
        GC_yx_Threshold(Length_GC_i,i)=GC_yx_Distri_Sort(round(0.99*trails));
    end
    
%     if (isnan(GC_xy_Threshold(Length_GC_i,i)))
%         
%         disp('Problem????')
%         
%     end
    Var(count,1) = GC*(GC>mean(GC_xy_Threshold(:,i)));% G_xy_Ftest_ori 
    Var_yx0(count,1) = GC_1*(GC_1>mean(GC_yx_Threshold(:,i)));
       
    G_xy_Ftest_Ori(i,1)=GC;
    G_xy_Ftest_Ori(i,2)=GC>mean(GC_xy_Threshold(:,i));
    if test_cond_num==1
    G_xy_Ftest_Ori(i,3)=GC_test;
    G_xy_Ftest_Ori(i,4)=cond(LX1'*LX1);
    G_xy_Ftest_Ori(i,5)=cond(LX2'*LX2);
    G_xy_Ftest_Ori(i,6)=Varxy;
    G_xy_Ftest_Ori(i,7)=Vary;
    G_xy_Ftest_Ori(i,8)=Varxy_test;
    G_xy_Ftest_Ori(i,9)=Vary_test;
    end

    G_yx_Ftest_Ori(i,1)=GC_1;
    G_yx_Ftest_Ori(i,2)=GC_1>mean(GC_yx_Threshold(:,i));
    if test_cond_num==1
    G_yx_Ftest_Ori(i,3)=GC_1_test;
    G_yx_Ftest_Ori(i,4)=cond(LX1'*LX1);
    G_yx_Ftest_Ori(i,5)=cond(LX2_1'*LX2_1);
    G_yx_Ftest_Ori(i,6)=Varyx;
    G_yx_Ftest_Ori(i,7)=Varx;
    G_yx_Ftest_Ori(i,8)=Varyx_test;
    G_yx_Ftest_Ori(i,9)=Varx_test;
    end
    
    
end
else
    Z_0_Reserve(1:length(pointset(:,referenceset(i))),i)=pointset(:,referenceset(i));
    Z_0_Reserve(1+length(pointset(:,referenceset(i))),i)=nan;
    Z_0_Reserve(2+length(pointset(:,referenceset(i))),i)=referenceset(i);
end
end
Var_0=Var(1:count,1);
% INDEX1=Var_0>=0;
% INDEX2=Var_0<=1;
% INDEX=INDEX1.*INDEX2;
% Var_xy=Var_0.*INDEX;
Var_xy=Var_0(find(Var_0>=0 & Var_0<=1));
Var_00=Var_yx0(1:count,1);
% INDEX10=Var_00>=0;
% INDEX20=Var_00<=1;
% INDEX0=INDEX10.*INDEX20;
% Var_yx=Var_00.*INDEX0;
Var_yx=Var_00(find(Var_0>=0 & Var_0<=1));

Prb_xy=mean(Temp_prb(:,1));
Prb_yx=mean(Temp_prb(:,2));
global_GC_Threshold=mean(G_yx_Ftest_Ori(:,3));
global_prb_gc=[Prb_xy Prb_yx];
%% Reserve the part of results
outfile=num2str(detal);
outfile2=num2str(Ref_Number);
outfile3=num2str(Couple_ie);

OutFile=[Resultsmat outfile outfile2 outfile3 'G_yx_Ftest_Ori' '.mat'];
OutFile1=[Resultsmat outfile outfile2 outfile3 'G_xy_Ftest_Ori' '.mat'];
OutFile2=[Resultsmat outfile outfile2 outfile3 'Z_0_Reserve' '.mat'];
OutFile3=[Resultsmat outfile outfile2 outfile3 'count_r' '.mat'];
OutFile4=[Resultsmat outfile outfile2 outfile3 'Var_xy' '.mat'];
OutFile5=[Resultsmat outfile outfile2 outfile3 'Var_yx' '.mat'];
OutFile6=[Resultsmat outfile outfile2 outfile3 'GC_xy_Threshold' '.mat'];
OutFile7=[Resultsmat outfile outfile2 outfile3 'GC_yx_Threshold' '.mat'];

save(OutFile,'G_yx_Ftest_Ori');
save(OutFile1,'G_xy_Ftest_Ori');
save(OutFile2,'Z_0_Reserve');
save(OutFile3,'count_r');
save(OutFile4,'Var_xy');
save(OutFile5,'Var_yx');
save(OutFile6,'GC_xy_Threshold');
save(OutFile7,'GC_yx_Threshold');

hh=length(Var_xy);
% saveas(gcf,['Overlap the attractor, neighborhood is  ' , num2str(detal) ],'png');
end
