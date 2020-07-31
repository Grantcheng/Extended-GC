% CEGCI for f1=0.280 f2=0.282 f3=0.290 Amp1=Amp2=10; % x->y,x->z
clc
clear
% close all
tic;
display('The HH Three Neurons All Couple Program ......   ');
global Threshold_Spike_x;
global Threshold_Spike_y;
global Threshold_Spike_z;
global Ref_Number;
global test_cond_num;
global Couple_ie;
global Reference_Points_Choose;
global Ration_S_NS;
global global_prb_gc;
global global_GC_Threshold;
global PVAL;

addpath utilities/;
settspath;
dim_hope=30;
information_dimension=0;
% N = input('input N  ');
% DetalN=input('input DetalN  ');
epstot=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16]
Dynamic_Data=100000;
time=0;
GC_Matrix_xy=zeros(1,2);
GC_Matrix_xz=zeros(1,2);
GC_Matrix_yz=zeros(1,2);
PVAL=input('The p-value of the GC test is:  ');
choose_same_tau=input('Whether use same tau cases:  Yes =  1,  No = 0  ');
data_chose=280;
Choose_Wolf_No_Sampling=0;
Chose_One_Couple=input('Whether use all ie  Yes =  1,  No = 0  ');
Noise_Add=input('Whether Noise_Add: Yes =1 , No = 0 ');
if Noise_Add==1
Realization=input('Input the realization .... ');
else
 Realization=1;   
end
neuro_type_choose=input('Which neuron type case choose: Case_1 = 1, Case_2 = 2 ');
test_cond_num=input('The test_cond_num is : Yes = 1, No = 0 ');
Reference_Points_Choose=input('Reference_Points_Choose is : 0 - Random;  1 - Threshold_Spike_x;  2 - Step_Fix;  3 - Threshold_Spike_x and Threshold_Spike_y  ');

if Noise_Add==1
    Noise_Type_Choose=input('The Noise_Type_Choose : 1 - SNR; 2 - DirectNoise   ');
if Noise_Type_Choose==1 
   SNR_Determined=input('The SNR is Determined :  ');
if SNR_Determined==1
    x_sigma=input('The x_sigma is :  ');
    y_sigma=input('The y_sigma is :  ');
    z_sigma=input('The z_sigma is :  ');
elseif SNR_Determined==0
sigma_white=input('The variance of white noise is :  ');
end
elseif Noise_Type_Choose==2
    sigma_white=input('The variance of white noise x is :  ');
    sigma_white_Noise1=input('The variance of white noise y is :  ');
    sigma_white_Noise2=input('The variance of white noise z is :  ');
end
end

if (Reference_Points_Choose==1 || Reference_Points_Choose==3)
Ration_S_NS=input('Ration_S_NS =  ');
end

if Chose_One_Couple==0
ie_3=input('The Couple Index is = ');
Leng_IE=length(ie_3);
Leng_IE_total=ie_3;
elseif Chose_One_Couple==1
Leng_IE=length(epstot);
    ie_3=input('The Start of Couple Index is = ');
    Leng_IE_total=Leng_IE;
end

    mx_record=zeros(Leng_IE,Realization);
    my_record=zeros(Leng_IE,Realization);
    mz_record=zeros(Leng_IE,Realization);
    taux_record=zeros(Leng_IE,Realization);
    tauy_record=zeros(Leng_IE,Realization);
    tauz_record=zeros(Leng_IE,Realization);
    
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
Var_xz=zeros(length(M),Realization);RSS_zx=zeros(length(M),Realization);
Var_yz=zeros(length(M),Realization);RSS_zy=zeros(length(M),Realization);
Preprocessing=0;

ie_1=input('The mx and my reset = ');
if ie_1==0
mx_vect=9;
my_vect=7;
mz_vect=7;
Tau_xy=29;

CGC_xy=zeros(Leng_IE,Realization);%
CGC_yx=zeros(Leng_IE,Realization);%

CGC_xz=zeros(Leng_IE,Realization);%
CGC_zx=zeros(Leng_IE,Realization);%

CGC_yz=zeros(Leng_IE,Realization);%
CGC_zy=zeros(Leng_IE,Realization);%

elseif ie_1==1
mx_vect=input('The mx is = ');
my_vect=input('The my is = ');
mz_vect=input('The mz is = ');
Tau_xy=input('Tau_xy is = ');

CGC_xy=zeros(Leng_IE,Realization);%
CGC_yx=zeros(Leng_IE,Realization);%

CGC_xz=zeros(Leng_IE,Realization);%
CGC_zx=zeros(Leng_IE,Realization);%

CGC_yz=zeros(Leng_IE,Realization);%
CGC_zy=zeros(Leng_IE,Realization);%

elseif ie_1==2
mx_vect=[6 7 8 9 10 11 12 13];
my_vect=[6 7 8 9 10 11 12 13];
mz_vect=[6 7 8 9 10 11 12 13];
Tau_xy=29;
Leng_vector_xyz=length(mx_vect)*length(my_vect)*length(mz_vect);

m_record=zeros(Leng_vector_xyz,3);
CGC_xy=zeros(Leng_vector_xyz,Realization);%
CGC_yx=zeros(Leng_vector_xyz,Realization);%

CGC_xz=zeros(Leng_vector_xyz,Realization);%
CGC_zx=zeros(Leng_vector_xyz,Realization);%

CGC_yz=zeros(Leng_vector_xyz,Realization);%
CGC_zy=zeros(Leng_vector_xyz,Realization);%

elseif ie_1==3
    mxmymz_reset=input('Whether mx my mz independent process = ');
    if mxmymz_reset==0
mx_vect=[6 7 8 9 10 11 12 13];
my_vect=[6 7 8 9 10 11 12 13];
mz_vect=[6 7 8 9 10 11 12 13];
    elseif mxmymz_reset==1
     mxmymz_reset_i=input('Choose the mx my mz reset process = ');   
        if mxmymz_reset_i==1
          mx_vect=[6 7];
          my_vect=[6 7 8 9 10 11 12 13];
          mz_vect=[6 7 8 9 10 11 12 13];  
        elseif mxmymz_reset_i==2
          mx_vect=[8 9];
          my_vect=[6 7 8 9 10 11 12 13];
          mz_vect=[6 7 8 9 10 11 12 13];           
        elseif mxmymz_reset_i==3
          mx_vect=[10 11];
          my_vect=[6 7 8 9 10 11 12 13];
          mz_vect=[6 7 8 9 10 11 12 13]; 
        elseif mxmymz_reset_i==4
          mx_vect=[12 13];
          my_vect=[6 7 8 9 10 11 12 13];
          mz_vect=[6 7 8 9 10 11 12 13];    
        end
    end
    
Leng_vector_xyz=length(mx_vect)*length(my_vect)*length(mz_vect);

m_record=zeros(Leng_vector_xyz,3);
CGC_xy=zeros(Leng_vector_xyz,Realization);%
CGC_yx=zeros(Leng_vector_xyz,Realization);%

CGC_xz=zeros(Leng_vector_xyz,Realization);%
CGC_zx=zeros(Leng_vector_xyz,Realization);%

CGC_yz=zeros(Leng_vector_xyz,Realization);%
CGC_zy=zeros(Leng_vector_xyz,Realization);%
end
    
Auto_Tau=input('Auto_Tau is = ');
Ref_Number=input('Ref_Number is = ');
display('The New Code of three neurons network In EGCI ');
display(['The Chose_One_Couple is: ' num2str(Chose_One_Couple)]);
display(['The data_chose is: ' num2str(data_chose)]);
Taken_dim=zeros(Leng_IE,1);

for ie=ie_3:Leng_IE_total
    mx=0;
    my=0;
    mz=0;
    taux=0;
    tauy=0;
    tauz=0;
    
    C=epstot(ie);
    Couple_ie=C;
    display(['The Strength is : ' num2str(C)]);
  if neuro_type_choose==1      
   if ie==1
    V1=load('NeuroType_Case/Case_1/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0/eig_21_10000000.txt');
    
   elseif ie==2
    V1=load('NeuroType_Case/Case_1/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.02/eig_21_10000000.txt');
    
    elseif ie==3
    V1=load('NeuroType_Case/Case_1/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.04/eig_21_10000000.txt');
    
    elseif ie==4
    V1=load('NeuroType_Case/Case_1/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.06/eig_21_10000000.txt');
    
    elseif ie==5
    V1=load('NeuroType_Case/Case_1/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.08/eig_21_10000000.txt');
    
    elseif ie==6
    V1=load('NeuroType_Case/Case_1/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.10/eig_21_10000000.txt');
     
    elseif ie==7
    V1=load('NeuroType_Case/Case_1/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.11/eig_21_10000000.txt');
     
    elseif ie==8
    V1=load('NeuroType_Case/Case_1/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.12/eig_21_10000000.txt');
        
    elseif ie==9
    V1=load('NeuroType_Case/Case_1/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.13/eig_21_10000000.txt');
    
    elseif ie==10
    V1=load('NeuroType_Case/Case_1/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.14/eig_21_10000000.txt');
        
    elseif ie==11
    V1=load('NeuroType_Case/Case_1/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.15/eig_21_10000000.txt');
     
    elseif ie==12
    V1=load('NeuroType_Case/Case_1/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('NeuroType_Case/Case_1/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('NeuroType_Case/Case_1/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    eigvalue=load('NeuroType_Case/Case_1/Couple=0.16/eig_21_10000000.txt');
   end
  elseif neuro_type_choose==2
      if ie==1
    V1=load('NeuroType_Case/Case_2/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0/eig_21_10000000.txt');
    
   elseif ie==2
    V1=load('NeuroType_Case/Case_2/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.02/eig_21_10000000.txt');
    
    elseif ie==3
    V1=load('NeuroType_Case/Case_2/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.04/eig_21_10000000.txt');
    
    elseif ie==4
    V1=load('NeuroType_Case/Case_2/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.06/eig_21_10000000.txt');
    
    elseif ie==5
    V1=load('NeuroType_Case/Case_2/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.08/eig_21_10000000.txt');
    
    elseif ie==6
    V1=load('NeuroType_Case/Case_2/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.10/eig_21_10000000.txt');
     
    elseif ie==7
    V1=load('NeuroType_Case/Case_2/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.11/eig_21_10000000.txt');
     
    elseif ie==8
    V1=load('NeuroType_Case/Case_2/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.12/eig_21_10000000.txt');
        
    elseif ie==9
    V1=load('NeuroType_Case/Case_2/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.13/eig_21_10000000.txt');
    
    elseif ie==10
    V1=load('NeuroType_Case/Case_2/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.14/eig_21_10000000.txt');
        
    elseif ie==11
    V1=load('NeuroType_Case/Case_2/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.15/eig_21_10000000.txt');
     
    elseif ie==12
    V1=load('NeuroType_Case/Case_2/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('NeuroType_Case/Case_2/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('NeuroType_Case/Case_2/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    eigvalue=load('NeuroType_Case/Case_2/Couple=0.16/eig_21_10000000.txt');
      end
  end
%% compute the information dimension and Taken's theorem
sum_eig=0;sum_eig_back=0;
for eig_i=1:length(eigvalue)-1
    sum_eig_back=sum_eig;
    sum_eig=sum_eig+eigvalue(eig_i);
    if (sum_eig_back>0)&&(sum_eig<0)
        information_dimension=(eig_i-1)+(sum_eig_back/abs(eigvalue(eig_i)));
    end
end   
Taken_dim(ie,1)=2*information_dimension+1;
%%------------------------------------
Sampling=0.01;
stv=0.01;% the sample interval
dt=0.01;% the data interval
stv_L=stv/dt;
N1=length(V1(1,:));
len_data=N1-Dynamic_Data;
X=[V1(1,Dynamic_Data:end);V2(1,Dynamic_Data:end);V3(1,Dynamic_Data:end)];
N=length(X(1,:));

display(['The N is: ' num2str(N)]);
clear V1;
clear V2;
clear V3;
x = X(1,:);
y = X(2,:);
z = X(3,:);

Threshold_Spike_x=(-55-min(x))./(max(x)-min(x));
Threshold_Spike_y=(-55-min(y))./(max(y)-min(y));
Threshold_Spike_z=(-55-min(z))./(max(z)-min(z));

x=(x-min(x))./(max(x)-min(x));
y=(y-min(y))./(max(y)-min(y));
z=(z-min(z))./(max(z)-min(z));
X=[x;y;z];

display(['The Auto_Tau is: ' num2str(Auto_Tau)]);
Leng_vector_xyz_i=0;
for mx_i=1:length(mx_vect)
    for my_i=1:length(my_vect)
        for mz_i=1:length(mz_vect)
            Leng_vector_xyz_i=Leng_vector_xyz_i+1;
for r_i=1:Realization
 x=X(1,:);
 y=X(2,:);
 z=X(3,:);
max_k=length(x);
if Noise_Add==1
    if Noise_Type_Choose==1
if  SNR_Determined==1
 sigma_white=sqrt(var(x))/10^(x_sigma)
 sigma_white_Noise1=sqrt(var(y))/10^(y_sigma)
 sigma_white_Noise2=sqrt(var(z))/10^(z_sigma)
end
White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white_Noise1.*randn(1,max_k);
White_Noise2=sigma_white_Noise2.*randn(1,max_k);
    elseif Noise_Type_Choose==2
White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white_Noise1.*randn(1,max_k);
White_Noise2=sigma_white_Noise2.*randn(1,max_k);       
    end
x=White_Noise+x;
y=White_Noise1+y;
z=White_Noise2+z;
end
%% ------------------len=length(X(1,:));
mm=0;
    %% -------Set the embedding dimension and time delay
  if Auto_Tau==1
     taux=0;
     tauy=0;
     tauz=0;
s1=signal(x(1,1:len_data)');
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

s2=signal(y(1,1:len_data)');
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

s3=signal(z(1,1:len_data)');
a3=amutual2(s3,250);
a3_1=data(a3);
for i=1:length(a3_1)-2
if (a3_1(i+1)<a3_1(i)) && (a3_1(i+1)<a3_1(i+2))
tauz=i
break;
end
end
if tauz==0
a3=amutual2(s3,350);
a3_1=data(a3);
for i=1:length(a3_1)-2
if (a3_1(i+1)<a3_1(i)) && (a3_1(i+1)<a3_1(i+2))
tauz=i
break;
end
end    
end

c1 = cao(s1,21,taux,3,1000);
c2 = cao(s2,21,tauy,3,1000);
c3 = cao(s3,21,tauz,3,1000);

c11=data(c1);
c21=data(c2);
c31=data(c3);

for i=1:length(c11)-1
EE(i,1)=var(c11(i:end));
end
for i=1:length(c21)-1
EE2(i,1)=var(c21(i:end));
end
for i=1:length(c31)-1
EE3(i,1)=var(c31(i:end));
end

HH=diff(EE);
HH2=diff(EE2);
HH3=diff(EE3);

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

for i=1:length(HH3)
   if abs(HH3(i))<10^(-4)
       mz=i+1
       break;
   end
end

clear EE;
clear EE2;
clear EE3;

if mx==0
c1 = cao(s1,dim_hope,taux,3,1000);
c11=data(c1);
for i=1:length(c11)-1
EE(i,1)=var(c11(i:end));
end
HH=diff(EE);
for i=1:length(HH)
   if abs(HH(i))<10^(-4)
       mx=i+1
       break;
   end
end
end

if my==0
c2 = cao(s2,dim_hope,tauy,3,1000);
c21=data(c2);
for i=1:length(c21)-1
EE2(i,1)=var(c21(i:end));
end
HH2=diff(EE2);
for i=1:length(HH2)
   if abs(HH2(i))<10^(-4)
       my=i+1
       break;
   end
end
end

if mz==0
c3 = cao(s3,dim_hope,tauz,3,1000);
c31=data(c3);
for i=1:length(c31)-1
EE3(i,1)=var(c31(i:end));
end
HH3=diff(EE3);
for i=1:length(HH3)
   if abs(HH3(i))<10^(-4)
       mz=i+1
       break;
   end
end
end
clear s1;
clear s2;
clear s3;
clear c1;
clear c2;
clear c3;
if mx<Taken_dim(ie,1)
mx=1+ceil(Taken_dim(ie,1));
end

if my<Taken_dim(ie,1)
my=1+ceil(Taken_dim(ie,1));
end

if mz<Taken_dim(ie,1)
mz=1+ceil(Taken_dim(ie,1));
end
 elseif Auto_Tau==2
    taux=Tau_xy
    mx=mx_vect(mx_i)
    my=my_vect(my_i)
    mz=mz_vect(mz_i)
    tauy=taux % 41
    tauz=taux % 41
 elseif Auto_Tau==3
     mx=mx_vect(mx_i)
    my=my_vect(my_i)
    mz=mz_vect(mz_i)
     taux=0;
     tauy=0;
     tauz=0;
s1=signal(x(1,1:len_data)');
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

s2=signal(y(1,1:len_data)');
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

s3=signal(z(1,1:len_data)');
a3=amutual2(s3,250);
a3_1=data(a3);
for i=1:length(a3_1)-2
if (a3_1(i+1)<a3_1(i)) && (a3_1(i+1)<a3_1(i+2))
tauz=i
break;
end
end
if tauz==0
a3=amutual2(s3,350);
a3_1=data(a3);
for i=1:length(a3_1)-2
if (a3_1(i+1)<a3_1(i)) && (a3_1(i+1)<a3_1(i+2))
tauz=i
break;
end
end    
end
 end
m=mx+my+mz;
if choose_same_tau==1
tau_final=min([taux,tauy,tauz]);
taux=tau_final;
tauy=tau_final;
tauz=tau_final;
end
%--------------------------------------------------
% (phase space reconstruction)
[xn,L1] = PhaSpaRecon(x,taux,mx);
[yn,L2] = PhaSpaRecon(y,tauy,my);
[zn,L3] = PhaSpaRecon(z,tauz,mz);
%%
global Resultsmat;
if data_chose==280
     if Noise_Add==1
         if neuro_type_choose==1     
         if choose_same_tau==1
         str_name=(['NeuroType_Case/Case_1/Noise_Results_All/mx_',num2str(mx),' my_',num2str(my),' mz_',num2str(mz),'Sigma_',num2str(sigma_white),' Results mat']);  
         if Noise_Type_Choose==1
         str_name_1=(['NeuroType_Case/Case_1/Noise_Results_All/SNR_X=',num2str(x_sigma),' SNR_Y=',num2str(y_sigma),' SNR_Z=',num2str(z_sigma),'_']);  
         elseif Noise_Type_Choose==2
             str_name_1=(['NeuroType_Case/Case_1/Noise_Results_All/Noise_X=',num2str(sigma_white),' Noise_Y=',num2str(sigma_white_Noise1),' Noise_Z=',num2str(sigma_white_Noise1),'_']); 
         end
         else
         str_name=(['NeuroType_Case/Case_1/Noise_Results_All/mx_',num2str(mx),' my_',num2str(my),' mz_',num2str(mz),'Sigma_',num2str(sigma_white),' Results mat']);  
         if Noise_Type_Choose==1
         str_name_1=(['NeuroType_Case/Case_1/Noise_Results_All/SNR_X=',num2str(x_sigma),' SNR_Y=',num2str(y_sigma),' SNR_Z=',num2str(z_sigma),'Difftau_']);  
         elseif Noise_Type_Choose==2
         str_name_1=(['NeuroType_Case/Case_1/Noise_Results_All/Noise_X=',num2str(sigma_white),' Noise_Y=',num2str(sigma_white_Noise1),' Noise_Z=',num2str(sigma_white_Noise1),'Difftau_']); 
         end
         end
         elseif neuro_type_choose==2
          if choose_same_tau==1
         str_name=(['NeuroType_Case/Case_2/Noise_Results_All/mx_',num2str(mx),' my_',num2str(my),' mz_',num2str(mz),'Sigma_',num2str(sigma_white),' Results mat']);  
         if Noise_Type_Choose==1
         str_name_1=(['NeuroType_Case/Case_2/Noise_Results_All/SNR_X=',num2str(x_sigma),' SNR_Y=',num2str(y_sigma),' SNR_Z=',num2str(z_sigma),'_']);  
         elseif Noise_Type_Choose==2
             str_name_1=(['NeuroType_Case/Case_2/Noise_Results_All/Noise_X=',num2str(sigma_white),' Noise_Y=',num2str(sigma_white_Noise1),' Noise_Z=',num2str(sigma_white_Noise1),'_']); 
         end
         else
         str_name=(['NeuroType_Case/Case_2/Noise_Results_All/mx_',num2str(mx),' my_',num2str(my),' mz_',num2str(mz),'Sigma_',num2str(sigma_white),' Results mat']);  
         if Noise_Type_Choose==1
         str_name_1=(['NeuroType_Case/Case_2/Noise_Results_All/SNR_X=',num2str(x_sigma),' SNR_Y=',num2str(y_sigma),' SNR_Z=',num2str(z_sigma),'Difftau_']);  
         elseif Noise_Type_Choose==2
         str_name_1=(['NeuroType_Case/Case_2/Noise_Results_All/Noise_X=',num2str(sigma_white),' Noise_Y=',num2str(sigma_white_Noise1),' Noise_Z=',num2str(sigma_white_Noise1),'Difftau_']); 
         end
         end            
         end
     else
         if neuro_type_choose==1 
         str_name=(['NeuroType_Case/Case_1/Results_All/mx_',num2str(mx),' my_',num2str(my),' mz_',num2str(mz),' Results mat']);
         elseif neuro_type_choose==2
         str_name=(['NeuroType_Case/Case_2/Results_All/mx_',num2str(mx),' my_',num2str(my),' mz_',num2str(mz),' Results mat']);    
         end
     end
if Noise_Add==1
Resultsmat=[pwd,'/',str_name];
Resultsmat_1=[pwd,'/',str_name_1];
else
Resultsmat=[pwd,'/',str_name];  
end
end

%% Reconstruct phase attractor
L=min([L1,L2,L3]);
%% ------------------------------------------ 
ZX_xy=zeros(mx+my+mz,L);
ZX_xz=zeros(mx+my+mz,L);
ZX_yz=zeros(mx+my+mz,L);
for i=1:L
   ZX_xy(:,i)=[xn(:,i);yn(:,i);zn(:,i)];%
   ZX_xz(:,i)=[xn(:,i);zn(:,i);yn(:,i)];%
   ZX_yz(:,i)=[yn(:,i);zn(:,i);xn(:,i)];%
end
clear xn;
clear yn;
clear zn;
%% -------------------------------------------------------
for u=1:length(M)
    display(['The Detal is: ' num2str(M(u))]);
%-------------------------------------------------------
[count_r,neighbors]=CGC_Three_Neighbors(ZX_xy,M(u),L);
%-------------------------------------------------------
    tic;
    [Varxy,Varyx,Lxy]=Conditional_Granger_Causality_Three_Neurons(ZX_xy,M(u),max_k,mx,my,mz,taux,tauy,1,count_r,neighbors);%x->y,y->x conditional on z
    global_prb_gc_xy=global_prb_gc;
    global_GC_Threshold_xy=global_GC_Threshold;
    global_prb_gc=nan;
    global_GC_Threshold=nan;
    
    [Varxz,Varzx,Lxz]=Conditional_Granger_Causality_Three_Neurons(ZX_xz,M(u),max_k,mx,mz,my,taux,tauz,1,count_r,neighbors);%x->z,z->x conditional on y
    global_prb_gc_xz=global_prb_gc;
    global_GC_Threshold_xz=global_GC_Threshold;
    global_prb_gc=nan;
    global_GC_Threshold=nan;
    
    [Varyz,Varzy,Lyz]=Conditional_Granger_Causality_Three_Neurons(ZX_yz,M(u),max_k,my,mz,mx,tauy,tauz,1,count_r,neighbors);%y->z,z->y conditional on x
    global_prb_gc_yz=global_prb_gc;
    global_GC_Threshold_yz=global_GC_Threshold;
    global_prb_gc=nan;
    global_GC_Threshold=nan;
    
    Lyx=Lxy;
    Lzx=Lxz;
    Lzy=Lyz;
    toc;
    time=toc+time;

clear ZX_xz;
clear ZX_xy;
clear ZX_yz;
%% -----------------------------------------
Var2=Varxy(1:Lxy,1);%
RSS2=Varyx(1:Lyx,1);%

Var2_xz=Varxz(1:Lxz,1);%
RSS2_zx=Varzx(1:Lzx,1);%

Var2_yz=Varyz(1:Lyz,1);%
RSS2_zy=Varzy(1:Lzy,1);%

Var(u,r_i)=mean(Var2);%
RSS(u,r_i)=mean(RSS2);%

Var_xz(u,r_i)=mean(Var2_xz);%
RSS_zx(u,r_i)=mean(RSS2_zx);%

Var_yz(u,r_i)=mean(Var2_yz);%
RSS_zy(u,r_i)=mean(RSS2_zy);%

GC_Matrix_xy(global_prb_gc_xy<global_GC_Threshold_xy) = 1;
GC_Matrix_xz(global_prb_gc_xz<global_GC_Threshold_xz) = 1;
GC_Matrix_yz(global_prb_gc_yz<global_GC_Threshold_yz) = 1;

Var(u,r_i)=GC_Matrix_xy(1,1)*Var(u,r_i);%
RSS(u,r_i)=GC_Matrix_xy(1,2)*RSS(u,r_i);%

Var_xz(u,r_i)=GC_Matrix_xz(1,1)*Var_xz(u,r_i);%
RSS_zx(u,r_i)=GC_Matrix_xz(1,2)*RSS_zx(u,r_i);%

Var_yz(u,r_i)=GC_Matrix_yz(1,1)*Var_yz(u,r_i);%
RSS_zy(u,r_i)=GC_Matrix_yz(1,2)*RSS_zy(u,r_i);%
%% ----------------------------------------
display(['GC_xy = ',num2str(Var(u,r_i))]);
display(['GC_yx = ',num2str(RSS(u,r_i))]);

display(['GC_xz = ',num2str(Var_xz(u,r_i))]);
display(['GC_zx = ',num2str(RSS_zx(u,r_i))]);

display(['GC_yz = ',num2str(Var_yz(u,r_i))]);
display(['GC_zy = ',num2str(RSS_zy(u,r_i))]);

if Noise_Add==1
display(['SNR of X is = ',num2str(log10(sqrt(var(x))/sigma_white))]);
display(['SNR of Y is = ',num2str(log10(sqrt(var(y))/sigma_white_Noise1))]);
display(['SNR of Z is = ',num2str(log10(sqrt(var(z))/sigma_white_Noise2))]);
end
end
time=toc+time;
display('The total time cost is :');
time
mx_record(ie,r_i)=mx;
my_record(ie,r_i)=my;
mz_record(ie,r_i)=mz;
taux_record(ie,r_i)=taux;
tauy_record(ie,r_i)=tauy;
tauz_record(ie,r_i)=tauz;
end
%-------------------------------------------------
%  fid=fopen('HH_280/HH_280_EGCI_Results.txt','a+');
if data_chose==280
    if Noise_Add==1
        if neuro_type_choose==1     
        OutFile_xy=[Resultsmat_1 'HH_Three_NeuroType_Case_1_EGCI_Noise_Results_xy' '.txt'];
        OutFile_xz=[Resultsmat_1 'HH_Three_NeuroType_Case_1_EGCI_Noise_Results_xz' '.txt'];
        OutFile_yz=[Resultsmat_1 'HH_Three_NeuroType_Case_1_EGCI_Noise_Results_yz' '.txt'];
        elseif neuro_type_choose==2
        OutFile_xy=[Resultsmat_1 'HH_Three_NeuroType_Case_2_EGCI_Noise_Results_xy' '.txt'];
        OutFile_xz=[Resultsmat_1 'HH_Three_NeuroType_Case_2_EGCI_Noise_Results_xz' '.txt'];
        OutFile_yz=[Resultsmat_1 'HH_Three_NeuroType_Case_2_EGCI_Noise_Results_yz' '.txt'];            
        end
        fid=fopen(OutFile_xy,'a+');
        fid_xz=fopen(OutFile_xz,'a+');
        fid_yz=fopen(OutFile_yz,'a+');
    else
        if neuro_type_choose==1 
        fid=fopen('NeuroType_Case/Case_1/HH_Three_NeuroType_Case_1_EGCI_Results_xy.txt','a+');
        fid_xz=fopen('NeuroType_Case/Case_1/HH_Three_NeuroType_Case_1_EGCI_Results_xz.txt','a+');
        fid_yz=fopen('NeuroType_Case/Case_1/HH_Three_NeuroType_Case_1_EGCI_Results_yz.txt','a+'); 
      elseif neuro_type_choose==2
        fid=fopen('NeuroType_Case/Case_1/HH_Three_NeuroType_Case_2_EGCI_Results_xy.txt','a+');
        fid_xz=fopen('NeuroType_Case/Case_1/HH_Three_NeuroType_Case_2_EGCI_Results_xz.txt','a+');
        fid_yz=fopen('NeuroType_Case/Case_1/HH_Three_NeuroType_Case_2_EGCI_Results_yz.txt','a+');           
        end
    end
end
 fprintf(fid,'%s \t','Couple');
 fprintf(fid,'%s \t','mx');
 fprintf(fid,'%s \t','my');
 fprintf(fid,'%s \t','mz');
 fprintf(fid,'%s \t','T_x');
 fprintf(fid,'%s \t','T_y');
 fprintf(fid,'%s \n','T_z');
 fprintf(fid,'\r\n');
 fprintf(fid,'%d \t',C);
 fprintf(fid,'%d \t',mx);
 fprintf(fid,'%d \t',my);
 fprintf(fid,'%d \t',mz);
 fprintf(fid,'%d \t',taux);
 fprintf(fid,'%d \t',tauy);
 fprintf(fid,'%d \n',tauz);
 fprintf(fid,'\r\n');
  
 fprintf(fid_xz,'%s \t','Couple');
 fprintf(fid_xz,'%s \t','mx');
 fprintf(fid_xz,'%s \t','my');
 fprintf(fid_xz,'%s \t','mz');
 fprintf(fid_xz,'%s \t','T_x');
 fprintf(fid_xz,'%s \t','T_y');
 fprintf(fid_xz,'%s \n','T_z');
 fprintf(fid_xz,'\r\n');
 fprintf(fid_xz,'%d \t',C);
 fprintf(fid_xz,'%d \t',mx);
 fprintf(fid_xz,'%d \t',my);
 fprintf(fid_xz,'%d \t',mz);
 fprintf(fid_xz,'%d \t',taux);
 fprintf(fid_xz,'%d \t',tauy);
 fprintf(fid_xz,'%d \n',tauz);
 fprintf(fid_xz,'\r\n');
  
 fprintf(fid_yz,'%s \t','Couple');
 fprintf(fid_yz,'%s \t','mx');
 fprintf(fid_yz,'%s \t','my');
 fprintf(fid_yz,'%s \t','mz');
 fprintf(fid_yz,'%s \t','T_x');
 fprintf(fid_yz,'%s \t','T_y');
 fprintf(fid_yz,'%s \n','T_z');
 fprintf(fid_yz,'\r\n');
 fprintf(fid_yz,'%d \t',C);
 fprintf(fid_yz,'%d \t',mx);
 fprintf(fid_yz,'%d \t',my);
 fprintf(fid_yz,'%d \t',mz);
 fprintf(fid_yz,'%d \t',taux);
 fprintf(fid_yz,'%d \t',tauy);
 fprintf(fid_yz,'%d \n',tauz);
 fprintf(fid_yz,'\r\n'); 
%%--------------------------------------------------------------------xy 
 fprintf(fid,'%s \t','Var');
 for Var_ii=1:length(Var)
     if Var_ii==length(Var)
          fprintf(fid,'%d \n',Var(Var_ii));
     else
          fprintf(fid,'%d \t',Var(Var_ii));
     end
     
 end
 fprintf(fid,'\r\n');
 fprintf(fid,'%s \t','RSS');
  for Var_ii=1:length(RSS)
     if Var_ii==length(RSS)
          fprintf(fid,'%d \n',RSS(Var_ii));
     else
          fprintf(fid,'%d \t',RSS(Var_ii));
     end
  end
  fprintf(fid,'%s \n','---------------------------------------------');
  fprintf(fid,'\r\n');
 fclose(fid);
%%--------------------------------------------------------------------xz  
  fprintf(fid_xz,'%s \t','Var');
 for Var_ii=1:length(Var_xz)
     if Var_ii==length(Var_xz)
          fprintf(fid_xz,'%d \n',Var_xz(Var_ii));
     else
          fprintf(fid_xz,'%d \t',Var_xz(Var_ii));
     end
 end
 fprintf(fid_xz,'\r\n');
 fprintf(fid_xz,'%s \t','RSS');
  for Var_ii=1:length(RSS_zx)
      if Var_ii==length(RSS_zx)
          fprintf(fid_xz,'%d \n',RSS_zx(Var_ii));
     else
         fprintf(fid_xz,'%d \t',RSS_zx(Var_ii));
     end
  end
  fprintf(fid_xz,'%s \n','---------------------------------------');
  fprintf(fid_xz,'\r\n');
 fclose(fid_xz);
%%--------------------------------------------------------------------yz  
  fprintf(fid_yz,'%s \t','Var');
 for Var_ii=1:length(Var_yz)
     if Var_ii==length(Var_yz)
          fprintf(fid_yz,'%d \n',Var_yz(Var_ii));
     else
          fprintf(fid_yz,'%d \t',Var_yz(Var_ii));
     end
 end
 fprintf(fid_yz,'\r\n');
 fprintf(fid_yz,'%s \t','RSS');
  for Var_ii=1:length(RSS_zy)
     if Var_ii==length(RSS_zy)
          fprintf(fid_yz,'%d \n',RSS_zy(Var_ii));
     else
         fprintf(fid_yz,'%d \t',RSS_zy(Var_ii));
     end
  end
  fprintf(fid_yz,'%s \n','--------------------------------------------------');
  fprintf(fid_yz,'\r\n');
  fclose(fid_yz);
 
if Chose_One_Couple==1
CGC_xy(ie,1:Realization)=Var;%
CGC_yx(ie,1:Realization)=RSS;%

CGC_xz(ie,1:Realization)=Var_xz;%
CGC_zx(ie,1:Realization)=RSS_zx;%

CGC_yz(ie,1:Realization)=Var_yz;%
CGC_zy(ie,1:Realization)=RSS_zy;%

if Noise_Add==1
SNR(ie,1:Realization)=log10(sqrt(var(X(2,:)))/sigma_white);
end
Tau(ie,1:Realization)=taux;  
else
CGC_xy(Leng_vector_xyz_i,1:Realization)=Var;%
CGC_yx(Leng_vector_xyz_i,1:Realization)=RSS;%

CGC_xz(Leng_vector_xyz_i,1:Realization)=Var_xz;%
CGC_zx(Leng_vector_xyz_i,1:Realization)=RSS_zx;%

CGC_yz(Leng_vector_xyz_i,1:Realization)=Var_yz;%
CGC_zy(Leng_vector_xyz_i,1:Realization)=RSS_zy;%

m_record(Leng_vector_xyz_i,:)=[mx my mz];

if Noise_Add==1
SNR(1,1:Realization)=log10(sqrt(var(X(2,:)))/sigma_white);
end
Tau(Leng_vector_xyz_i,1:Realization)=taux; 
end
display(['The Strength is : ' num2str(C)]);
display('The next Couple');
        end
    end
end
end

display(['Noise_Add: ' num2str(Noise_Add)]);

display(['Noise_Add: ' num2str(Noise_Add)]);

