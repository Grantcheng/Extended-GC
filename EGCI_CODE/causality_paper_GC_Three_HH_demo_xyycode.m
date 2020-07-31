function causality_paper_GC_Three_HH_demo_xyycode(DEMO)
clear;
clc;
%-----------------------------------------------------------------------
%   demo parameters
PVAL    =   0.001;       % probability threshold for Granger causality significance
NLAGS   =   -1;         % if -1, best model order is assessed automatically

GC_CAL_HOME = fileparts(mfilename('fullpath'));

addpath(GC_CAL_HOME);
addpath([GC_CAL_HOME,'/GC_clean-mainline/GCcal']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/GCcal_spectrum']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/tools_and_utilities']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/prj_neuron_gc']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/prj_neuron_gc/scan_worker_template']);
addpath([GC_CAL_HOME,'/gc_time_frequency']);
addpath([GC_CAL_HOME,'/gc_time_frequency/utilities/']);

disp('GC calculation package path added.');
disp(' ');
%%
epstot=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];

%Cases_choose=input('Which case choose: Fig_3 = 1, Fig_4 = 2, Fig_5 = 3, Fig_11 = 4 ');

%% Noise
Noise_Add=input('Whether Noise_Add: Yes =1 , No = 0 ');
if Noise_Add==1
SNR_Determined=input('The SNR is Determined :  ');
if SNR_Determined==1
    x_sigma=input('The x_sigma is :  ');
    y_sigma=input('The y_sigma is :  ');
    z_sigma=input('The z_sigma is :  ');
elseif SNR_Determined==0
sigma_white=input('The variance of white noise is :  ');
end
end

%%
% acquire demo data
if nargin == 0,
    DEMO = 1;
end
Num=1000000;

for Cases_choose=1:5
GC_Values_Save = zeros(3,3,12);    
GC2_Values_Save = zeros(3,3,12);   
gc_zero_line_Save=zeros(1,12);
for ie=1:12
%% import data: Three Neurons Network Case
    C=epstot(ie);
    Couple_ie=C;
    display(['The Cases_choose is: ' num2str(Cases_choose)]);
    display(['The Strength is : ' num2str(C)]);
  if Cases_choose==1
   if ie==1
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
      
   elseif ie==2
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
     
    elseif ie==7
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
           
    elseif ie==11
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
        
    elseif ie==12
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.282w3_0.29.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.282w3_0.29.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig8_Case1/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.282w3_0.29.txt');
   end
  elseif Cases_choose==2
      if ie==1
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig6_Case2/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
      end
    elseif Cases_choose==3
       if ie==1
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig7_Case3/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
      end       
    elseif Cases_choose==4
    if ie==1
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig9_Case4/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    end      
    elseif Cases_choose==5
    if ie==1
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.02/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.02/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.02/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.04/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.04/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.04/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.06/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.06/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.06/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.08/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.08/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.08/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.10/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.10/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.10/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.11/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.11/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.11/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.12/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.12/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.12/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.13/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.13/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.13/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.14/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.14/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.14/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.15/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.15/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.15/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.16/Solution/Three_HH_Solution0w1_0.28w2_0.29w3_0.371.txt');
    V2=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.16/Solution/Three_HH_Solution1w1_0.28w2_0.29w3_0.371.txt');
    V3=load('Causality_Paper_Data/Three_Neuron_Network/Fig10_Case6/Couple=0.16/Solution/Three_HH_Solution2w1_0.28w2_0.29w3_0.371.txt');
    end
  end

  x=V1(1,end-Num+1:end);
  y=V2(1,end-Num+1:end);
  z=V3(1,end-Num+1:end);
  
  x=(x-min(x))./(max(x)-min(x));
  y=(y-min(y))./(max(y)-min(y));
  z=(z-min(z))./(max(z)-min(z));
  max_k=length(x);

if Noise_Add==1
if  SNR_Determined==1
 sigma_white=sqrt(var(x))/10^(x_sigma);
 sigma_white_Noise1=sqrt(var(y))/10^(y_sigma);
 sigma_white_Noise2=sqrt(var(z))/10^(z_sigma);

end
White_Noise=sigma_white.*randn(1,max_k);
White_Noise1=sigma_white_Noise1.*randn(1,max_k);
White_Noise2=sigma_white_Noise2.*randn(1,max_k);

x=White_Noise+x;
y=White_Noise1+y;
z=White_Noise2+z;
end

  X=[x;y;z];
  
 clear V1;
 clear V2;
 clear V3;

sfile = ['ccademo.',num2str(DEMO),'.net'];
nvar = size(X,1);

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
    X=diff(X,1,2);
    uroot_2 = cca_check_cov_stat(X,10);
    inx_2 = find(uroot_2);
    if sum(uroot_2) == 0,
        disp('OK, data is covariance stationary by ADF');
else
    disp('WARNING, data is NOT covariance stationary by ADF');
    disp(['unit roots found in variables: ',num2str(inx_2(1,1))]);  
    end
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
    
len=length(X(1,:));
% Choose a fitting order for GC.
od_max = 100;
[od_joint, od_vec] = chooseOrderFull(X, 'AICc', od_max);
m = max([od_joint, od_vec]);
% For a fast schematic test, use:
% m = chooseOrderAuto(X, 'AICc')

% The Granger Causality value (in matrix form).
GC = nGrangerTfast(X, m);

% Significance test: Non-zero probably based on 0-hypothesis (GC==0).
%p_nonzero = gc_prob_nonzero(GC, m, len);

% % The connectivity matrix.
 p_value = 0.0001;
% net_adjacency = p_nonzero > 1 - p_value;

% This should give the same result.
gc_zero_line = chi2inv(1-p_value, m)/len;
net_adjacency2 = GC > gc_zero_line;
gc_zero_line_Save(1,ie)=gc_zero_line;
% extract the significant causal interactions only
GC2 = GC.*net_adjacency2;
GC2_Values_Save(:,:,ie)=GC2;
GC_Values_Save(:,:,ie)=GC;
end

if Noise_Add==1
 if Cases_choose==1
save('Fig8_Case1_GC_noise.mat','GC_Values_Save')
save('Fig8_Case1_GC2_noise.mat','GC2_Values_Save')
save('Fig8_HH_gc_zero_line_xyy_noise.mat','gc_zero_line_Save')
elseif Cases_choose==2
save('Fig6_Case2_GC_noise.mat','GC_Values_Save')   
save('Fig6_Case2_GC2_noise.mat','GC2_Values_Save')   
save('Fig6_HH_gc_zero_line_xyy_noise.mat','gc_zero_line_Save')
elseif Cases_choose==3
save('Fig7_Case3_GC_noise.mat','GC_Values_Save')   
save('Fig7_Case3_GC2_noise.mat','GC2_Values_Save')   
save('Fig7_HH_gc_zero_line_xyy_noise.mat','gc_zero_line_Save')
elseif Cases_choose==4
save('Fig9_Case4_GC_noise.mat','GC_Values_Save')
save('Fig9_Case4_GC2_noise.mat','GC2_Values_Save')
save('Fig9_HH_gc_zero_line_xyy_noise.mat','gc_zero_line_Save')
elseif Cases_choose==5
save('Fig10_Case6_GC_noise.mat','GC_Values_Save')
save('Fig10_Case6_GC2_noise.mat','GC2_Values_Save')
save('Fig10_HH_gc_zero_line_xyy_noise.mat','gc_zero_line_Save')
end   
else
if Cases_choose==1
save('Fig8_Case1_GC.mat','GC_Values_Save')
save('Fig8_Case1_GC2.mat','GC2_Values_Save')
save('Fig8_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==2
save('Fig6_Case2_GC.mat','GC_Values_Save')   
save('Fig6_Case2_GC2.mat','GC2_Values_Save')   
save('Fig6_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==3
save('Fig7_Case3_GC.mat','GC_Values_Save')   
save('Fig7_Case3_GC2.mat','GC2_Values_Save')   
save('Fig7_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==4
save('Fig9_Case4_GC.mat','GC_Values_Save')
save('Fig9_Case4_GC2.mat','GC2_Values_Save')
save('Fig9_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==5
save('Fig10_Case6_GC.mat','GC_Values_Save')
save('Fig10_Case6_GC2.mat','GC2_Values_Save')
save('Fig10_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
end
end
clear gc_zero_line_Save;
clear GC_Values_Save;
end

