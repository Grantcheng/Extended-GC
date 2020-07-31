function causality_paper_GC_Three_HH_demo(DEMO)
clear;
clc;
%ccaStartup;
%-----------------------------------------------------------------------
% FUNCTION: cca_demo.m
% PURPOSE:  demonstrate causal connectivity analysis toolbox
%
% INPUTS:   1   -   dataset from Baccala, L. A. & Sameshima, K. 2001
%           2   -   dataset from Schelter et al. 2006
%
% OUTPUT:   graphical output
%           file 'demo[1,2].net' for Pajek program
%           (see vlado.fmf.uni-lj.si/pub/networks/pajek/)
%
%           Written by Anil K Seth, December 2005
%           Updated April 2009 AKS
% COPYRIGHT NOTICE AT BOTTOM
%-----------------------------------------------------------------------

%   demo parameters
N       =   2000;       % number of observations
PVAL    =   0.01;       % probability threshold for Granger causality significance
NLAGS   =   -1;         % if -1, best model order is assessed automatically

GC_CAL_HOME = fileparts(mfilename('fullpath'));

addpath(GC_CAL_HOME);
addpath([GC_CAL_HOME,'/GC_clean-mainline/GCcal']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/GCcal_spectrum']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/tools_and_utilities']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/prj_neuron_gc']);
addpath([GC_CAL_HOME,'/GC_clean-mainline/prj_neuron_gc/scan_worker_template']);
addpath([GC_CAL_HOME,'/gc_time_frequency']);
addpath([GC_CAL_HOME,'/gc_time_frequency/utilities']);

disp('GC calculation package path added.');
disp(' ');
%%
epstot=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];

%CCases_choose=input('Which case choose: Case_1 = 1, Case_2 = 2, Case_3 = 3, Case_4 = 4, Case_6 = 6 ');

%%
% acquire demo data
if nargin == 0,
    DEMO = 1;
end
Num=1000000;
% X = cca_testData(N,DEMO);
for Cases_choose=1:5
GC_Values_Save = zeros(3,3,12);    
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
  X=[V1(1,end-Num+1:end);V2(1,end-Num+1:end);V3(1,end-Num+1:end)];
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
    
% find best model order
if NLAGS == -1,
    disp('finding best model order ...');
    [bic,aic] = cca_find_model_order(X,2,12);
    disp(['best model order by Bayesian Information Criterion = ',num2str(bic)]);
    disp(['best model order by Aikaike Information Criterion = ',num2str(aic)]);
    NLAGS = aic;
end

%-------------------------------------------------------------------------
% analyze time-domain granger

% find time-domain conditional Granger causalities [THIS IS THE KEY FUNCTION]
disp('finding conditional Granger causalities ...');
ret = cca_granger_regress(X,NLAGS,1);   % STATFLAG = 1 i.e. compute stats

% check that residuals are white
dwthresh = 0.05/nvar;    % critical threshold, Bonferroni corrected
waut = zeros(1,nvar);
for ii=1:nvar,
    if ret.waut<dwthresh,
        waut(ii)=1;
    end
end
inx = find(waut==1);
if isempty(inx),
    disp('All residuals are white by corrected Durbin-Watson test');
else
    disp(['WARNING, autocorrelated residuals in variables: ',num2str(inx)]);
end

% check model consistency, ie. proportion of correlation structure of the
% data accounted for by the MVAR model
if ret.cons>=80,
    disp(['Model consistency is OK (>80%), value=',num2str(ret.cons)]);
else
    disp(['Model consistency is <80%, value=',num2str(ret.cons)]);
end
        
% analyze adjusted r-square to check that model accounts for the data (2nd
% check)
rss = ret.rss_adj;
inx = find(rss<0.3);
if isempty(inx)
    disp(['Adjusted r-square is OK: >0.3 of variance is accounted for by model, val=',num2str(mean(rss))]);
else
    disp(['WARNING, low (<0.3) adjusted r-square values for variables: ',num2str(inx)]);
    disp(['corresponding values are ',num2str(rss(inx))]);
    disp('try a different model order');
end

% find significant Granger causality interactions (Bonferonni correction)
[PR,q] = cca_findsignificance(ret,PVAL,1);
disp(['testing significance at P < ',num2str(PVAL), ', corrected P-val = ',num2str(q)]);

% extract the significant causal interactions only
GC = ret.gc;
GC2 = GC.*PR;
GC_Values_Save(:,:,ie)=GC2;
% calculate causal connectivity statistics
disp('calculating causal connectivity statistics');
causd = cca_causaldensity(GC,PR);


disp(['time-domain causal density = ',num2str(causd.cd)]);
disp(['time-domain causal density (weighted) = ',num2str(causd.cdw)]);

% create Pajek readable file
cca_pajek(PR,GC,sfile);
%-------------------------------------------------------------------------
% plot time-domain granger results
h0=figure(1); clf reset;
FSIZE = 8;
colormap(flipud(bone));

% plot raw time series
for i=2:nvar,
    X(i,:) = X(i,:)+(10*(i-1));
end
subplot(2,3,1);
set(gca,'FontSize',FSIZE);
plot(X');
axis('square');
set(gca,'Box','off');
xlabel('time');
set(gca,'YTick',[]);
xlim([0 N]);
title('Time Series ');

% plot granger causalities as matrix
subplot(2,3,2);
set(gca,'FontSize',FSIZE);
imagesc(GC2);
axis('square');
set(gca,'Box','off');
title(['Granger causality, p<',num2str(PVAL)]);
xlabel('from');
ylabel('to');
set(gca,'XTick',[1:N]);
set(gca,'XTickLabel',1:N);
set(gca,'YTick',[1:N]);
set(gca,'YTickLabel',1:N);

% plot granger causalities as a network
subplot(2,3,3);
cca_plotcausality(GC2,[],5);

if Cases_choose==1
    print(h0,'-depsc2','-r300',['Fig8_Case1_',num2str(ie),'.eps']) 
    elseif Cases_choose==2
    print(h0,'-depsc2','-r300',['Fig6_Case2_',num2str(ie),'.eps'])   
    elseif Cases_choose==3  
    print(h0,'-depsc2','-r300',['Fig7_Case3_',num2str(ie),'.eps'])   
    elseif Cases_choose==4
    print(h0,'-depsc2','-r300',['Fig9_Case4_',num2str(ie),'.eps'])     
    elseif Cases_choose==5
    print(h0,'-depsc2','-r300',['Fig10_Case6_',num2str(ie),'.eps'])
end

end

if Cases_choose==1
save('Fig8_Case1.mat','GC_Values_Save')
elseif Cases_choose==2
save('Fig6_Case2.mat','GC_Values_Save')    
elseif Cases_choose==3
save('Fig7_Case3.mat','GC_Values_Save')        
elseif Cases_choose==4
save('Fig9_Case4.mat','GC_Values_Save')
elseif Cases_choose==5
save('Fig10_Case6.mat','GC_Values_Save')
end

clear GC_Values_Save;
end

