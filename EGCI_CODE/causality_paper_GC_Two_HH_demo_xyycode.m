function causality_paper_GC_Two_HH_demo_xyycode(DEMO)
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

%%
% acquire demo data
if nargin == 0,
    DEMO = 1;
end
Num=1000000;

% X = cca_testData(N,DEMO);
for Cases_choose=2:4
GC_Values_Save = zeros(2,2,12);  
GC2_Values_Save = zeros(2,2,12);  
gc_zero_line_Save=zeros(1,12);
for ie=1:12
%% import data: Three Neurons Network Case
    C=epstot(ie);
    Couple_ie=C;
    display(['The Cases_choose is: ' num2str(Cases_choose)]);
    display(['The Strength is : ' num2str(C)]);
  if Cases_choose==1
   if ie==1
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0/HH_Solution_w1_0.28_1.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.02/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.02/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
   elseif ie==3
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.04/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.04/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.06/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.06/HH_Solution_w1_0.28_w2_0.282_1.txt');
  
    elseif ie==5
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.08/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.08/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.10/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.10/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==7
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.11/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.11/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.12/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.12/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.13/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.13/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.14/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.14/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==11
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.15/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.15/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==12
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.16/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig3_HH/Couple=0.16/HH_Solution_w1_0.28_1.txt');

   end
  elseif Cases_choose==2
      if ie==1
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0/HH_Solution_w1_0.28_1.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.02/HH_Solution_w1_0.28_w2_0_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.02/HH_Solution_w1_0.28_w2_0_1.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.04/HH_Solution_w1_0.28_w2_0_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.04/HH_Solution_w1_0.28_w2_0_1.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.06/HH_Solution_w1_0.28_w2_0_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.06/HH_Solution_w1_0.28_w2_0_1.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.08/HH_Solution_w1_0.28_w2_0_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.08/HH_Solution_w1_0.28_w2_0_1.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.10/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.10/HH_Solution_w1_0.28_1.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.11/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.11/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.12/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.12/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==9
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.13/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.13/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.14/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.14/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.15/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.15/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.16/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig4_HH/Couple=0.16/HH_Solution_w1_0.28_1.txt');

      end
    elseif Cases_choose==3
       if ie==1
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0/HH_Solution_w1_0.371_1.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.02/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.02/HH_Solution_w1_0.371_1.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.04/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.04/HH_Solution_w1_0.371_1.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.06/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.06/HH_Solution_w1_0.371_1.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.08/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.08/HH_Solution_w1_0.371_1.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.10/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.10/HH_Solution_w1_0.371_1.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.11/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.11/HH_Solution_w1_0.371_1.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.12/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.12/HH_Solution_w1_0.371_1.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.13/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.13/HH_Solution_w1_0.371_1.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.14/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.14/HH_Solution_w1_0.371_1.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.15/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.15/HH_Solution_w1_0.371_1.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.16/HH_Solution_w1_0.371_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig5_HH/Couple=0.16/HH_Solution_w1_0.371_1.txt');

       end       
    elseif Cases_choose==4
    if ie==1
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0/HH_Solution_w1_0.28_1.txt');
    
   elseif ie==2
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.02/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.02/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
    elseif ie==3
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.04/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.04/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
    elseif ie==4
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.06/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.06/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
    elseif ie==5
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.08/HH_Solution_w1_0.28_w2_0.282_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.08/HH_Solution_w1_0.28_w2_0.282_1.txt');
    
    elseif ie==6
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.10/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.10/HH_Solution_w1_0.28_1.txt');
         
    elseif ie==7
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.11/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.11/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==8
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.12/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.12/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==9
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.13/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.13/HH_Solution_w1_0.28_1.txt');
    
    elseif ie==10
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.14/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.14/HH_Solution_w1_0.28_1.txt');
        
    elseif ie==11
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.15/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.15/HH_Solution_w1_0.28_1.txt');
     
    elseif ie==12
    V1=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.16/HH_Solution_w1_0.28_0.txt');
    V2=load('Causality_Paper_Data/Two_Neuron_Network/Fig11_HH/Couple=0.16/HH_Solution_w1_0.28_1.txt');
    end      
  end
  
  x=V1(1,end-Num+1:end);
  y=V2(1,end-Num+1:end);
    
%   x=(x-min(x))./(max(x)-min(x));
%   y=(y-min(y))./(max(y)-min(y));
  
  X=[x;y];
  
 clear V1;
 clear V2;

sfile = ['ccademo.',num2str(DEMO),'.net'];
nvar = size(X,1);
Stationary=0;
% detrend and demean data
disp('detrending and demeaning data');
X = cca_detrend(X);
X = cca_rm_temporalmean(X);
X=X+20;
% check covariance stationarity
disp('checking for covariance stationarity ...');
uroot = cca_check_cov_stat(X,10);
inx = find(uroot);
if sum(uroot) == 0,
    disp('OK, data is covariance stationary by ADF');
else
    disp('WARNING, data is NOT covariance stationary by ADF');
    disp(['unit roots found in variables: ',num2str(inx(1,1))]); 
    Y0=log(X(1,:));
    Y1=log(X(2,:));
    count=0;
 while Stationary<1   
   % X=diff(log(X),1,2);
    Y=diff(Y0);
    Y2=diff(Y1);
    X=[Y;Y2];
    uroot_2 = cca_check_cov_stat(X,10);
    inx_2 = find(uroot_2);
    if sum(uroot_2) == 0,
        disp('OK, data is covariance stationary by ADF');
        Stationary=1;
    else
    Y0=Y;
    Y1=Y2;
    count=count+1;
    disp('WARNING, data is NOT covariance stationary by ADF');
    disp(['unit roots found in variables: ',num2str(inx_2(1,1))]);  
    if count>10
        break;        
    end
    end
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

if Cases_choose==1
save('Fig3_HH_GC_xyy.mat','GC_Values_Save')
save('Fig3_HH_GC2_xyy.mat','GC2_Values_Save')
save('Fig3_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==2
save('Fig4_HH_GC_xyy.mat','GC_Values_Save')  
save('Fig4_HH_GC2_xyy.mat','GC2_Values_Save')  
save('Fig4_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==3
save('Fig5_HH_GC_xyy.mat','GC_Values_Save')   
save('Fig5_HH_GC2_xyy.mat','GC2_Values_Save')  
save('Fig5_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
elseif Cases_choose==4
save('Fig11_HH_GC_xyy.mat','GC_Values_Save')  
save('Fig11_HH_GC2_xyy.mat','GC2_Values_Save')  
save('Fig11_HH_gc_zero_line_xyy.mat','gc_zero_line_Save')
end
clear gc_zero_line_Save;
clear GC_Values_Save;
clear GC2_Values_Save;
end

