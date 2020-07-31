function [count_r,neighbors]=CGC_Three_Neighbors(pointset,detal,N)

global Threshold_Spike_x;
global Ref_Number
global Reference_Points_Choose; % The Reference points are chose at spike region
global Ration_S_NS;

ref=Ref_Number;
disp(['The Ref_Number is: ',num2str(ref)])

S_NS=Ration_S_NS*ref
%% -------The number of reference points around the attractor is ref
if Reference_Points_Choose==0
K=randperm(length(pointset'));
pointset_1=pointset(:,1:N-1)';
atria = nn_prepare(pointset_1, 'euclidian');
referenceset=K(1:ref);
elseif Reference_Points_Choose==1
HH=pointset(1,:)>Threshold_Spike_x;

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