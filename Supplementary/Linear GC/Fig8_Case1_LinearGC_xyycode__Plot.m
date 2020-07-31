%% Fig8_Case1_HH_Three_Neurons_LinearGC
% x->y, x->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
gc_zero_line=1.6132e-04;
noise=1;
significant=1;
if significant==0
%% GC results no significant test
if noise==0
GC_xy=[0.5110e-04 0.3545e-05 0.1029e-04 0.0022e-04 0.0021e-03 0.0015e-03 0.0005e-03 0.0205e-04 0.0149e-04 0.0017e-03 0.0018e-03 0.0019e-03];  % 
GC_xz=[0.3154e-04 0.5016e-05 0.0313e-04 0.0006e-04 0.0004e-03 0.0012e-03 0.0008e-03 0.0107e-04 0.0168e-04 0.0017e-03 0.0015e-03 0.0017e-03];
GC_yx=[0.5389e-04 0.1537e-05 0.1092e-04 0.0037e-04 0.0015e-03 0.0010e-03 0.0002e-03 0.0183e-04 0.0109e-04 0.0013e-03 0.0015e-03 0.0015e-03];  %
GC_yz=[0.0416e-04 0.0189e-05 0.0715e-04 0.3137e-04 0.0270e-03 0.0000e-03 0.0017e-03 0.0921e-04 0.3809e-04 0.0172e-03 0.0216e-03 0.0372e-03];
GC_zx=[0.3556e-04 0.7907e-05 0.0648e-04 0.0010e-04 0.0003e-03 0.0009e-03 0.0008e-03 0.0059e-04 0.0129e-04 0.0013e-03 0.0010e-03 0.0012e-03];
GC_zy=[0.0356e-04 0.1114e-05 0.4543e-04 0.0125e-04 0.2061e-03 0.1354e-03 0.1066e-03 0.9450e-04 0.6358e-04 0.1198e-03 0.1256e-03 0.1067e-03];
else
GC_xy=[0.1717e-03 0.1839e-03 0.1262e-03 0.1481e-03 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];  % 
GC_xz=[0.1610e-03 0.1343e-03 0.1375e-03 0.1020e-03 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
GC_yx=[0.2634e-03 0.1766e-03 0.1471e-03 0.1295e-03 0.0001 0.0001 0.0002 0.0001 0.0001 0.0001 0.0001 0.0001];  %
GC_yz=[0.2090e-03 0.1320e-03 0.3290e-03 0.8426e-03 0.0012 0.0018 0.0021 0.0023 0.0026 0.0027 0.0030 0.0033];
GC_zx=[0.2315e-03 0.1998e-03 0.1775e-03 0.1348e-03 0.0002 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
GC_zy=[0.2315e-03 0.2154e-03 0.3793e-03 0.6637e-03 0.0002 0.0018 0.0020 0.0022 0.0024 0.0026 0.0028 0.0030];
end
else
%% GC2 results have significant test
if noise==0
GC_xy=[0.5110e-04 0 0 0 0 0 0 0 0 0 0 0];  % 
GC_xz=[0.3154e-04 0 0 0 0 0 0 0 0 0 0 0];
GC_yx=[0.5389e-04 0 0 0 0 0 0 0 0 0 0 0];  %
GC_yz=[0 0 0 0.3137e-04 0.0270e-03 0 0 0 0.3809e-04 0.0172e-03 0.0216e-03 0.0372e-03];
GC_zx=[0.3556e-04 0 0 0 0 0 0 0 0 0 0 0];
GC_zy=[0 0 0.4543e-04 0 0.2061e-03 0.1354e-03 0.1066e-03 0.9450e-04 0.6358e-04 0.1198e-03 0.1256e-03 0.1067e-03];
else
GC_xy=[0.1717e-03 0.1839e-03 0 0 0 0 0 0 0 0 0 0];  % 
GC_xz=[0 0 0 0 0 0 0 0 0 0 0 0];
GC_yx=[0.2634e-03 0.1766e-03 0 0 0 0 0 0 0 0 0 0];  %
GC_yz=[0.2090e-03 0 0.3290e-03 0.8426e-03 0.0012 0.0018 0.0021 0.0023 0.0026 0.0027 0.0030 0.0033];
GC_zx=[0.2315e-03 0.1998e-03 0.1775e-03 0 0 0 0 0 0 0 0 0];
GC_zy=[0.2315e-03 0.2154e-03 0.3793e-03 0.6637e-03 0.0012 0.0018 0.0020 0.0022 0.0024 0.0026 0.0028 0.0030];    
end
end

% Plot
h0=figure;
clf;
plot(C,GC_xy,'MarkerSize',10,'Marker','s','LineWidth',2,...
'Color',[0 0 1]);hold on
plot(C,GC_yx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 1 0]);
plot(C,GC_xz,'MarkerSize',10,'Marker','s','LineWidth',2,...
'Color',[1 0 1]);
plot(C,GC_zx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 0 0]);
plot(C,GC_yz,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[0 1 0]);
plot(C,GC_zy,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[0 0 0]);
% Legend 
% Legend 
hleg1 =legend('x->y','y->x','x->z','z->x','y->z','z->y','northeast');
set(hleg1,'Location','northeast') 
set(hleg1,'FontName','Times New Roman','FontSize',18,'FontWeight','normal')
 
axis([0 0.16 0 5e-03])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
if noise==0
print(h0,'-depsc2','-r300','LineraGC_HH_Three_Neurons_Fig8.eps')
else
print(h0,'-depsc2','-r300','LineraGC_Noise_HH_Three_Neurons_Fig8.eps')
end