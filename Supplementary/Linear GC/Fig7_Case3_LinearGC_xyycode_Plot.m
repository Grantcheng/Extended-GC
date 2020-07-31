%% Fig7_Case3_HH_Three_Neurons_LinearGC
% x->y, y->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
gc_zero_line=1.6132e-04;
noise=1;
significant=1;
if significant==0
%% GC results no significant test
if noise==0
GC_xy=[0.3191e-04 0.5076e-05 0.0307e-04 0.0072e-04 0.0045e-04 0.1234e-05 0.0077e-04 0.1029e-05 0.1589e-05 0.1632e-05 0.1388e-05 0.1555e-05];  % 
GC_xz=[0.1228e-04 0.0291e-05 0.0254e-04 0.3406e-04 0.0556e-04 0.3802e-05 0.0860e-04 0.2096e-05 0.0303e-05 0.0001e-05 0.1405e-05 0.7082e-05];
GC_yx=[0.3581e-04 0.7777e-05 0.0520e-04 0.0041e-04 0.0065e-04 0.1608e-05 0.0119e-04 0.1623e-05 0.2509e-05 0.2709e-05 0.2444e-05 0.2773e-05];  %
GC_yz=[0.1337e-04 0.7193e-05 0.1026e-04 0.3225e-04 0.0001e-04 0.0068e-05 0.0050e-04 0.0234e-05 0.0052e-05 0.0048e-05 0.0133e-05 0.0175e-05];
GC_zx=[0.1411e-04 0.0249e-05 0.0119e-04 0.9478e-04 0.1205e-04 0.8463e-05 0.1542e-04 0.4889e-05 0.1481e-05 0.0222e-05 0.1018e-05 0.7309e-05];
GC_zy=[0.1509e-04 0.5715e-05 0.0611e-04 0.4871e-04 0.0000e-04 0.0037e-05 0.0052e-04 0.0306e-05 0.0018e-05 0.0117e-05 0.0258e-05 0.0349e-05];
else
GC_xy=[0.1318e-03 0.1126e-03 0.1354e-03 0.1223e-03 0.1032e-03 0.1177e-03 0.0717e-03 0.1083e-03 0.0893e-03 0.0800e-03 0.0795e-03 0.0947e-03];  % 
GC_xz=[0.1779e-03 0.2132e-03 0.1385e-03 0.1181e-03 0.1532e-03 0.1165e-03 0.1161e-03 0.1118e-03 0.0819e-03 0.1106e-03 0.1259e-03 0.1268e-03];
GC_yx=[0.2158e-03 0.1918e-03 0.1510e-03 0.1571e-03 0.1495e-03 0.1203e-03 0.1108e-03 0.1512e-03 0.1024e-03 0.1729e-03 0.1277e-03 0.1111e-03];  %
GC_yz=[0.1996e-03 0.1644e-03 0.1289e-03 0.1182e-03 0.1018e-03 0.0888e-03 0.0846e-03 0.0870e-03 0.0920e-03 0.0692e-03 0.0984e-03 0.0803e-03];
GC_zx=[0.2602e-03 0.2384e-03 0.2100e-03 0.2045e-03 0.1288e-03 0.1774e-03 0.1504e-03 0.1500e-03 0.1500e-03 0.1514e-03 0.1509e-03 0.1466e-03];
GC_zy=[0.2759e-03 0.2225e-03 0.1869e-03 0.1409e-03 0.1063e-03 0.0750e-03 0.1031e-03 0.1340e-03 0.0898e-03 0.0904e-03 0.1091e-03 0.0839];
end
else
%% GC2 results have significant test
if noise==0
GC_xy=[0.3191e-04 0 0 0 0 0 0 0 0 0 0 0];  % 
GC_xz=[0 0 0 0 0 0 0 0 0 0 0 0];
GC_yx=[0.3581e-04 0 0 0 0 0 0 0 0 0 0 0];  %
GC_yz=[0 0 0 0 0 0 0 0 0 0 0 0];
GC_zx=[0 0 0 0 0 0 0.1542e-04 0 0 0 0 0];
GC_zy=[0 0 0 0 0 0 0 0 0 0 0 0];
else
GC_xy=[0 0 0 0 0 0 0 0 0 0 0 0];  % 
GC_xz=[0.1779e-03 0.2132e-03 0 0 0 0 0 0 0 0 0 0];
GC_yx=[0.2158e-03 0.1918e-03 0 0 0 0 0 0 0 0 0 0];  %
GC_yz=[0.1996e-03 0.1644e-03 0 0 0 0 0 0 0 0.1729e-03 0 0];
GC_zx=[0.2602e-03 0.2384e-03 0.2100e-03 0.2045e-03 0 0.1774e-03 0 0 0 0 0 0];
GC_zy=[0.2759e-03 0.2225e-03 0.1869e-03 0 0 0 0 0 0 0 0 0];    
end
end

% Plot
h0=figure;
clf;
plot(C,GC_xy,'MarkerSize',10,'Marker','s','LineWidth',2,...
'Color',[0 0 1]);hold on
plot(C,GC_yx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 1 0]);
plot(C,GC_xz,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 0 1]);
plot(C,GC_zx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 0 0]);
plot(C,GC_yz,'MarkerSize',10,'Marker','s','LineWidth',2,...
'Color',[0 1 0]);
plot(C,GC_zy,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[0 0 0]);
% Legend 
% Legend 
hleg1 =legend('x->y','y->x','x->z','z->x','y->z','z->y','northeast');
set(hleg1,'Location','northeast') 
set(hleg1,'FontName','Times New Roman','FontSize',18,'FontWeight','normal')
 
axis([0 0.16 0 0.5e-03])
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
print(h0,'-depsc2','-r300','LinearGC_HH_Three_Neurons_Fig7.eps')
else
print(h0,'-depsc2','-r300','LinearGC_Noise_HH_Three_Neurons_Fig7.eps')
end