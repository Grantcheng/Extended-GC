%% Fig10_Case6_HH_Three_Neurons_KGC
% x->y, y->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.1221 0.0352 0.1156 0.1532 0.0855 0.2913 0.2537 0.1472 0.1016 0.0890 0.0945 0.1134];  % 
GC_xz=[0.0247 0.0114 0.0019 0.0332 0.0019 0.0110 0.0247 0.0314 0.0364 0.0392 0.0406 0.0414];
GC_yx=[0.1036 0.1940 0.0373 0.1456 0.0769 0.0893 0.0942 0.0965 0.1036 0.1133 0.1244 0.1341];  %
GC_yz=[0.0260 0.4750 0.0211 0.1645 0.0168 0.0156 0.0152 0.0172 0.0188 0.0199 0.0210 0.0221];
GC_zx=[0.0440 0.5521 0.0233 0.0173 0.0344 0.0371 0.0276 0.0304 0.0318 0.0330 0.0344 0.0382];
GC_zy=[0.0564 0.2519 0.0236 0.0880 0.0333 0.0408 0.0455 0.0453 0.0441 0.0426 0.0438 0.0467];

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
 
axis([0 0.16 0 1.00])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','KGC_HH_Three_Neurons_Fig10.eps')