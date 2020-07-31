%% Fig9_Case4_HH_Three_Neurons_KGC
% x->y, x->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0437 0.0464 0.3058 0.2880 0.1708 0.2682 0.1608 0.2592 0.1809 0.1650 0.1580 0.1777];  % 
GC_xz=[0.4269 0.4176 0.1101 0.1060 0.1485 0.1637 0.1757 0.1709 0.1465 0.1156 0.0626 0.0358];
GC_yx=[0.1729 0.1674 0.1660 0.1553 0.0852 0.1267 0.0769 0.1111 0.0734 0.0761 0.0784 0.1441];  %
GC_yz=[0.0200 0.0140 0.0975 0.1085 0.1105 0.1178 0.1179 0.1847 0.1046 0.1193 0.1037 0.0351];
GC_zx=[0.0410 0.0277 0.1318 0.1624 0.1507 0.1095 0.0815 0.0495 0.1449 0.1375 0.1870 0.1181];
GC_zy=[0.0178 0.0182 0.1082 0.0950 0.0573 0.0797 0.0566 0.0520 0.0567 0.0554 0.0805 0.0513];

% Plot
h0=figure;
clf;
plot(C,GC_xy,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[0 0 1]);hold on
plot(C,GC_yx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 1 0]);
plot(C,GC_xz,'MarkerSize',10,'Marker','s','LineWidth',2,...
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
print(h0,'-depsc2','-r300','KGC_HH_Three_Neurons_Fig9.eps')