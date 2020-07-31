%% Fig7_Case3_HH_Three_Neurons_LinearGC
% x->y, y->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0077 0.0065 0.0047 0.0031 0.0022 0.0017 0.0014 0.0012 0.0010 0.0009 0.0008 0.0007];  % 
GC_xz=[0.0166 0.0124 0.0115 0.0059 0.0043 0.0032 0.0032 0.0028 0.0027 0.0025 0.0022 0.0021];
GC_yx=[0.0063 0.0060 0.0045 0.0035 0.0030 0.0025 0.0023 0.0020 0.0019 0.0017 0.0016 0.0015];  %
GC_yz=[0.0129 0.0086 0.0052 0.0020 0.0010 0.0005 0.0004 0.0003 0.0002 0.0002 0.0001 0.0001];
GC_zx=[0.0111 0.0085 0.0085 0.0053 0.0043 0.0034 0.0034 0.0032 0.0031 0.0028 0.0025 0.0026];
GC_zy=[0.0092 0.0063 0.0039 0.0016 0.0007 0.0004 0.0003 0.0002 0.0002 0.0001 0.0001 0.0001];

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
 
axis([0 0.16 0 0.02])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','LineraGC_HH_Three_Neurons_Fig7.eps')