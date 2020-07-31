%% Fig9_Case4_HH_Three_Neurons_LinearGC
% x->y, x->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0077 0.0081 0.0085 0.0086 0.0087 0.0088 0.0088 0.0088 0.0089 0.0089 0.0089 0.0089];  % 
GC_xz=[0.0166 0.0068 0.0036 0.0015 0.0009 0.0007 0.0003 0.0003 0.0003 0.0003 0.0002 0.0002];
GC_yx=[0.0063 0.0069 0.0072 0.0074 0.0075 0.0075 0.0076 0.0076 0.0076 0.0076 0.0077 0.0077];  %
GC_yz=[0.0129 0.0053 0.0032 0.0013 0.0007 0.0003 0.0004 0.0003 0.0002 0.0002 0.0003 0.0003];
GC_zx=[0.0111 0.0060 0.0040 0.0025 0.0019 0.0016 0.0010 0.0011 0.0008 0.0008 0.0008 0.0005];
GC_zy=[0.0092 0.0048 0.0038 0.0021 0.0017 0.0010 0.0012 0.0010 0.0007 0.0006 0.0009 0.0008];

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
print(h0,'-depsc2','-r300','LineraGC_HH_Three_Neurons_Fig9.eps')