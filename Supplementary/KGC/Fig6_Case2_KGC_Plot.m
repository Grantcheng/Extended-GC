%% Fig6_Case2_HH_Three_Neurons_KGC
% y->z, x->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0437 0.2975 0.2205 0.2157 0.2197 0.2277 0.2312 0.2343 0.2391 0.2403 0.2410 0.2412];  % 
GC_xz=[0.4269 0.1121 0.0368 0.0224 0.0260 0.0322 0.0348 0.0373 0.0400 0.0423 0.0443 0.0473];
GC_yx=[0.1729 0.1786 0.4018 0.2889 0.2160 0.1711 0.1564 0.1444 0.1397 0.1317 0.1247 0.1177];  %
GC_yz=[0.0200 0.1356 0.1393 0.0299 0.0262 0.0222 0.0202 0.0184 0.0254 0.0243 0.0225 0.0211];
GC_zx=[0.0410 0.1602 0.0381 0.0998 0.1088 0.1103 0.1109 0.1117 0.1164 0.1166 0.1159 0.1165];
GC_zy=[0.0178 0.1361 0.1650 0.1298 0.1180 0.1089 0.1078 0.1087 0.1246 0.1256 0.1258 0.1262];

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
 
axis([0 0.16 0 0.6])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','KGC_HH_Three_Neurons_Fig6.eps')