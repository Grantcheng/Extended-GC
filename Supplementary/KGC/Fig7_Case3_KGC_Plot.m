%% Fig7_Case3_HH_Three_Neurons_KGC
% x->y, y->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0437 0.0265 0.1344 0.0624 0.0960 0.0203 0.0102 0.0091 0.0088 0.0101 0.0120 0.0149];  % 
GC_xz=[0.4269 0.1204 0.0785 0.0661 0.1802 0.0899 0.1145 0.1645 0.2192 0.2743 0.3295 0.3854];
GC_yx=[0.1729 0.1541 0.1680 0.0908 0.0203 0.0377 0.0419 0.0337 0.0274 0.0252 0.0283 0.0355];  %
GC_yz=[0.0200 0.0859 0.0566 0.8141 0.0075 0.0453 0.0219 0.0180 0.0205 0.0228 0.0230 0.0212];
GC_zx=[0.0410 0.1278 0.0353 0.0201 0.1884 0.0056 0.0091 0.0092 0.0113 0.0114 0.0130 0.0146];
GC_zy=[0.0178 0.0289 0.0409 0.0371 0.0087 0.0066 0.0145 0.0134 0.0115 0.0099 0.0087 0.0093];

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
 
axis([0 0.16 0 1.0])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','KGC_HH_Three_Neurons_Fig7.eps')