%% Fig8_Case1_HH_Three_Neurons_KGC
% x->y, x->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0256 0.0389 0.0383 0.0157 0.0180 0.0195 0.0265 0.0280 0.0282 0.0300 0.0313 0.0334];  % 
GC_xz=[0.2647 0.0264 0.3366 0.0335 0.0319 0.0328 0.0287 0.0248 0.0229 0.0216 0.0214 0.0218];
GC_yx=[0.0292 0.0199 0.0161 0.0087 0.0083 0.0723 0.0166 0.0093 0.0057 0.0071 0.0062 0.0067];  %
GC_yz=[0.0239 0.0916 0.0999 0.4261 0.3227 0.2907 0.2119 0.2423 0.2801 0.3166 0.3417 0.3541];
GC_zx=[0.1675 0.2262 0.2156 0.1417 0.1280 0.0418 0.0566 0.0588 0.0561 0.0534 0.0505 0.0480];
GC_zy=[0.1655 0.0913 0.0241 0.1173 0.0864 0.1589 0.1396 0.1165 0.0971 0.0855 0.0827 0.0877];

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
print(h0,'-depsc2','-r300','KGC_HH_Three_Neurons_Fig8.eps')