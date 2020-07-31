%% Fig10_Case6_HH_Three_Neurons_LinearGC
% x->y, y->z
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[0.0334 0.0323 0.0336 0.0331 0.0328 0.0327 0.0338 0.0333 0.0337 0.0331 0.0342 0.0340];  % 
GC_xz=[0.0374 0.0370 0.0369 0.0373 0.0371 0.0375 0.0370 0.0367 0.0364 0.0357 0.0355 0.0355];
GC_yx=[0.0330 0.0316 0.0333 0.0307 0.0303 0.0297 0.0303 0.0295 0.0301 0.0312 0.0325 0.0320];  %
GC_yz=[0.0375 0.0371 0.0366 0.0341 0.0332 0.0311 0.0304 0.0299 0.0299 0.0302 0.0313 0.0295];
GC_zx=[0.0348 0.0342 0.0338 0.0341 0.0337 0.0338 0.0330 0.0332 0.0322 0.0324 0.0320 0.0321];
GC_zy=[0.0355 0.0357 0.0342 0.0331 0.0325 0.0315 0.0304 0.0300 0.0294 0.0303 0.0303 0.0298];

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
 
axis([0 0.16 0 0.04])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','LineraGC_HH_Three_Neurons_Fig10.eps')