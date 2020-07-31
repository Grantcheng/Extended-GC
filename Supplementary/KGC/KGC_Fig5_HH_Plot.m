%% Fig5_HH_TWO_Neurons_KGC
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_yx=[0 0 0 0 0 0 0 0 0 0 0.1463 0.1476];
GC_xy=[0 0 0 0 0 0 0 0 0 0 0.0984 0.0998];

% Plot
h0=figure;
clf;
plot(C,GC_xy,'MarkerSize',10,'Marker','s','LineWidth',2,...
'Color',[0 0 1]);hold on
plot(C,GC_yx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 0 1]);
% Legend 
% Legend 
hleg1 =legend('Causality x->y','Causality y->x','northeast');
set(hleg1,'Location','northeast') 
set(hleg1,'FontName','Times New Roman','FontSize',18,'FontWeight','normal')
 
axis([0 0.16 0 0.40])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','KGC_HH_Two_Neurons_Fig5.eps')