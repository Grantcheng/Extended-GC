%% Fig4_HH_TWO_Neurons_LinearGC
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
GC_xy=[6.229e-05 0.6472 0.6377 0.6295 0.6236 5.572e-05 0 0 0 0 0 0];
GC_yx=[5.207e-05 0 0 0 0 8.958e-05 1.184e-04 1.287e-04 4.856e-05 1.563e-04 1.599e-04 2.527e-04];

% Plot
h0=figure;
clf;
plot(C,GC_xy,'MarkerSize',10,'Marker','s','LineWidth',2,...
'Color',[0 0 1]);hold on
plot(C,GC_yx,'MarkerSize',10,'Marker','p','LineWidth',2,...
'Color',[1 0 1]);
% Legend 
hleg1 =legend('x->y','y->x','northeast');
set(hleg1,'Location','northeast') 
set(hleg1,'FontName','Times New Roman','FontSize',18,'FontWeight','normal')
axis([0 0.16 0 0.65])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','LineraGC_HH_Two_Neurons_Fig4.eps')