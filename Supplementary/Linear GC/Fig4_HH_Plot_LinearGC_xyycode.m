%% Fig4_HH_TWO_Neurons_LinearGC
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
gc_zero_line=1.6132e-04;
significant=1;
if significant==0
%% GC results no significant test
GC_xy=[0.1076e-03 0.8366 0.6870 0.5281 0.4048 0.0325 0.0299 0.0288 0.0318 0.0326 0.0372 0.0343];
GC_yx=[0.1108e-03 0.0175 0.0157 0.0124 0.0102 0.0520 0.0503 0.0512 0.0533 0.0529 0.0607 0.0576];
else
%% GC2 results have significant test
GC_xy=[0 0.8366 0.6870 0.5281 0.4048 0.0325 0.0299 0.0288 0.0318 0.0326 0.0372 0.0343];
GC_yx=[0 0.0175 0.0157 0.0124 0.0102 0.0520 0.0503 0.0512 0.0533 0.0529 0.0607 0.0576];
end

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
print(h0,'-depsc2','-r300','LineraGC_HH_Two_Neurons_Fig4.eps')