%% Fig11_HH_TWO_Neurons_LinearGC
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
gc_zero_line=1.6132e-04;
significant=1;
if significant==0
%% GC results no significant test
GC_xy=[0.0314 0.3391 0.2854 0.2157 0.1582 0.5005e-03 0.0013 0.0019 0.2821e-03 0.0017 0.0014 0.6019e-03];
GC_yx=[0.0158 0.0755 0.0601 0.0471 0.0436 0.1370e-03 0.0006 0.0009 0.3385e-03 0.0008 0.0006 0.1905e-03];
else
%% GC2 results have significant test
GC_xy=[0.0314 0.3391 0.2854 0.2157 0.1582 0.5005e-03 0.0013 0.0019 0.2821e-03 0.0017 0.0014 0.6019e-03];
GC_yx=[0.0158 0.0755 0.0601 0.0471 0.0436 0.1370e-03 0.0006 0.0009 0.3385e-03 0.0008 0.0006 0.1905e-03];
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
axis([0 0.16 0 0.5])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','LineraGC_HH_Two_Neurons_Fig11.eps')