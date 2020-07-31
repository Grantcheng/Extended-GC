%% Fig5_HH_TWO_Neurons_LinearGC
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
gc_zero_line=1.6132e-04;
significant=1;
if significant==0
%% GC results no significant test
GC_xy=[0.5226e-03 0.0237 0.0275 0.0297 0.0322 0.0228 0.0325 0.0340 0.0359 0.0344 0.0352 -0.0220];
GC_yx=[0.6009e-03 0.0550 0.0552 0.0553 0.0554 0.0555 0.0556 0.0556 0.0555 0.0555 0.0555 0.0554];
else
%% GC2 results have significant test
GC_xy=[0 0.0237 0.0275 0.0297 0.0322 0.0228 0.0325 0.0340 0.0359 0.0344 0.0352 0];
GC_yx=[0 0.0550 0.0552 0.0553 0.0554 0.0555 0.0556 0.0556 0.0555 0.0555 0.0555 0.0554];
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
axis([0 0.16 0 0.1])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','LineraGC_HH_Two_Neurons_Fig5.eps')