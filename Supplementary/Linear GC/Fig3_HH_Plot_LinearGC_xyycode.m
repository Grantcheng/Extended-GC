%% Fig3_HH_TWO_Neurons_LinearGC
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
gc_zero_line=1.6132e-04;
significant=1;
if significant==0
%% GC results no significant test
GC_xy=[0.1114e-03 0.2634 0.2329 0.1926 0.1560 0.0109 0.0106 0.0119 0.0123 0.0108 0.0125 0.0126];
GC_yx=[0.1015e-03 0.1734 0.1407 0.1098 0.0938 0.0259 0.0235 0.0266 0.0265 0.0251 0.0262 0.0247];
else
%% GC2 results have significant test
GC_xy=[0 0.2634 0.2329 0.1926 0.1560 0.0109 0.0106 0.0119 0.0123 0.0108 0.0125 0.0126];
GC_yx=[0 0.1734 0.1407 0.1098 0.0938 0.0259 0.0235 0.0266 0.0265 0.0251 0.0262 0.0247];
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
 
axis([0 0.16 0 0.3])
set(gca, 'XTick', 0:0.02:0.16); 
% set(gca, 'YTick',0:0.02:0.65);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'GC'},'FontSize',20,'FontName','Arial');

%%
print(h0,'-depsc2','-r300','LineraGC_HH_Two_Neurons_Fig3.eps')