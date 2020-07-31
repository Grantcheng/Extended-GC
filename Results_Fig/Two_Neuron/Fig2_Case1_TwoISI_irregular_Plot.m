% Plot EGCI and NEGCI Figure
EGCI_Data =[
         0         0    0.0056         0    0.0017
    0.0200    0.0009    0.0630    0.0012    0.0133
    0.0400    0.0034    0.0801    0.0007    0.0218
    0.0600    0.0035    0.0773    0.0006    0.0245
    0.0800    0.0047    0.0962    0.0028    0.0209
    0.1000    0.0046    0.0999    0.0008    0.0093
    0.1100    0.0045    0.1258    0.0003    0.0135
    0.1200    0.0062    0.1300    0.0023    0.0309
    0.1300    0.0054    0.1139    0.0019    0.0188
    0.1400    0.0060    0.1187    0.0011    0.0163
    0.1500    0.0063    0.1247    0.0011    0.0159
    0.1600    0.0068    0.1114    0.0016    0.0175];
    
NEGCI_Data =[
         0         0         0         0         0
    0.0200    0.0073         0    0.0009         0
    0.0400    0.0089         0    0.0019         0
    0.0600    0.0083         0    0.0006         0
    0.0800    0.0093         0    0.0007         0
    0.1000    0.0111         0    0.0003         0
    0.1100    0.0139         0    0.0025         0
    0.1200    0.0137         0    0.0035         0
    0.1300    0.0131         0    0.0009         0
    0.1400    0.0123         0    0.0012         0
    0.1500    0.0122         0    0.0030         0
    0.1600    0.0145         0    0.0027         0];


%%%%%%%%%%%%%%%%%%%% parameters of figures used %%%%%%%%%%%%%%%%%%%%%
LEGEND_LOCATION = ['NorthWest';'NorthEast';'SouthWest';'SouthEast']; %% legend position
LINE_STYLE = ['- ';'--';': ';'-.'];
LINE_WIDTH = [0.25:0.25:10.0]; %% the width of drawing lines
MARKER_STYLE = ['.'; '*'; 'o'; '+'; 's'; 'x'; 'd'; '^'; 'v'; '>'; '<'; 'p'; 'h'];
MARKER_SIZE = [1:25];
MARKER_FACE_COLOR = ['r'; 'y'; 'g'; 'c'; 'b'; 'm'; 'k'; 'w'];
COLOR_STYLE = ['r'; 'y'; 'g'; 'c'; 'b'; 'm'; 'k'; 'w'];
FONT_SIZE = 4:4:60;
LEGEND_LOCATION = {'Northeast','Northwest','Southeast','Southwest'};


% %% Plot EGCI
% h1=figure;
% errorbar(EGCI_Data(:,1),EGCI_Data(:,2),EGCI_Data(:,4),'MarkerSize',10,'Marker','*','LineWidth',2,...
% 'Color',[0 0 1]);hold on
% errorbar(EGCI_Data(:,1),EGCI_Data(:,3),EGCI_Data(:,5),'MarkerSize',18,'Marker','.','LineWidth',2,...
% 'Color',[1 0 0]);
% 
% axis([0 0.16 0 0.2])
% set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.05:0.2);
% 
% % Create axis font size
% set(gca,'FontName','Arial','FontSize',14)
% % Create xlabel
% xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% % Create ylabel
% ylabel({'EGCI'},'FontSize',20,'FontName','Arial');

%% Plot NEGCI
h2=figure;
hold on;
errorbar(NEGCI_Data(:,1),NEGCI_Data(:,2),NEGCI_Data(:,4),'markersize',MARKER_SIZE(12),'Marker','*','LineWidth',LINE_WIDTH(6),'Color',[1 0 0]);
errorbar(NEGCI_Data(:,1),NEGCI_Data(:,3),NEGCI_Data(:,5),'markersize',MARKER_SIZE(12),'Marker','.','LineWidth',LINE_WIDTH(6),'Color',[0 0 1]);

axis([0 0.16 0 0.02]);
set(gca,'xtick',[0,0.04,0.08,0.12,0.16],'xticklabel',{'0','0.04','0.08','0.12','0.16'},'fontsize',FONT_SIZE(8));
set(gca,'ytick',[0,0.005,0.01,0.015,0.02],'yticklabel',{'0','0.005','0.01','0.015','0.02'},'fontsize',FONT_SIZE(8));

% Create xlabel
xlabel({'S'},'FontSize',FONT_SIZE(10),'FontName','Arial');
% Create ylabel
ylabel({'Extended GC'},'FontSize',FONT_SIZE(10),'FontName','Arial');

%%
% print(h1,'-depsc2','-r300','HH_28401_All_Couple_Noaddnoise.eps')
print(h2,'-depsc2','-r300','HH_28401_All_Couple_addnoise.eps')