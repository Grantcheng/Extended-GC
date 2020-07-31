% Plot EGCI and NEGCI Figure
EGCI_Data =[0    0.0244    0.0249    0.0040    0.0038
    0.0200    0.0489    0.0319    0.0052    0.0084
    0.0400    0.0591    0.0180    0.0146    0.0019
    0.0600    0.0601    0.0311    0.0023    0.0104
    0.0800    0.0576    0.0322    0.0253    0.0215
    0.1000    0.0872    0.0580    0.0131    0.0191
    0.1100    0.0911    0.0402    0.0136    0.0117
    0.1200    0.0790    0.0630    0.0160    0.0128
    0.1300    0.0831    0.0404    0.0144    0.0002
    0.1400    0.0712    0.0346    0.0058    0.0059
    0.1500    0.0677    0.0331    0.0054    0.0072
    0.1600    0.1116    0.0623    0.0106    0.0106];
    
NEGCI_Data =[0         0         0         0         0
    0.0200    0.0137         0    0.0028         0
    0.0400    0.0343         0    0.0028         0
    0.0600    0.0313         0    0.0026         0
    0.0800    0.0595         0    0.0065         0
    0.1000    0.0461         0    0.0064         0
    0.1100    0.0410         0    0.0013         0
    0.1200    0.0525         0    0.0044         0
    0.1300    0.0682         0    0.0200         0
    0.1400    0.0653         0    0.0059         0
    0.1500    0.0539         0    0.0046         0
    0.1600    0.0610         0    0.0047         0];


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
% 'Color',[1 0 0]);hold on
% errorbar(EGCI_Data(:,1),EGCI_Data(:,3),EGCI_Data(:,5),'MarkerSize',18,'Marker','.','LineWidth',2,...
% 'Color',[0 0 1]);
% 
% axis([0 0.16 0 0.14])
% set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.04:0.14);
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

axis([0 0.16 0 0.09])
set(gca,'xtick',[0,0.04,0.08,0.12,0.16],'xticklabel',{'0','0.04','0.08','0.12','0.16'},'fontsize',FONT_SIZE(8));
set(gca,'ytick',[0,0.03,0.06,0.09],'yticklabel',{'0','0.03','0.06','0.09'},'fontsize',FONT_SIZE(8));

% Create xlabel
xlabel({'S'},'FontSize',FONT_SIZE(10),'FontName','Arial');
% Create ylabel
ylabel({'Extended GC'},'FontSize',FONT_SIZE(10),'FontName','Arial');

%%
% print(h1,'-depsc2','-r300','HH_371_All_Couple_Noaddnoise.eps')
print(h2,'-depsc2','-r300','HH_371_All_Couple_addnoise.eps')