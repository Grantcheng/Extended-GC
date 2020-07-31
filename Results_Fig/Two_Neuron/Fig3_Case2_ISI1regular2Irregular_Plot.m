% Plot EGCI and NEGCI Figure
EGCI_Data =[0    0.0143         0    0.0083         0
    0.0200    0.0547         0    0.0144         0
    0.0400    0.0902         0    0.0097         0
    0.0600    0.0894         0    0.0230         0
    0.0800    0.1135         0    0.0128         0
    0.1000    0.1283         0    0.0206         0
    0.1100    0.1137         0    0.0305         0
    0.1200    0.1033    0.0016    0.0078    0.0012
    0.1300    0.1213    0.0018    0.0144    0.0002
    0.1400    0.1278    0.0024    0.0145    0.0003
    0.1500    0.1093    0.0025    0.0175    0.0012
    0.1600    0.1399    0.0038    0.0172    0.0017];
    
NEGCI_Data =[
         0         0         0         0         0
    0.0200    0.0243         0    0.0167         0
    0.0400    0.0551         0    0.0067         0
    0.0600    0.0658         0    0.0155         0
    0.0800    0.0907         0    0.0061         0
    0.1000    0.0854         0    0.0229         0
    0.1100    0.0781         0    0.0120         0
    0.1200    0.0806         0    0.0062         0
    0.1300    0.1008         0    0.0190         0
    0.1400    0.0969         0    0.0100         0
    0.1500    0.0827         0    0.0122         0
    0.1600    0.1059         0    0.0107         0];


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

axis([0 0.16 0 0.13]);
set(gca,'xtick',[0,0.04,0.08,0.12,0.16],'xticklabel',{'0','0.04','0.08','0.12','0.16'},'fontsize',FONT_SIZE(8));
set(gca,'ytick',[0,0.04,0.08,0.12],'yticklabel',{'0','0.04','0.08','0.12'},'fontsize',FONT_SIZE(8));


% Create xlabel
xlabel({'S'},'FontSize',FONT_SIZE(10),'FontName','Arial');
% Create ylabel
ylabel({'Extended GC'},'FontSize',FONT_SIZE(10),'FontName','Arial');

%%
% print(h1,'-depsc2','-r300','HH_280_All_Couple_Noaddnoise.eps')
print(h2,'-depsc2','-r300','HH_280_All_Couple_addnoise.eps')