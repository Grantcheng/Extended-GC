% Plot EGCI and NEGCI Figure
EGCI_Data =[0         0    0.1220         0    0.0115
    0.0200    0.2283    0.1842    0.0282    0.0294
    0.0400    0.2297    0.2016    0.0152    0.0128
    0.0600    0.2825    0.2117    0.0518    0.0051
    0.0800    0.2503    0.1903    0.0185    0.0280
    0.1000    0.2263    0.1953    0.0244    0.0268
    0.1100    0.2258    0.2090    0.0220    0.0415
    0.1200    0.2522    0.2281    0.0562    0.0252
    0.1300    0.2609    0.2064    0.0269    0.0183
    0.1400    0.2609    0.2375    0.0237    0.0576
    0.1500    0.2686    0.2168    0.0173    0.0116
    0.1600    0.2835    0.2289    0.0242    0.0319];
    
NEGCI_Data =[0         0              0.0477         0         0.0076
    0.0200    0.0288    0.0457    0.0105    0.0078
    0.0400    0.0737    0.0566    0.0029    0.0113
    0.0600    0.0815    0.0573    0.0156    0.0076
    0.0800    0.0924    0.0501    0.0203    0.0177
    0.1000    0.1082    0.0550    0.0096    0.0148
    0.1100    0.1028    0.0609    0.0242    0.0088
    0.1200    0.1224    0.0569    0.0269    0.0132
    0.1300    0.1306    0.0646    0.0162    0.0055
    0.1400    0.1119    0.0583    0.0382    0.0201
    0.1500    0.1251    0.0561    0.0185    0.0085
    0.1600    0.1447    0.0503    0.0136    0.0120
    0.1700    0.1278    0.0675    0.0270    0.0269
    0.1800    0.1335    0.0629    0.0311    0.0094
    0.1900    0.1425    0.0616    0.0243    0.0113
    0.2000    0.1301    0.0598    0.0214    0.0155];


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
% axis([0 0.16 0 0.35])
% set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.07:0.35);
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

axis([0 0.16 0 0.16])
set(gca,'xtick',[0,0.04,0.08,0.12,0.16],'xticklabel',{'0','0.04','0.08','0.12','0.16'},'fontsize',FONT_SIZE(8));
set(gca,'ytick',[0,0.04,0.08,0.12,0.16],'yticklabel',{'0','0.04','0.08','0.12','0.16'},'fontsize',FONT_SIZE(8));

% Create xlabel
xlabel({'S'},'FontSize',FONT_SIZE(10),'FontName','Arial');
% Create ylabel
ylabel({'Extended GC'},'FontSize',FONT_SIZE(10),'FontName','Arial');

%%
% print(h1,'-depsc2','-r300','HH_3766_All_Couple_Noaddnoise.eps')
print(h2,'-depsc2','-r300','HH_3766_All_Couple_addnoise.eps')