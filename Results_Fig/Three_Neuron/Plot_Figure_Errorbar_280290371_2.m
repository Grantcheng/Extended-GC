%% SNR_X=4;SNR_Y=4.1;SNR_Z=4.2
% f1=0.280 f2=0.290 f3=0.371
% X->Y; Y->Z
% Results_280290371_2
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
CGC_xy =[0 0 0 0 0 
0.0434    0.0422    0.0372    0.0352    0.0424    
0.0975    0.1078    0.1052    0.0986    0.1040
0.1320    0.1131    0.1176    0.1327    0.1125
0.1150    0.1204    0.1160    0.1132    0.1246
0.1578    0.1857    0.1752    0.1730    0.1897
0.1836    0.1640    0.1869    0.1806    0.2003
0.2021    0.1969    0.2020    0.1906    0.1919
0.1927    0.1898    0.2195    0.2164    0.2112
0.1850    0.1918    0.2092    0.2095    0.1786
0.1975    0.1891    0.2200    0.1924    0.2008
0.1846    0.1716    0.1875    0.2070    0.2131];

CGC_yx =[0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0];

CGC_xz =[0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0];

CGC_zx =[0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0];

CGC_yz =[0 0 0 0 0
0.0731    0.0721    0.0674    0.0684    0.0717    
0.0754    0.0719    0.0904    0.0810    0.0753    
0.1277    0.1176    0.1419    0.1315    0.1314
0.1160    0.1283    0.1219    0.1356    0.1466
0.1731    0.1796    0.1510    0.1490    0.1922
0.1615    0.1428    0.1500    0.1564    0.1334
0.1987    0.1713    0.1955    0.1363    0.1490
0.1632    0.1400    0.1581    0.1854    0.1838
0.1824    0.1563    0.1698    0.1574    0.1568
0.1781    0.2077    0.2041    0.2171    0.2110
0.1352    0.1838    0.1515    0.1494    0.1315];

CGC_zy =[0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0];

% The mean of the CGC 
CGC_xy_mean=mean(CGC_xy,2);
CGC_yx_mean=mean(CGC_yx,2);
CGC_xz_mean=mean(CGC_xz,2);
CGC_zx_mean=mean(CGC_zx,2);
CGC_yz_mean=mean(CGC_yz,2);
CGC_zy_mean=mean(CGC_zy,2);
% The std of the CGC
CGC_xy_std=std(CGC_xy,0,2);
CGC_yx_std=std(CGC_yx,0,2);
CGC_xz_std=std(CGC_xz,0,2);
CGC_zx_std=std(CGC_zx,0,2);
CGC_yz_std=std(CGC_yz,0,2);
CGC_zy_std=std(CGC_zy,0,2);
% Plot
h0=figure;
clf;
errorbar(C,CGC_xy_mean,CGC_xy_std,'MarkerSize',10,'Marker','*','LineWidth',2,...
'Color',[0 0 1]);hold on
errorbar(C,CGC_yz_mean,CGC_yz_std,'MarkerSize',10,'Marker','*','LineWidth',2,...
'Color',[1 0 0]);
errorbar(C,CGC_xz_mean,CGC_xz_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[1 1 0]);
errorbar(C,CGC_yx_mean,CGC_yx_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[0 1 0]);
errorbar(C,CGC_zx_mean,CGC_zx_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[1 0 1]);
errorbar(C,CGC_zy_mean,CGC_zy_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[0 0 0]);
% Legend 
hleg1 =legend('x->y','y->z','x->z','y->x','z->x','z->y','NorthWest');
set(hleg1,'Location','NorthWest') 
axis([0 0.16 0 0.24])
set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.04:0.24);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'Conditional Extended GC'},'FontSize',20,'FontName','Arial');


%%
print(h0,'-depsc2','-r300','HH_Three_Neurons_280282290_2_Addnoise_SNR4_4.1_4.2_ErrorBar.eps')