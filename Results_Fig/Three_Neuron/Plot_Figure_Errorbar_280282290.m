%% SNR_X=4;SNR_Y=4;SNR_Z=4.2
% f1=0.280 f2=0.282 f3=0.290
% X->Y; X->Z
% Results_280282290
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
CGC_xy =[0         0         0         0         0
    0.0345    0.0385    0.0388    0.0421    0.0392
    0.0713    0.0772	0.0732	  0.0714	0.0691
    0.1098    0.0763    0.0999    0.0899    0.0678
    0.0677    0.0533    0.0890    0.0684    0.0656
    0.1076    0.1114    0.1213    0.1818    0.1352
    0.0949    0.0819    0.0636    0.0470    0.0886
    0.0898    0.0750    0.1189    0.1043    0.0939
    0.1241    0.1158    0.1150    0.1200    0.0874
    0.1366    0.0947    0.1068    0.1408    0.1056
    0.0915    0.1197    0.1606    0.1132    0.0971
    0.1289    0.1388    0.1479    0.1573    0.0779];

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

CGC_xz =[0         0         0         0         0
    0.0405    0.0373    0.0498    0.0548    0.0485
    0.1119	  0.1020    0.1001	  0.1117	0.1018
    0.1076    0.0995    0.1041    0.1601    0.0784
    0.0805    0.0904    0.1139    0.0892    0.0845
    0.1644    0.1868    0.1291    0.2246    0.1498
    0.1229    0.1313    0.1351    0.1256    0.1252
    0.1320    0.1388    0.1434    0.1752    0.1269
    0.1636    0.1308    0.1461    0.1853    0.1428
    0.1315    0.1596    0.1569    0.1811    0.1557
    0.1499    0.1506    0.1832    0.1439    0.1510
    0.1761    0.1701    0.1525    0.1871    0.1539];

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

CGC_yz =[0     0     0     0     0
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
errorbar(C,CGC_xz_mean,CGC_xz_std,'MarkerSize',10,'Marker','*','LineWidth',2,...
'Color',[0 0 1]);hold on
errorbar(C,CGC_yz_mean,CGC_yz_std,'MarkerSize',10,'Marker','.','LineWidth',2,...
'Color',[1 1 0]);
errorbar(C,CGC_xy_mean,CGC_xy_std,'MarkerSize',18,'Marker','*','LineWidth',2,...
'Color',[1 0 0]);
errorbar(C,CGC_yx_mean,CGC_yx_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[0 1 0]);
errorbar(C,CGC_zx_mean,CGC_zx_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[1 0 1]);
errorbar(C,CGC_zy_mean,CGC_zy_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[0 0 0]);
% Legend 
hleg1 =legend('x->z','y->z','x->y','y->x','z->x','z->y','NorthWest');
set(hleg1,'Location','NorthWest') 
axis([0 0.16 0 0.22])
set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.02:0.22);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'Conditional Extended GC'},'FontSize',20,'FontName','Arial');


%%
print(h0,'-depsc2','-r300','HH_Three_Neurons_280282290_Addnoise_SNR4_4_4.2_ErrorBar.eps')