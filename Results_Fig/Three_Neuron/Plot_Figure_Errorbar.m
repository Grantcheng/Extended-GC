%% SNR_X=4;SNR_Y=4;SNR_Z=4.2
% f1=0.280 f2=0.290 f3=0.371
% X->Z; Y->Z
% Results_280290371
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
CGC_xy =[0     0     0     0     0
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
    0.0688    0.0585    0.0603    0.0543    0.0520
    0.1083    0.0844    0.0561    0.0657    0.0864
    0.0851    0.1214    0.0927    0.1065    0.1441
    0.1327    0.1027    0.0936    0.1094    0.0901
    0.1212    0.1261    0.1173    0.1232    0.1337
    0.1126    0.1317    0.1052    0.1352    0.1042
    0.1149    0.1095    0.1509    0.1201    0.1142
    0.1388    0.0727    0.1410    0.1360    0.1043
    0.0878    0.1149    0.1571    0.1046    0.1015
    0.1107    0.1169    0.1347    0.1420    0.1415
    0.1054    0.1411    0.1600    0.1281    0.1644];

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

CGC_yz =[0         0         0         0         0
    0.0523    0.0780    0.0628    0.0584    0.0595
    0.0695    0.0395    0.0756    0.0655    0.0754
    0.0995    0.0888    0.0640    0.0976    0.1051
    0.0867    0.0932    0.1061    0.1133    0.0736
    0.0925    0.1360    0.1161    0.1120    0.1064
    0.1457    0.0888    0.1224    0.1278    0.0792
    0.0897    0.1236    0.1075    0.1335    0.1441
    0.1378    0.1184    0.1469    0.1212    0.1470
    0.0811    0.1465    0.1697    0.1140    0.1223
    0.1416    0.1492    0.1076    0.1438    0.1187
    0.0947    0.1266    0.1245    0.1349    0.1463];

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
errorbar(C,CGC_yz_mean,CGC_yz_std,'MarkerSize',10,'Marker','*','LineWidth',2,...
'Color',[1 0 0]);
errorbar(C,CGC_xy_mean,CGC_xy_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[1 1 0]);
errorbar(C,CGC_yx_mean,CGC_yx_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[0 1 0]);
errorbar(C,CGC_zx_mean,CGC_zx_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[1 0 1]);
errorbar(C,CGC_zy_mean,CGC_zy_std,'MarkerSize',18,'Marker','.','LineWidth',2,...
'Color',[0 0 0]);
% Legend 
hleg1 =legend('x->z','y->z','x->y','y->x','z->x','z->y','NorthWest');
set(hleg1,'Location','NorthWest') 
axis([0 0.16 0 0.2])
set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.02:0.2);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'Conditional Extended GC'},'FontSize',20,'FontName','Arial');


%%
print(h0,'-depsc2','-r300','HH_Three_Neurons_280290371_Addnoise_SNR4_4_4.2_ErrorBar.eps')




