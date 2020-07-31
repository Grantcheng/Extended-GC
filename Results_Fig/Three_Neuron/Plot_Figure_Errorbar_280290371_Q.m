%% SNR_X=4;SNR_Y=4.1;SNR_Z=4.2
% f1=0.280 f2=0.290 f3=0.371
% NeuronType: 1  1  0
% X->Z; Y->Z
% Results_280290371_Q
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

CGC_xz =[0     0     0     0     0
    0.0572    0.0638    0.0564    0.0404    0.0489
    0.0568    0.0496    0.0646    0.0574    0.0712
    0.0482    0.0689    0.0626    0.0458    0.0567
    0.0487    0.0654    0.0637    0.0640    0.0758
    0.0793    0.0646    0.0613    0.0800    0.0915
    0.0791    0.0473    0.0720    0.0740    0.0666
    0.0675    0.0509    0.0722    0.0820    0.0656
    0.0904    0.0496    0.0712    0.0908    0.0790
    0.0814    0.0998    0.0540    0.0935    0.0799
    0.0728    0.0912    0.0801    0.0788    0.0821
    0.0958    0.0906    0.1292    0.0841    0.0878];

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
    0.0607    0.0757    0.0597    0.0420    0.0408
    0.0720    0.0648    0.0738    0.0596    0.0671
    0.0827    0.0469    0.0609    0.0526    0.0412
    0.0675    0.0793    0.0762    0.0546    0.0916
    0.0831    0.0724    0.0739    0.0593    0.0693
    0.0849    0.0744    0.0863    0.0626    0.0624
    0.0699    0.0708    0.0512    0.0861    0.0654
    0.0889    0.0619    0.0715    0.0764    0.0614
    0.0818    0.0735    0.0812    0.0969    0.0764
    0.0532    0.0876    0.0980    0.1017    0.0933
    0.1137    0.1152    0.0795    0.0737    0.0619];

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
print(h0,'-depsc2','-r300','HH_Three_Neurons_280282290_Q_Addnoise_12_12_25_ErrorBar.eps')