%% 1) Neuron Type = 1 
% f1=0.280 f2=0.282 f3=0.290
% SNR_X=;SNR_Y=;SNR_Z=
% X->Y; X->Z
%% -----------------------------------------------------------
%% 2) Neuron Type = 2 
% f1=0.280 f2=0.290 f3=0.371
% SNR_X=5.9;SNR_Y=5.9;SNR_Z=5.8
% X->Y; Y->Z
%% -------------------------------------------------------------
neuro_type_choose=input('Which neuron type case choose: Case_1 = 1, Case_2 = 2   ');
C=[0 0.02 0.04 0.06 0.08 0.1 0.11 0.12 0.13 0.14 0.15 0.16];
if neuro_type_choose==1
    
elseif neuro_type_choose==2   
CGC_xy =[0         0         0         0         0
         0.0272    0.0307    0.0298    0.0347    0.0315    % mean=0.0308
         0.0395    0.0237    0.0377    0.0397    0.0273    % mean=0.0336
         0.0542    0.0361    0.0254    0.0442    0.0736   % mean=0.0467
    0.0477    0.0733    0.0468    0.0633    0.0342
    0.0737    0.0534    0.0625    0.0467    0.0503
    0.0623    0.0418    0.0532    0.0638    0.0544
    0.0768    0.0839    0.0620    0.0704    0.0474
    0.0569    0.0555    0.0522    0.0699    0.0772
    0.0790    0.0482    0.0525    0.0805    0.0666
    0.0810    0.0612    0.0666    0.0788    0.0458
    0.0765    0.1229    0.0841    0.0826    0.0470];

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

CGC_yz =[0         0         0         0         0
    0.0439    0.0433    0.0283    0.0753    0.0507 % mean=0.0483
    0.0544    0.0449    0.0810    0.0675    0.0472 % mean=0.0590
    0.0698    0.0559    0.0729    0.0857    0.0710 % mean=0.0711
    0.0776    0.0763    0.1006    0.0919    0.0740
    0.0886    0.0847    0.0977    0.0942    0.0912
    0.1100    0.0968    0.0736    0.0724    0.0841
    0.1064    0.0774    0.1426    0.0965    0.0830
    0.0947    0.1072    0.0848    0.0940    0.1182
    0.0893    0.1066    0.0912    0.1036    0.1038
    0.1132    0.0813    0.1376    0.0959    0.0870
    0.1430    0.0924    0.0793    0.1048    0.0813];

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
end
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
axis([0 0.16 0 0.16])
set(gca, 'XTick', 0:0.02:0.16);  set(gca, 'YTick',0:0.04:0.16);

% Create axis font size
set(gca,'FontName','Arial','FontSize',14)
% Create xlabel
xlabel({'S^C'},'FontSize',20,'FontName','Arial');
% Create ylabel
ylabel({'Conditional Extended GC'},'FontSize',20,'FontName','Arial');


%%
if neuro_type_choose==1
print(h0,'-depsc2','-r300','HH_Three_NeuType_Case1_Addnoise_SNR5.9_5.9_5.8_ErrorBar.eps')
elseif neuro_type_choose==2
print(h0,'-depsc2','-r300','HH_Three_NeuType_Case2_Addnoise_SNR5.9_5.9_5.8_ErrorBar.eps')    
end