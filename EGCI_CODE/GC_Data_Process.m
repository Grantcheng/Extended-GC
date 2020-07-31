Zo=load('mx_10 my_7 Results mat0.11000.17Z_0_Reserve.mat');
Z0=Zo.Z_0_Reserve;
clear Zo
Var=load('mx_10 my_7 Results mat0.11000.17Var_xy.mat');
RSS=load('mx_10 my_7 Results mat0.11000.17Var_yx.mat');
Var=Var.Var_xy;
RSS=RSS.Var_yx;

Count=load('mx_10 my_7 Results mat0.11000.17count_r.mat');
Count=Count.count_r;
G_xy=load('mx_10 my_7 Results mat0.11000.17G_xy_Ftest_Ori.mat');
G_xy=G_xy.G_xy_Ftest_Ori;
G_yx=load('mx_10 my_7 Results mat0.11000.17G_yx_Ftest_Ori.mat');
G_yx=G_yx.G_yx_Ftest_Ori;
clear G_xy_Ft
clc
[A Index]=sort(Count);

G_xy_1=G_xy(:,1);
G_yx_1=G_yx(:,1);

G_xy_Th=load('mx_10 my_7 Results mat0.11000.17GC_xy_Threshold.mat');
G_yx_Th=load('mx_10 my_7 Results mat0.11000.17GC_yx_Threshold.mat');

G_xy_Th=G_xy(:,2);
G_yx_Th=G_yx(:,2);


G_xy_1(find(G_xy_1<0))=100; % the case of GC is small than zero
G_xy_2=G_xy_1(find(G_xy_1>0 & G_xy_1<=1)); % the zero term is that the number of points in this cycle is small than 10*m
G_xy_th2=G_xy_Th(find(G_xy_1>0 & G_xy_1<=1)); % obtain the GC_xy threshold except the GC<0 term
G_xy_Final=G_xy_2(find(G_xy_2>G_xy_th2));% obtain the GC_xy value except the GC<0 term

% G_yx_1(find(G_yx_1<0))=100; % the case of GC_yx is small than zero

G_yx_2=G_yx_1(find(G_xy_1>0 & G_xy_1<=1)); % the zero term is that the number of points in this cycle is small than 10*m
G_yx_th2=G_xy_Th(find(G_xy_1>0 & G_xy_1<=1)); % obtain the GC_yx threshold except the GC_yx<0 term
G_yx_Final=G_yx_2(find(G_yx_2>G_yx_th2));% obtain the GC_yx value except the GC_yx<0 term

%% Plot
G_xy_Small_0=find(G_xy_1==100);% the GC is small than zero , even if the points is large than 10*m
G_xy_No=find(G_xy_1==0);% The number of points is small than 10*m in some cycle

G_yx_Small_0=find(G_yx_1<0);% the GC is small than zero , even if the points is large than 10*m
G_yx_No=find(G_yx_1==0);% The number of points is small than 10*m in some cycle

Whole_Cycle=1:100;

Part_Cycle_Except=union(G_xy_Small_0,G_xy_No);
G_xy_Normal=setdiff(Whole_Cycle,Part_Cycle_Except);% the cycle is normal
G_xy_1(G_xy_1==100)=0;
G_yx_1(find(G_yx_1<0))=0;
Part_Cycle_Except_1=union(G_yx_Small_0,G_yx_No);
G_yx_Normal=setdiff(Whole_Cycle,Part_Cycle_Except_1);% the cycle is normal
G_yx_1(G_yx_1==100)=0;

% ------------- G_xy

figure;[ax,h1,h2]=plotyy(1:100,G_xy_1(Index),1:100,A);
set(get(ax(1),'Ylabel'),'string','GC\_xy','color','k') %y1
set(get(ax(2),'Ylabel'),'string','Number of Points','color','k') %y2

% set(ax(1),'ylim',[-2,2],'ytick',[-2,-1-0,1,1.5,2]);  %����ķ�Χ
% set(ax(2),'ylim',[-5,5],'ytick',[-5:2:5]);  %����ķ�Χ


set(h1,'Marker','*')


xlabel('region Position')
hold on
% plot(1:100,G_xy_Th(Index),'m');
% method = 1 :
index_gc=find(G_xy_Th==0);
[p11 p12 p13]=intersect(index_gc,Index);
 hold on;plot(p13,G_xy_1(index_gc),'m*');


%
hold on;
[p1 p2 p3]=intersect(G_xy_Small_0,Index);
Number_Points_G_xy_Small_0=A(p3);  % The number of points in those GC_xy<0 region;
plot((p3),G_xy_Th(G_xy_Small_0),'r*');
% ------------- G_yx

figure;[ax2,h21,h22]=plotyy(1:100,G_yx_1(Index),1:100,A);
set(get(ax2(1),'Ylabel'),'string','GC\_yx','color','k') %y1
set(get(ax2(2),'Ylabel'),'string','Number of Points','color','k') %y2


% set(ax2(1),'ylim',[-2,2],'ytick',[-2,-1-0,1,1.5,2]);  %����ķ�Χ
% set(ax2(2),'ylim',[-5,5],'ytick',[-5:2:5]);  %����ķ�Χ

set(h21,'Marker','*')


xlabel('region Position')
hold on
plot(1:100,G_yx_Th(Index),'m');
hold on;
[p21 p22 p23]=intersect(G_yx_Small_0,Index);
Number_Points_G_yx_Small_0=A(p23);  % The number of points in those GC_xy<0 region;
plot((p23),G_yx_Th(G_yx_Small_0),'r*');
%% Plot The Attractor

V1=load('HH_28401/HH_Solution_w1_0.28_0.txt');
V2=load('HH_28401/HH_Solution_w1_0.28_1.txt');
Sampling=0.01;
stv=0.01;% the sample interval
dt=0.01;% the data interval
stv_L=stv/dt;
X=[V1;V2];
N=length(X(1,:));
clear V1;
clear V2;
x = X(1,:);
y = X(2,:);
Threshold_Spike_x=(-55-min(x))./(max(x)-min(x));
Threshold_Spike_y=(-55-min(y))./(max(y)-min(y));
x=(x-min(x))./(max(x)-min(x));
y=(y-min(y))./(max(y)-min(y));
X=[x;y];
taux=29;
mx=10;my=7;
tauy=taux % 41
m=mx+my;
[xn,dn1,L1] = PhaSpaRecon(x,taux,mx);
[yn,dn2,L2] = PhaSpaRecon(y,tauy,my);
figure;plot(xn(1,:),xn(2,:),'k.','markersize',1)
hold on;
plot(Z0(1,G_xy_Small_0),Z0(2,G_xy_Small_0),'r.');
hold on;
plot(Z0(1,G_xy_No),Z0(2,G_xy_No),'g*');
hold on;
plot(Z0(1,G_xy_Normal),Z0(2,G_xy_Normal),'b.');
xlabel('x(n)');
ylabel('x(n+1)');
% Z0_index_1=find(Z0(1,:)>0.4 & Z0(1,:)<0.7 & Z0(2,:)>0.8)
% Z0_index_2=find(Z0(1,:)>0.4 & Z0(1,:)<0.7 & Z0(2,:)<0.5 & Z0(2,:)>0.2)
% Z0_index_3=find(Z0(1,:)>0.15 & Z0(1,:)<0.32 & Z0(2,:)<0.5 & Z0(2,:)>0.15)
% Z0_index=find(Z0(1,:)>0.85)
% G_xy_1(Z0_index)
% G_yx(Z0_index,1)
% H_0=Count(Z0_index)
% text(Z0(1,Z0_index(1,2)),Z0(2,Z0_index(1,2)),num2str(18),'color','r');
%
% for t_i=1:length(Z0_index)
% plot(Z0(1,Z0_index(1,t_i)),Z0(2,Z0_index(1,t_i)),'b*','markersize',10);
% text(Z0(1,Z0_index(1,t_i)),Z0(2,Z0_index(1,t_i)),num2str(Z0_index(1,t_i)),'color','r');
% end
%
%%
disp(['Number_Points_G_xy_Small_0 = :'])
 Number_Points_G_xy_Small_0

 disp(['Number_Points_G_yx_Small_0 = :'])
Number_Points_G_yx_Small_0

disp(['G_xy_Final = ',num2str( mean(G_xy_Final))]);
disp(['G_yx_Final = ',num2str( mean(G_yx_Final))]);
disp(' ');
disp(['G_xy_0 = ',num2str( mean(Var))]);
disp(['G_yx_0 = ',num2str( mean(RSS))]);

%%
% X1=x(1,:);
% X1=X1>Threshold_Spike_x;
% aa=find(diff(X1)==1);
% SpikeTiming=aa+1;
% Spike_Vol=x(1,SpikeTiming);
% [Z0_data,Z0_index]=sort(Z0(end,:));
% Z0_cha(1,:)=Z0(1,:);
% Z0_cha(2,:)=Z0(end,:);
% for spike_i=1:length(Z0_index)
% spike_Near(spike_i)=find(Z0_data(1,spike_i)<SpikeTiming,1)-1;
% end
% 
% Spike_Cycle_Inter=Z0_data-SpikeTiming(spike_Near);
% [Spike_Cycle_Inter_data,Spike_Cycle_Inter_index]=sort(Spike_Cycle_Inter);
% 
% Paixu=Z0_index(Spike_Cycle_Inter_index);
% for t_i=1:length(Z0_index)
% hold on;plot(Z0(1,Paixu(1,t_i)),Z0(2,Paixu(1,t_i)),'b*','markersize',10);
% end

%%
y_y=Z0(2,:)<0.9.*Z0(1,:);
y_y1=find(y_y==1);
y_y0=find(y_y==0);
[y_y_y,y_y_index]=sort(Z0(1,y_y1));
% [y_y0_y,y_y0_index]=sort(Z0(2,y_y0),'descend');
[y_y0_y,y_y0_index]=sort(Z0(1,y_y0),'descend');

y_y_all=[y_y1(y_y_index) y_y0(y_y0_index)];
G_xy_all=G_xy(y_y_all,:);
G_yx_all=G_yx(y_y_all,:);
GC_Fina=[G_xy_all Count(y_y_all,1) G_yx_all];

figure;plot(xn(1,:),xn(2,:),'k.','markersize',1)
xlabel('x(n)');
ylabel('x(n+1)');
hold on;plot(Z0(1,y_y1(y_y_index(1:100))),Z0(2,y_y1(y_y_index(1:100))),'b*');
hold on;plot(Z0(1,y_y1(y_y_index(101:200))),Z0(2,y_y1(y_y_index(101:200))),'r*');
hold on;plot(Z0(1,y_y1(y_y_index(201:300))),Z0(2,y_y1(y_y_index(201:300))),'c*');
hold on;plot(Z0(1,y_y1(y_y_index(301:400))),Z0(2,y_y1(y_y_index(301:400))),'m*');
hold on;plot(Z0(1,y_y1(y_y_index(401:end))),Z0(2,y_y1(y_y_index(401:end))),'y*');
hold on;plot(Z0(1,y_y0(y_y0_index(1:100))),Z0(2,y_y0(y_y0_index(1:100))),'g*');
hold on;plot(Z0(1,y_y0(y_y0_index(101:200))),Z0(2,y_y0(y_y0_index(101:200))),'b*');
hold on;plot(Z0(1,y_y0(y_y0_index(201:300))),Z0(2,y_y0(y_y0_index(201:300))),'r*');
hold on;plot(Z0(1,y_y0(y_y0_index(301:end))),Z0(2,y_y0(y_y0_index(301:end))),'g*');

hold on;plot(Z0(1,y_y1(y_y_index(1:14))),Z0(2,y_y1(y_y_index(1:14))),'b*');
hold on;plot(Z0(1,y_y0(y_y0_index(1:10))),Z0(2,y_y0(y_y0_index(1:10))),'g*');
hold on;plot(Z0(1,y_y0(y_y0_index(15:40))),Z0(2,y_y0(y_y0_index(15:40))),'r*');
hold on;plot(Z0(1,y_y0(y_y0_index(41:64))),Z0(2,y_y0(y_y0_index(41:64))),'c*');
hold on;plot(Z0(1,y_y0(y_y0_index(65:86))),Z0(2,y_y0(y_y0_index(65:86))),'m*');

%%

subplot(3,3,1)
plot(1:100,G_yx_all(1:100,1),'b*','markersize',2)
subplot(3,3,2)
plot(101:200,G_yx_all(101:200,1),'b*','markersize',2)
subplot(3,3,3)
plot(201:300,G_yx_all(201:300,1),'b*','markersize',2)
subplot(3,3,4)
plot(301:400,G_yx_all(301:400,1),'b*','markersize',2)
subplot(3,3,5)
plot(401:550,G_yx_all(401:550,1),'b*','markersize',2)
subplot(3,3,6)
plot(551:650,G_yx_all(551:650,1),'b*','markersize',2)
subplot(3,3,7)
plot(651:750,G_yx_all(651:750,1),'b*','markersize',2)
subplot(3,3,8)
plot(751:850,G_yx_all(751:850,1),'b*','markersize',2)
subplot(3,3,9)
plot(851:1000,G_yx_all(851:1000,1),'b*','markersize',2)
%%
figure;[ax,h1,h2]=plotyy(1:100,GC_Fina(:,5),1:100,GC_Fina(:,15)); % Condition number
set(get(ax(1),'Ylabel'),'string','CondNum_{xy}','color','k') %y1
set(get(ax(2),'Ylabel'),'string','CondNum_{yx}','color','k') %y2


figure;[ax,h1,h2]=plotyy(1:100,GC_Fina(:,1),1:100,GC_Fina(:,11));
set(get(ax(1),'Ylabel'),'string','GC_{xy}','color','k') %y1
set(get(ax(2),'Ylabel'),'string','GC_{yx}','color','k') %y2

figure;plot(1:100,GC_Fina(:,5),1:100,GC_Fina(:,15));
figure;plot(1:100,GC_Fina(:,1),1:100,GC_Fina(:,11));



%%
%%
clear;clc;Zo=load('mx_10 my_13System_Results mat0.11000.17Z_0_Reserve.mat');
Z0=Zo.Z_0_Reserve;
clear Zo
Var=load('mx_10 my_13System_Results mat0.11000.17Var_xy.mat');
RSS=load('mx_10 my_13System_Results mat0.11000.17Var_yx.mat');
Var=Var.Var_xy;
RSS=RSS.Var_yx;
Count=load('mx_10 my_13System_Results mat0.11000.17count_r.mat');
Count=Count.count_r;
G_xy=load('mx_10 my_13System_Results mat0.11000.17G_xy_Ftest_Ori.mat');
G_xy=G_xy.G_xy_Ftest_Ori;
G_yx=load('mx_10 my_13System_Results mat0.11000.17G_yx_Ftest_Ori.mat');
G_yx=G_yx.G_yx_Ftest_Ori;
clear G_xy_Ft
clc
[A Index]=sort(Count);
G_xy_1=G_xy(:,1);
G_yx_1=G_yx(:,1);
G_xy_Th=load('mx_10 my_13System_Results mat0.11000.17GC_xy_Threshold.mat');
G_yx_Th=load('mx_10 my_13System_Results mat0.11000.17GC_yx_Threshold.mat');
G_xy_Th=G_xy(:,2);
G_yx_Th=G_yx(:,2);
G_xy_1(find(G_xy_1<0))=100; % the case of GC is small than zero
G_xy_2=G_xy_1(find(G_xy_1>0 & G_xy_1<=1)); % the zero term is that the number of points in this cycle is small than 10*m
G_xy_th2=G_xy_Th(find(G_xy_1>0 & G_xy_1<=1)); % obtain the GC_xy threshold except the GC<0 term
G_xy_Final=G_xy_2(find(G_xy_2>G_xy_th2));% obtain the GC_xy value except the GC<0 term
% G_yx_1(find(G_yx_1<0))=100; % the case of GC_yx is small than zero
G_yx_2=G_yx_1(find(G_xy_1>0 & G_xy_1<=1)); % the zero term is that the number of points in this cycle is small than 10*m
G_yx_th2=G_xy_Th(find(G_xy_1>0 & G_xy_1<=1)); % obtain the GC_yx threshold except the GC_yx<0 term
G_yx_Final=G_yx_2(find(G_yx_2>G_yx_th2));
y_y=Z0(2,:)<0.9.*Z0(1,:);
y_y1=find(y_y==1);
y_y0=find(y_y==0);
[y_y_y,y_y_index]=sort(Z0(1,y_y1));
% [y_y0_y,y_y0_index]=sort(Z0(2,y_y0),'descend');
[y_y0_y,y_y0_index]=sort(Z0(1,y_y0),'descend');
y_y_all=[y_y1(y_y_index) y_y0(y_y0_index)];
G_xy_all=G_xy(y_y_all,:);
G_yx_all=G_yx(y_y_all,:);
GC_Fina=[G_xy_all Count(y_y_all,1) G_yx_all];

figure;[ax,h1,h2]=plotyy(1:100,log10(GC_Fina(:,5)),1:100,log10(GC_Fina(:,6))); % Condition number
set(get(ax(1),'Ylabel'),'string','CondNum_{xy} JR(log10)','color','k') %y1
set(get(ax(2),'Ylabel'),'string','CondNum_{yx} AR(log10)','color','k') %y2

%%
global_prb_gc=[mean(GC_Fina(:,4)) mean(GC_Fina(:,16))];
global_GC_Threshold=mean(GC_Fina(:,3));

Var=mean(GC_Fina(:,1));
RSS=mean(GC_Fina(:,13));
GC_Matrix=zeros(1,2);

GC_Matrix(global_prb_gc<global_GC_Threshold) = 1;

GC_xy=GC_Matrix(1,1)*Var;%
GC_yx=GC_Matrix(1,2)*RSS;%




