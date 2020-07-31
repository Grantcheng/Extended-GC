function [t,y]=rk4(fun,t0,y0,h,N)
%标准四阶Runge_Kutta公式，其中fun为微分方程，t0为初始点，y0为初始向量（列向量）
%h为区间步长，N为区间的个数，x为Xn构成的向量，y为Yn构成的矩阵
t=zeros(1,N+1);y=zeros(length(y0),N+1);
t(1)=t0;
y(:,1)=y0;
for n=1:N
    t(n+1)=t(n)+h;
    k1=feval(fun,t(n),y(:,n));
    k2=feval(fun,t(n)+1/2*h,y(:,n)+h/2*k1);
    k3=feval(fun,t(n)+1/2*h,y(:,n)+h/2*k2);
    k4=feval(fun,t(n)+h,y(:,n)+h*k3);
    y(:,n+1)=y(:,n)+h/6*(k1+2*k2+2*k3+k4);
end