function f = Attractor1(x,y,N,m)
figure(1);
plot(x(2:end),x(1:end-1),'k.','MarkerSize',3); 
title(['Driving attractor Example1 N= ',num2str(N),' m= ',num2str(m)]);%��ͼ�α��� 
xlabel('x(n)');            %��X��˵��
ylabel('x(n-1)');            %��Y��˵��
figure(2);
plot(y(2:end),y(1:end-1),'k.','MarkerSize',3); 
title(['linearly  Driving attractor Example1 N= ',num2str(N),' m= ',num2str(m)]);%��ͼ�α��� 
xlabel('y(n)');            %��X��˵��
ylabel('y(n-1)');            %��Y��˵��