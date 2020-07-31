function f = Attractor(x,y,yy,N,m)
figure;
plot(x(2:end),x(1:end-1),'k.','MarkerSize',3); 
title(['Driving attractor Example1 N= ',num2str(N),' m= ',num2str(m)]);%��ͼ�α��� 
xlabel('x(n)');            %��X��˵��
ylabel('x(n-1)');            %��Y��˵��

figure;
plot(y(2:end),y(1:end-1),'k.','MarkerSize',3); 
title(['linearly  Driving attractor Example1 N= ',num2str(N),' m= ',num2str(m)]);%��ͼ�α��� 
xlabel('y(n)');            %��X��˵��
ylabel('y(n-1)');            %��Y��˵��

figure;
plot(yy(2:end),yy(1:end-1),'k.','MarkerSize',3); 
title(['nolinearly  Driving attractor Example1 N= ',num2str(N),' m= ',num2str(m)]);%��ͼ�α���
xlabel('yy(n)');            %��X��˵��
ylabel('yy(n-1)');            %��Y��˵��
    end
