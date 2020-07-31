% Bad GC example
% Reverse direction

len = 1e6;

A = zeros(2, 2*9);
[a, de] = gen_hfreq_coef(0.8, 0.1, 8);
A(1, 1:2:2*8) = a;
A(2, 1:2:2*9) = 100*[1 a];
A(2, 2:2:2*8) = a;
De = [de, 0.0; 0.0, 1.0];
X0 = gendata_linear(A, De, len+1e4);
X0 = X0(:, 1e4+1:end);

X0(1,:)=(X0(1,:)-min(X0(1,:)))./(max(X0(1,:))-min(X0(1,:)));
X0(2,:)=(X0(2,:)-min(X0(2,:)))./(max(X0(2,:))-min(X0(2,:)));

x = 3.98.*X0(1,:).*(1-X0(1,:));
y = X0(2,:);

X = [x;y];

GC_show_common;

