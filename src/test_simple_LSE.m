close all;
L = 256;
K = 12;
t = (-K*L:K*L)/L;
Lt = length(t);
C = 3 + 2i;
z = 1 + 1i;
M = exp(-z*t.^2);
x = C*M;
n =  randn(1, Lt) + 1i*randn(1, Lt);
y = x + 0.1*n;

C_LSE = dot(M, y)/dot(M, M);
LSE = C_LSE*M;

R = y(round(L/2))*M;

figure;
plot(abs(y));

figure;
hold on;
plot(abs(x), 'k');
plot(abs(R), 'r--');
plot(abs(LSE), 'g--');