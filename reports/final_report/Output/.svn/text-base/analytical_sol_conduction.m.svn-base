clear all;
close all;
clc;

a = 2;
alpha = 1;
T_initial = 1000;

x = [0:0.01:a];

time = [0.01 0.1 1];

for t=1:length(time)

T_init = ones(numel(x),1)*T_initial;

T = zeros(numel(x),1);

for i=2:numel(x)-1
   for n=1:10000
        B(n) = -T_init(i)*2*(-1+(-1)^n)/n/pi;
        T(i) = T(i) + B(n)*sin(n*pi*x(i)/a)*exp(-n^2*pi^2*alpha*time(t)/a^2);
   end
end

f=figure(1);
plot(x,T,'r');
hold on;
file = ['1D.' num2str(time(t))];
sol = load([file '.csv']);
plot(sol(:,3),sol(:,1),'+');
xlabel('x at y=0.5');
ylabel('Temperature (K)');
title('Temperature Profile over x at fixed y=0.5');
legend('Analytical Solution','Unstructured FEM Solver');
end

text(0.6,980,['t = ' num2str(time(1)) ' s']);
text(1,900,['t = ' num2str(time(2)) ' s']);
text(1,130,['t = ' num2str(time(3)) ' s']);
saveas(f,'1D.png');