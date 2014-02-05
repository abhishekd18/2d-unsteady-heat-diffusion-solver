clear all;
close all;
clc;

a = 2;
b = 1;

lambda = [0:0.1:10];
fun = lambda.*sin(lambda*a)-cos(lambda*a);
l_n = [];
for i=1:numel(lambda)-1
    if(fun(i)*fun(i+1)<0)
        l_n = [l_n, (lambda(i)+lambda(i+1))/2];
    end
end

l_n = unique(l_n);

beta = [0:0.1:10];
fun = beta.*sin(beta*b)-cos(beta*b);
b_m = [];
for i=1:numel(beta)-1
    if(fun(i)*fun(i+1)<0)
        b_m = [b_m, (beta(i)+beta(i+1))/2];
    end
end

b_m = unique(b_m);

T_initial = 1000;
T_inf = 300;

x = [0:0.01/b:a];
y = [0:0.01/a:b];

time = [0.01 0.1 1];

for t=1:length(time)

T_init = ones(numel(x),numel(y))*T_initial;

theta_init = T_init - T_inf;

theta = zeros(numel(x),numel(y));

for i=1:numel(x)
    for j=1:numel(y)
       for n=1:numel(l_n)
           for m=1:numel(b_m)
        theta(i,j) = theta(i,j) + theta_init(i,j)*4*exp(-(l_n(n)^2+b_m(m)^2)*time(t))*...
                     sin(l_n(n)*a)*cos(l_n(n)*x(i))*sin(b_m(m)*b)*cos(b_m(m)*y(j))/...
       (l_n(n)*a + sin(l_n(n)*a)*cos(l_n(n)*a))/(b_m(m)*b + sin(b_m(m)*b)*cos(b_m(m)*b));
   m
           end
       end
    end
end

T = theta + T_inf;

surf(x,y,T');


for i=1:numel(x)
    T_d(i) = T(i,i);
    len(i) = sqrt(x(i)^2 + y(i)^2);
end

f=figure(1);
plot(len,T_d,'r');
hold on;
sol = load(['billet.' num2str(time(t)) '.csv']);
plot(sol(:,3),sol(:,1),'+');
xlabel('diagonal line');
ylabel('Temperature (K)');
title(['Temperature Profile over diagonal at t = ' num2str(time(t)) ' s']);
legend('Analytical Solution','Unstructured FEM Solver');
saveas(f,['billet' num2str(time(t)) '.png'],'png');
hold off;
end


% Calculate coefficients lambda_n and beta_n
% imax = 20;
% imin = 0.5;
% fun = @(lambda) lambda*sin(lambda*a)-cos(lambda*a);
% l_n = [];
% aZero1 = 0;
% tempZero1 = imin;
% while(tempZero1 <= imax)
%     aZero1 = fzero(fun, tempZero1);
%     tempZero1 = aZero1+ 1;
%     l_n = [l_n, aZero1];
% end
% 
% l_n = unique(l_n)
% 
% imax = 10;
% imin = 0.5;
% fun = @(beta) beta*sin(beta*b)-cos(beta*b);
% b_m = [];
% aZero1 = 0;
% tempZero1 = imin;
% while(tempZero1 <= imax)
%     aZero1 = fzero(fun, tempZero1);
%     tempZero1 = aZero1+ 1;
%     b_m = [b_m, aZero1];
% end
% 
% b_m = unique(b_m)