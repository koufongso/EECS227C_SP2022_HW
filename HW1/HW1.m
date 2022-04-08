%% 1. Construct A,b
clear all;
close all;
% construct A,b
k = 5;
n = 10;
A = zeros(n,n,k);
b = zeros(n,k);
for index =1:k
    for i = 1:n
        b(i,index) = exp(i/index)*sin(i*index);
        for j = 1:n
            if(i<j)
                A(i,j,index) = exp(i/j)*cos(i*j)*sin(index);
                
            end
        end
    end
end

for index = 1:k
    A(:,:,index) = A(:,:,index) + A(:,:,index)';
end

d = zeros(size(A));
for index = 1:k
    for i=1:n
        A(i,i,index) = (i/10)*abs(sin(index))+sum(abs(A(i,:,index)))-abs(A(i,i,index));
    end
end




%% 1.(a)
x = ones(n,1);
y = f(x,A,b);
fprintf('(a) f(x) = %.3f\n',y);

%% 1.(b) find min f(x)

% param
T = 1000000; % pick a maximum iteration 
C = 1;
eta = C/sqrt(T);

% initial point and tract the optimum solution
x = ones(n,1);
x_opt = x; 
f_min = f(x,A,b);
i_opt = 0;
% find min f(x)
for i = 1: T
    x = x - eta*g(x,A,b,1);
    f_temp = f(x,A,b);
    if(f_temp<f_min)
        f_min = f_temp;
        x_opt = x;
        i_opt = i;
    end
end

fprintf('(b) f(x*) = %.3f\n',f_min);

%% 1.(b) plot gap
x = ones(n,1);
record = zeros(1,T);
for i = 1:T
    x = x - eta*g(x,A,b,1);
    record(i) = f(x,A,b) - f_min;
end

figure(1)
loglog(record,'DisplayName','$\eta = \frac{C}{\sqrt{T}}$');
grid on;
xlabel('iteration');
ylabel('suboptimality gap')




%% 1.(c)
eta =@(i,x) (f(x,A,b)-f_min)/norm(g(x,A,b,0))^2; % Polyak step size
x = ones(n,1);
record = zeros(1,T);

for i = 1:T
    x = x - eta(i,x)*g(x,A,b,1);
    record(i) = f(x,A,b) - f_min;
end

figure(1)
hold on
loglog(record,'DisplayName','Polyak step size');
grid on;
xlabel('iteration');
ylabel('suboptimality gap');
xline(i_opt,'DisplayName','iteration of optimum solution','LineStyle','--');
legend('Interpreter','latex');

fprintf('when away from the optimum solution, the Polyak step size has a higher convergence rate comparing to the fixed step size');







