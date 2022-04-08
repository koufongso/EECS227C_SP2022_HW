%% 1. Level method
%%Construct A,b
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

%%Level method
KMAX = 100000; % max iteration
x_ = Inf(n,1); % store previous x
x = ones(n,1); % current x , initial x
xx = NaN(n,KMAX); % store all x
ff = NaN(KMAX,1);
gg = NaN(n,KMAX);
ff_opt = NaN(KMAX,1);
xx(:,1) = x;
ff(1) = f(x,A,b);
gg(:,1) = g(x,A,b,0);
ff_opt(1) = f(x,A,b); %current optimum
k = 1;
fopt = min(ff_opt);
while (k<KMAX && norm(x-x_)>1e-5)
    x_ = x;
    %define
    fi = @(t) max(ff(1:k,1)+diag((t*ones(1,k)-xx(:,1:k))'*gg(:,1:k)));

    %step 1: solve min fi(x)
    fi_ = solvePiecewiseLinear(n,x,fi);
    if(isnan(fi_))
        fi_ = x;
    end
    li = 0.5*fi_+0.5*fopt;
    %step 2: project x
    x = levelSetProjection(n,x,fi,li);
 
    k = k+1;
    xx(:,k) = x;
    ff(k) = f(x,A,b);
    gg(:,k) = g(x,A,b,0);
    ff_opt(k) = f(x,A,b);
    fopt = min(ff_opt);
end

disp(fopt);

%%
f_opt_ = -0.481341;
plot(ff_opt-f_opt_);
xlabel('iteration');
ylabel('suboptimality gap');
grid on


%% 2. Convergence rate of GD and AGD for smooth optimization
A = randn(20);
f = @(x) (A*x)'*(A*x);
g = @(x) 2*(A')*A*x;
x0 = ones(20,1);

T = 1e4;
% gradient descent method
[xx,ff_gd] = gradientDescent(x0,f,g,T);
figure
plot(ff_gd);
hold on

% acclerated gradietn descent method
% r = @(x) x'*log(x); %use negative entropy mirror map
% dr = @(x) log(x)+1; % mirror map gradient
% dr_inv = @(dr) exp(dr-1);
% Dr = @(x,y) x'*log(x./y)-sum(x)+sum(y);

r = @(x) x'*x;
dr = @(x) x;
dr_inv = @(dr) dr;
Dr = @(x,y) (x-y)'*(x-y);
[yy,ff_agd] = accelGradientDescent(x0,f,g,Dr,dr,dr_inv,T);
plot(ff_agd);
grid on;
xlabel('iteration');
ylabel('suboptimality gap');
legend('gradient descent','accelerated gradietn descent');

%%
% zoom in
figure
plot(ff_gd);
grid on;
hold on
plot(ff_agd);
xlim([0,100]);
legend('gradient descent','accelerated gradietn descent');

%%
A = diag(2*ones(20,1))+diag(-1*ones(19,1),1)+diag(-1*ones(19,1),-1);
f = @(x) x'*A*x;
g = @(x) 2*A*x;
x0 = ones(20,1);

T = 1e4;
% gradient descent method
[xx,ff_gd] = gradientDescent(x0,f,g,T);
figure
plot(ff_gd);
xlabel('iteration');
ylabel('suboptimality gap');
hold on

% acclerated gradietn descent method
% Regulator 1:
% r = @(x) x'*log(x); %use negative entropy mirror map
% dr = @(x) log(x)+1; % mirror map gradient
% dr_inv = @(dr) exp(dr-1);
% Dr = @(x,y) x'*log(x./y)-sum(x)+sum(y);

% Regulator 2:
r = @(x) x'*x;
dr = @(x) x;
dr_inv = @(dr) dr;
Dr = @(x,y) (x-y)'*(x-y);

[yy,ff_agd] = accelGradientDescent(x0,f,g,Dr,dr,dr_inv,T);
plot(ff_agd);
legend('gradient descent','accelerated gradietn descent');


%%
% zoom in
figure
plot(ff_gd);
grid on;
hold on
plot(ff_agd);
xlim([0,100]);
legend('gradient descent','accelerated gradietn descent');


%% 
% In both cases, the accelerated gradient descent method's performance is
% better than gradietn descent method.

