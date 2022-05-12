%% Q1. Damped Newton
Mf = 20;
epsilons = [1;0.1;0.01;0.005];
k = length(epsilons);
ff = cell(k,1);
ll = cell(k,1);

for i=1:k
    epsilon = epsilons(i);
    count = 0;
    a = [1:10]'/epsilon;
    
    f = @(x) a'*x-ones(1,10)*log(1-x.^2);
    Df = @(x) a+2*x./(1-x.^2);
    Hf = @(x) diag((1+2*(x.^2))./((1-x.^2).^2));
    lambda = @(x) sqrt(Df(x)'*(Hf(x)\Df(x)));
    
    x = zeros(10,1);
    lambda_= lambda(x);
    ff{i} = [f(x)];
    ll{i} = [lambda_];
    
    while(lambda_>1e-6)
        x = x - (Hf(x)\Df(x))/(1+Mf*lambda_);
        ff{i} = [ff{i},f(x)];
        ll{i} = [ll{i},lambda_];
        lambda_= lambda(x);
        count = count+1;
    end
end

%% plot 
figure
hold on
grid on
xlabel('step')
ylabel('$f(x)$','interpreter','Latex')
for i = 1:k
    plot(ff{i},'displayName',['$\epsilon$=',num2str(epsilons(i))])
end
set(gca,'fontsize',12);
legend('interpreter','latex');

figure
hold on
grid on
xlabel('step')
ylabel('$\lambda$','interpreter','Latex')
for i = 1:k
    plot(ll{i},'displayName',['$\epsilon$=',num2str(epsilons(i))])
end
set(gca,'fontsize',12);
legend('interpreter','latex');
%%

% since the first derivate is a function of $\epsilon$ and the magnitude of
% the first derivate increases as epsilon decreses. As a result, the 
% inital $\lambda$ is bigger when the $\epsilon$ is smaller. The 
% convergence rate is bounded for each epsilon. Therefore, it takes more 
% steps for the funtion with small $\lambda$ to converge.


%% Q2 visualization
delta = 10;
Df = @(x) [2*x(1)/(delta^2-x(1)^2);2*x(2)/(1-x(1)^2)];
Hf = @(x) diag([1/(delta-x(1))^2+1/(delta+x(1))^2,1/(1-x(2))^2+1/(1+x(2))^2]);

beta = delta;
lhs = @(x) norm(-Df(x)+Hf(x)*x);
rhs = @(x) beta*norm(x)^2;
xx=[];
yy=[];
lhs_=[];
rhs_=[];
for x = -0.9:0.05:0.9
    for y = -0.9:0.05:0.9
        xx = [xx,x];
        yy = [yy,y];
        lhs_ = [lhs_,lhs([x;y])];
        rhs_ = [rhs_,rhs([x;y])];
    end
end
figure
hold on
[X,Y] = meshgrid(-0.9:0.05:0.9);

surf(X,Y,reshape(lhs_,size(X)));
surf(X,Y,reshape(rhs_,size(X)));
legend

%%
delta = 1;
f = @(x) abs(4*(x.^3)./(x.^2-delta^2).^2);
x = -1:0.01:1;
plot(x,f(x));
hold on

%%
delta = 100;
[X,Y] = meshgrid(-delta:0.05:delta,-1:0.05:1);
D1 =  X.^2./(X.^2+delta^2) + Y.^2./(Y.^2+1);
contour(X,Y,D1,[1/32,1/32],'r');
D2 = sqrt(X.^2+Y.^2);
hold on
r = sqrt(delta^2+1)-delta;
contour(X,Y,D2,[r,r],'b');



