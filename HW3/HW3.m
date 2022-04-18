%% LP Problem setup

%%
% min 1n'*x
% s.t. -1<=xi<=1 for i=1,...,n
% n = 20

n = 20;
% cost and constraints 
c = ones(n,1);
A = [eye(n,n);-1*eye(n,n)];
b = ones(2*n,1);
m = length(b);

% objective function
f = @(t,x) t*c'*x+F(x);             % modified objective function
F = @(x) -sum(log(b-A*x));          % logarithmic barrier function
DF = @(x) A'*((b-A*x).^(-1));       % gradient of F
HF = @(x) A'*diag((b-A*x).^(-2))*A;  % hessian of F
Df = @(t,x) t*c+DF(x);              % gradient of f

% helper function
xnorm_star = @(v,Hx) sqrt(v'*((Hx)\v));
nfunc = @(t,x) -inv(HF(x))*Df(t,x);

%% Central Path Following for LP (Vishnoi ver.)
KMAX = 100000;
% pick initial point, param
x = 0*ones(n,1); 
t = 0.1;
t0 =t;
if(xnorm_star(nfunc(t,x),HF(x))<=1/6)
    disp('inital condition satisfied');
end
% iteration
k=0;
while (k<KMAX && t0*(1+1/(20*sqrt(m)))^k<=(m/err)) 
    x = x+nfunc(t,x);
    t = t*(1+1/(20*sqrt(m)));
    k = k+1;
end
% two Newton steps
x = x+nfunc(t,x);
t = t*(1+1/(20*sqrt(m)));
x = x+nfunc(t,x);
t = t*(1+1/(20*sqrt(m)));

disp(x);


%%
x = 0.3;
diff = 1e-6;
KMAX = 10;
k = 0;
xx= [x];
x_old = x;
while (k<KMAX)
    x_old = x;
    x = -x*log(x);
    xx = [xx,x];
    k = k+1;
end
plot(xx);
yline(exp(-1));
