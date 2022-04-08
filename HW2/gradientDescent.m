function [xx,ff] = gradientDescent(x0,fh,dfh,T)
    xx = nan(size(x0,1),T);
    ff = nan(1,T);
    xx(:,1) = x0;
    ff(1) = fh(x0);
    eta = 0.001;
    for t = 1:T-1
        xx(:,t+1) = xx(:,t)-eta*dfh(xx(:,t))/norm(dfh(xx(:,t)));
        ff(t+1) = fh(xx(:,t+1));
    end
end