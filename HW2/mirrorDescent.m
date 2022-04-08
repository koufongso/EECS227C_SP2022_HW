function [xbar] = mirrorDescent(x0,fh,dfh,drh,drh_inv,Drh,eta,T)
    n = length(x0);
    xx = zeros(n,T);
    xx(:,1) = x0;
    for i=1:T
        x = xx(:,i);
        w = drh_inv(drh(x)-eta*dfh(x));
        cvx_begin quiet
        variable x_(n,1)
            minimize Drh(x_,w)
        cvx_end
        xx(:,i+1) = x_;
    end
    xbar = mean(xx,2);
end