function [fmin] = solvePiecewiseLinear(n,x0,fh)
    % solve 
    cvx_begin quiet
        variable x(n,1)
        minimize fh(x);
        subject to
            norm(x-x0)<=1;
    cvx_end
    fmin = fh(x);
end