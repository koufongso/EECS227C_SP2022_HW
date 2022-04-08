function [result] = fi(i,x,fh,gh)
    result = -inf;
    for j=1:i
        xj = x(:,j);
        y = fh(xj)+(x-xj)'*gh(xj);
        if(y>result)
            result = y;
        end
    end
end