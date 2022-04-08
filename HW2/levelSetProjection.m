function x = levelSetProjection(n,x0,fh,lv)
    cvx_begin quiet
        variable x(n,1)
        minimize norm(x-x0)
        subject to
            fh(x)<=lv;
            norm(x-x0)<=1;
    cvx_end
end