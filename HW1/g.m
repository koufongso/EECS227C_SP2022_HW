% function of subgradient
function [subgrad] = g(x,A,b,needNormalize)
    [~,active] = f(x,A,b);
    subgrad = 2*A(:,:,active)*x - b(:,active);
    if (needNormalize) 
        subgrad = subgrad/norm(subgrad,2);
    end
end