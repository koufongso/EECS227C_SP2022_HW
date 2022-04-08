%% f(x)
function [val,active_index] = f(x,A,b)
    k = size(A,3);
    temp = zeros(1,k);
    for index = 1:k
        temp(index) = x'*A(:,:,index)*x-b(:,index)'*x;
    end
    [val,active_index] = max(temp);
end

%% subgradient function
function [subgrad] = g(x,A,b,needNormalize)
    [~,active] = f(x,A,b);
    subgrad = 2*A(:,:,active)*x + b(:,active);
    if (needNormalize) 
        subgrad = subgrad/norm(subgrad,2);
    end
end