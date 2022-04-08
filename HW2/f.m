%% f(x)
function [val,active_index] = f(x,A,b)
    k = size(A,3);
    temp = zeros(1,k);
    for index = 1:k
        temp(index) = x'*A(:,:,index)*x-b(:,index)'*x;
    end
    [val,active_index] = max(temp);
end