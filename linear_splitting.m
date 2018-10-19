function [num, s] = linear_splitting(v)
% optimal position of the boundary given the principal eigenvector
% in case we are restricted to LINEAR divisions;
% s(i) = 1 or -1;
% "optimal" splitting maximizes (v*s)

    len = length(v);
    s = zeros(len, 1);
    
    % summ(j) = (v*s_j), where $s_j$ has the boundary at position $j$
    summ = zeros(len, 1);
    for j = 1:len
        summ(j) = abs(sum(v(1:j)) - sum(v(j+1:len)));
    end

    % $num$ is the position of the optimal boundary
    [~, num] = max(summ);
    
    if sum(v(1:num)) - sum(v(num+1:len)) > 0
        s(1:num) = 1;
        s(num+1:len) = -1;
    else
        s(1:num) = -1;
        s(num+1:len) = 1;
    end
    
end