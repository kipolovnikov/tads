function num = division(b, x, y, lenmin)
% division of the sub-network defined by positions of right $x$ and 
% left $y$ boundaries % of a larger network defined by 
% the modularity matrix $b$

    % tolerance
    eps = 10^(-10);

    % size of the sub-network
    n = y - x + 1;
    
    if n >= lenmin 
 

        % modularity matrix for the sub-network
        b1 = b(x:y, x:y);
        for i = 1:n
            b1(i, i) = b1(i, i) - sum(b1(i, :));
        end

        [v, ~] = eig(b1);

        % determining an optimal division $s$ and position of the boundary $num$
        % according to the principal eigenvector
        [num, s] = linear_splitting(v(:, n));

        % only if modularity difference is large enough, we adopt the splitting
        if s'*b1*s < eps || num < lenmin || n-num < lenmin 
           num = -1; 
        end

    else 
        num = -1;
        
    end
end