function [numc, lbord, rbord] = spectral_modularity(m, g)       

% spectral modularity algorithm run on a sparse matrix
% with a cutoff for the minimal cluster size (l >= lenmin)

    % minimal length of TAD; set to 2 for 20kb resolution and to 4 for 10kb resolution
    lenmin = 2;
    
    n = length(m);

    diag = 0;
    for i = 1:n
        diag = diag + m(i, i);
    end
    
    % strength of the weighted network
    w = (sum(sum(m)) + diag)/2;

    % vector with nodes degrees
    k = zeros(n, 1);
    for i = 1:n
        k(i) = sum(m(i, :));
    end
    
    % calculating modularity matrix
    b = m;
    for i = 1:n
        for j = i:n
            b(i, j) = m(i, j) - g*k(i)*k(j)/(2*w);
            b(j, i) = b(i, j);
        end
    end

    % hierarchical splitting till there is no more positive gain in modularity
    % score ($num$ equals -1)
    lbord = [1];
    rbord = [n];

    i = 1;
    while i <= length(lbord)
       num = division(b, lbord(i), rbord(i), lenmin);
       while num ~= -1
           if i ~= length(lbord)
               lbord_new = [lbord(1:i), lbord(i) + num, lbord(i+1:length(lbord))];
           else
               lbord_new = [lbord(1:i), lbord(i) + num];
           end

           if i ~= 1
               rbord_new = [rbord(1:i-1), lbord(i) + num - 1, rbord(i:length(rbord))];
           else
               rbord_new = [lbord(i) + num - 1, rbord(i:length(rbord))];
           end
           lbord = lbord_new;
           rbord = rbord_new;

           num = division(b, lbord(i), rbord(i), lenmin);
       end
       i = i + 1;
    end

    % calculating the number of clusters, the size of which is not less than $lenmin and 
    % which have more than two contacts inside
    numc = 0;
    for i = 1:length(lbord)
        len = rbord(i) - lbord(i) + 1;
        %if len >= lenmin && nnz(m(lbord(i):rbord(i), lbord(i):rbord(i))) >= 2
        if len >= lenmin
            numc = numc + 1;
        end
    end

    if numc ~= length(lbord)    
        "TOO SPARSE!!!"
    end
    
    
end
