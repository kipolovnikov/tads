function g_opt_corr = corrected_g_opt(m, g_opt)

% correction of g_opt:
% obtaining the corrected $g$ as a value at which the local
% slope of the total number of clusters is minimal

    gmin = round(g_opt-1);
    gmax = round(g_opt*2);
    numc = zeros(gmax-gmin+1, 1);

    % for each g splitting the matrix: number of clusters, left & right
    % boundaries
    for g = gmin:gmax
        [numc(g-gmin+1), lb, rb] = spectral_modularity(m, g);
    end

    % $ep$ is a vector of local differencies 
    ep = zeros(length(numc)-1, 1);
    for i = 1:length(numc)-1
        ep(i) = numc(i+1)-numc(i);
    end

    % epsum is a vector of smoothed local differencies
    delt = 3;
    epsum = zeros(length(numc)-delt, 1);
    for i = 1:length(numc)-delt
        for j = 0:delt-1
            epsum(i) = epsum(i) + ep(i+j);
        end
    end

    % g_opt_corr is a value at which the smoothed local difference is
    % minimal
    [val id] = min(epsum);
    g_opt_corr = gmin + id - 1 + round(delt/2);

end