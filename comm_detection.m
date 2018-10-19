function [num, lbord, rbord] = comm_detection(m)

    % calculating optimal value of g
    g_opt = optimal_g(m);

    % corrected value of g_opt
    g_opt_corr = corrected_g_opt(m, g_opt);

    % obtaining number of clusters & positions of boundaries
    [num, lbord, rbord] = spectral_modularity(m, g_opt_corr);
    
end