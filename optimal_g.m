function g_opt = optimal_g(m)

%determining optimal parameter of splitting for matrix m

    g0 = 1;
    eps = 1;
    it = 0;
    
    n = length(m);

    while eps > 0.01 && it < 20

        [num, lbord, rbord] = spectral_modularity(m, g0);
        kappa = zeros(length(lbord), 1);
        m_in = 0;
        for i = 1:length(lbord)
           diag = 0;
           for j = lbord(i):rbord(i)
               diag = diag + m(j, j);
               kappa(i) = kappa(i) + sum(m(j, :));
           end
           m_in = m_in + (sum(sum(m(lbord(i):rbord(i), lbord(i):rbord(i)))) + diag)/2;
        end

        diag = 0;
        for i = 1:n
            diag = diag + m(i, i);
        end
        m_total = (sum(sum(m)) + diag)/2;

        w_in = 4*m_total*m_in/sum(kappa.^2);
        w_out = (2*m_total-2*m_in)/(2*m_total-sum(kappa.^2)/(2*m_total));

        g = (w_in - w_out)/(log(w_in) - log(w_out));
        eps = abs(g0-g);
        g0=g;
        it = it + 1;
    
    end
    
    g_opt = g0;

end