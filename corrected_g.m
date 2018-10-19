
dens = dens_cells();
[i fsor] = sort(dens, 'descend');

%sasha
%2L: 116, 2R: 108, 3L: 140, 3R: 152, X: 89
%g; 2L: 15, 2R: 12, 3L: 22, 3R: 25, X: 12

numcells = 22;

cellnames = ["pop", "sum", "A2", "A3", "A5", "A6", "A8", "A9", "B3", "B6", "B15", "B16", "B19", "B26", "B31", "sc1", "sc16", "sc19", "sc21", "sc23", "sc24", "sc29"];
chrnames = ["chr2L", "chr2R", "chr3L", "chr3R", "chrX"];
chrlen = [2302, 2115, 2455, 2791, 2243];

gopt = zeros(22, 5);
g0_corr = zeros(22, 5);
win = zeros(22, 5);
wout = zeros(22, 5);
eeps = zeros(22, 5);

for p1 = 1:22
    for q = 1:5

        %p = fsor(p1);
        p = p1;
        
        cellnames(p)
        chrnames(q)

        lmin = 2;
        lmax = 100;
        n = chrlen(q);
        m = zeros(n, n);
        filename = "mtx."+cellnames(p)+".10."+chrnames(q)+".txt";
        if p ~= 1
            m = importdata(filename, ' ');
        else
            f = fopen(filename, 'r');
            for i = 1:n
                for j = 1:n
                    m(i, j) = fscanf(f, '%f', 1);
                end
            end
            fclose(f);
        end


        for i = 1:n
            m(i, i) = 0;
        end
        
        %for i = 1:n-1
        %    if m(i, i+1) == 0
        %        m(i, i+1) = 1;
        %        m(i+1, i) = 1;
        %    end
        %end

        %m = m(1:500, 1:500);
        n = length(m);

        w = sum(sum(m))/2;
        k = zeros(n, 1);
        for i = 1:n
            k(i) = sum(m(i, :));
        end

        %digits(100)
        g0 = 1;
        eps = 1;
        it = 0;
        garr1 = 0;
        w_in = 0;
        w_out = 0;
        flag = 0;

        while eps > 0.01 && it < 20
%        while it < 40

        %for i = 1:n
        %    k(i) = sum(m(i, :));
        %    m(i, i) = w_in;
        %end

        w = sum(sum(m))/2;
        for i = 1:n
            k(i) = sum(m(i, :));
        end

        b = m;
        for i = 1:n
            for j = i:n
                b(i, j) = m(i, j) - g0*k(i)*k(j)/(2*w);
                b(j, i) = b(i, j);
            end
        end

        lbord = [1];
        rbord = [n];

        i = 1;
        while i <= length(lbord)
           num = division(b, lbord(i), rbord(i));
           while num ~= -1
        %       lbord
        %       rbord
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

               num = division(b, lbord(i), rbord(i));
           end
           i = i + 1;
        end


        %lbord
        %rbord

        kappa = zeros(length(lbord), 1);
        sumc = zeros(length(lbord), 1);

        m_in = 0;
        for i = 1:length(lbord)
           for j = lbord(i):rbord(i)
               kappa(i) = kappa(i) + sum(m(j, :));
           end
           sumc(i) = (sum(sum(m(lbord(i):rbord(i), lbord(i):rbord(i)))))/2;
           m_in = m_in + sumc(i);
        end

        m_total = w;

        kp = 0;
        for dd = 1:length(kappa)
            kp = kp + kappa(dd)^2;
        end

        w_in = 4*m_total*m_in/kp;
        w_out = (2*m_total-2*m_in)/(2*m_total-kp/(2*m_total));

        g = (w_in - w_out)/(log(w_in) - log(w_out));
        eps = abs(g0-g);
        g0=g;

        it = it + 1;
        garr1(it) = g0;
        
        %for i = 1:length(lbord)
        %    for j = lbord(i):rbord(i)  
        %        m(j, j) = w_in;
        %    end
        %end
        
        end
        
        if eps > 0.01
            g0 = garr1(20)+(garr1(19)+garr1(18)+garr1(17))/4;
        end
        
        %if p == 1 || p == 2
        %    g0 = g0*2;
        %end

        gopt(p, q) = g0;
        win(p, q) = w_in;
        wout(p, q) = w_out;
        eeps(p, q) = eps;

        %g0 = g0*2;
        
        g0
        eps
        win(p, q)
        wout(p, q)    


    gmin = round(g0-1);
    gmax = round(g0*2);
    dltg = 1;
    mlen = zeros((gmax-gmin+1)/dltg + 1);
    numc = zeros((gmax-gmin+1)/dltg + 1);
    g1 = 0;
    for g = gmin:dltg:gmax
        g1 = g1 + 1;
    b = m;
    for i = 1:n
        for j = i:n
            b(i, j) = m(i, j) - g*k(i)*k(j)/(2*w);
            b(j, i) = b(i, j);
        end
    end

    lbord = [1];
    rbord = [n];

    i = 1;
    while i <= length(lbord)
       num = division(b, lbord(i), rbord(i));
       while num ~= -1
    %       lbord
    %       rbord
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

           num = division(b, lbord(i), rbord(i));
       end
       i = i + 1;
    end


    j = 0;
    lens = 0;
    for ii = 1:length(lbord)
    len = rbord(ii) - lbord(ii);
    if len > 3 && nnz(m(lbord(ii):rbord(ii), lbord(ii):rbord(ii))) > 2
        j = j + 1;
        numc(g1) = numc(g1) + 1;
    end
    end

    end
    
    numc = numc(1:nnz(numc));


 
    
    ep = zeros(length(numc)-1, 1);
    for i = 1:length(numc)-1
        ep(i) = numc(i+1)-numc(i);
    end
    
    delt = 3;
    epsum = zeros(length(numc)-delt, 1);
    for i = 1:length(numc)-delt
        for j = 0:delt-1
            epsum(i) = epsum(i) + ep(i+j);
        end
    end
    
    [val id] = min(epsum);
    g0_corrected = gmin + dltg*(id - 1) + round(delt/2)
    gs = gmin:dltg:gmax;
    %figure
    %if q ~= 10
    %    set(gcf,'Visible', 'off'); 
    %end
    %plot(gs, numc)
    %line([g0, g0], [min(numc), max(numc)], 'Color', 'red', 'LineStyle', '--')
    %line([g0_corrected, g0_corrected], [min(numc), max(numc)], 'Color', 'red', 'LineStyle', '-')
    %saveas(gcf, "C:\Users\kipol\Documents\mod_spectra\final\"+cellnames(p)+chrnames(q)+".png");
    %medlen(p, q) = mlen(g_opt);
    %numclst(p, q) = numc(g_opt);
    
    
%

   g0_corr(p, q) = g0_corrected;
   
   w = sum(sum(m))/2;

    %for i = 1:n
    %    k(i) = sum(m(i, :));
    %    m(i, i) = (w_in*k(i)^2)/(2*w-2*k(i)*w_in);
    %end

    %w = sum(sum(m))/2;

    for i = 1:n
        k(i) = sum(m(i, :));
    end

    b = m;
    for i = 1:n
        for j = i:n
            b(i, j) = m(i, j) - g0_corrected*k(i)*k(j)/(2*w);
            b(j, i) = b(i, j);
        end
    end

    lbord = [1];
    rbord = [n];

    i = 1;
    while i <= length(lbord)
       num = division(b, lbord(i), rbord(i));
       while num ~= -1
    %       lbord
    %       rbord
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

           num = division(b, lbord(i), rbord(i));
       end
       i = i + 1;
    end


    kappa = zeros(length(lbord), 1);
    sumc = zeros(length(lbord), 1);

    m_in = 0;
    for i = 1:length(lbord)
       for j = lbord(i):rbord(i)
           kappa(i) = kappa(i) + sum(m(j, :));
       end
       sumc(i) = (sum(sum(m(lbord(i):rbord(i), lbord(i):rbord(i)))))/2;
       m_in = m_in + sumc(i);
    end

    m_total = w;

    kp = 0;
    for dd = 1:length(kappa)
        kp = kp + kappa(dd)^2;
    end

    w_in = 4*m_total*m_in/kp;


    if q == 10
        numcom = length(lbord);
        lambda = zeros(numcom, 1);
        mreal = zeros(numcom, 1);
        for i = 1:numcom
           len = rbord(i) - lbord(i);
           if len > 3 && nnz(m(lbord(i):rbord(i), lbord(i):rbord(i))) > 2
               mreal(i) = sum(sum(m(lbord(i):rbord(i), lbord(i):rbord(i))));
               for d = lbord(i):rbord(i)
                       for s = lbord(i):rbord(i)   
                           lambda(i) = lambda(i) + k(d)*k(s)*w_in/(2*w);
                       end
               end
           end
        end

        figure
        plot(lambda)
        hold on
        plot(mreal)
        pbaspect([1 1 1])
    end


    num = 0;
    for ii = 1:length(lbord)
        len = rbord(ii) - lbord(ii);
        if len > 3 && nnz(m(lbord(ii):rbord(ii), lbord(ii):rbord(ii))) > 2
            num = num + 1;
        end
    end

    
    f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_"+chrnames(q)+"_lb_spmod_sn_r_opt1_delt2", 'w');
    f2 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_"+chrnames(q)+"_rb_spmod_sn_r_opt1_delt2", 'w');
    fprintf(f1, '%d\n', num);
    fprintf(f2, '%d\n', num);
    for ii = 1:length(lbord)
        len = rbord(ii) - lbord(ii);
        if len > 3 && nnz(m(lbord(ii):rbord(ii), lbord(ii):rbord(ii))) > 2
             fprintf(f1, '%d ', lbord(ii));
             fprintf(f2, '%d ', rbord(ii));
        end
    end
    fclose(f1);
    fclose(f2);
    
    end
    
    
%     
%     num = 10000;
%     f = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_ann_spmod_sn_r_opt1_delt2.2Dannot", 'w');
%     fprintf(f, "chr1	x1	x2	chr2	y1	y2	color	comment\n");
% 
%     lb = 0;
%     rb = 0;
%     %chr2L
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr2L"+"_lb_spmod_sn_r_opt1_delt2", 'r');
%     k = fscanf(f1, '%f', 1);
%     for i = 1:k
%         lb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr2L"+"_rb_spmod_sn_r_opt1_delt2", 'r');
%     k1 = fscanf(f1, '%f', 1);
%     for i = 1:k
%         rb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     for i = 1:k
%         fprintf(f, "chr2L	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"chr2L	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"0,0,255	TAD\n");
%     end
% 
%     lb = 0;
%     rb = 0;
%     %chr2R
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr2R"+"_lb_spmod_sn_r_opt1_delt2", 'r');
%     k = fscanf(f1, '%f', 1);
%     for i = 1:k
%         lb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr2R"+"_rb_spmod_sn_r_opt1_delt2", 'r');
%     k1 = fscanf(f1, '%f', 1);
%     for i = 1:k
%         rb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     for i = 1:k
%         fprintf(f, "chr2R	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"chr2R	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"0,0,255	TAD\n");
%     end
% 
%     lb = 0;
%     rb = 0;
%     %chr3L
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr3L"+"_lb_spmod_sn_r_opt1_delt2", 'r');
%     k = fscanf(f1, '%f', 1);
%     for i = 1:k
%         lb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr3L"+"_rb_spmod_sn_r_opt1_delt2", 'r');
%     k1 = fscanf(f1, '%f', 1);
%     for i = 1:k
%         rb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     for i = 1:k
%         fprintf(f, "chr3L	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"chr3L	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"0,0,255	TAD\n");
%     end
% 
%     lb = 0;
%     rb = 0;
%     %chr3R
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr3R"+"_lb_spmod_sn_r_opt1_delt2", 'r');
%     k = fscanf(f1, '%f', 1);
%     for i = 1:k
%         lb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chr3R"+"_rb_spmod_sn_r_opt1_delt2", 'r');
%     k1 = fscanf(f1, '%f', 1);
%     for i = 1:k
%         rb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     for i = 1:k
%         fprintf(f, "chr3R	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"chr3R	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"0,0,255	TAD\n");
%     end
% 
%     lb = 0;
%     rb = 0;
%     %chrX
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chrX"+"_lb_spmod_sn_r_opt1_delt2", 'r');
%     k = fscanf(f1, '%f', 1);
%     for i = 1:k
%         lb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     f1 = fopen("C:\Users\kipol\Documents\MATLAB\"+cellnames(p)+"_chrX"+"_rb_spmod_sn_r_opt1_delt2", 'r');
%     k1 = fscanf(f1, '%f', 1);
%     for i = 1:k
%         rb(i) = fscanf(f1, '%f', 1);
%     end
%     fclose(f1);
% 
%     for i = 1:k
%         fprintf(f, "chrX	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"chrX	"+(lb(i)-1)*num+"	"+(rb(i)-1)*num+"	"+"0,0,255	TAD\n");
%     end
% 
%     fclose(f);
%     
end


f = fopen("C:\Users\kipol\Documents\MATLAB\optimal_gammas", 'w');
for i = 1:22
    for j = 1:5
        fprintf(f, '%f ', gopt(i, j));
    end
    fprintf(f, '\n');
end
fclose(f);  

f = fopen("C:\Users\kipol\Documents\MATLAB\optimal_eps", 'w');
for i = 1:22
    for j = 1:5
        fprintf(f, '%f ', eeps(i, j));
    end
    fprintf(f, '\n');
end
fclose(f);  

f = fopen("C:\Users\kipol\Documents\MATLAB\optimal_win", 'w');
for i = 1:22
    for j = 1:5
        fprintf(f, '%f ', win(i, j));
    end
    fprintf(f, '\n');
end
fclose(f);  

f = fopen("C:\Users\kipol\Documents\MATLAB\optimal_wout", 'w');
for i = 1:22
    for j = 1:5
        fprintf(f, '%f ', wout(i, j));
    end
    fprintf(f, '\n');
end
fclose(f);  

f = fopen("C:\Users\kipol\Documents\MATLAB\optimal_g0_corr", 'w');
for i = 1:22
    for j = 1:5
        fprintf(f, '%f ', g0_corr(i, j));
    end
    fprintf(f, '\n');
end
fclose(f);  