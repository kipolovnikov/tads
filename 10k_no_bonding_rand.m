lmin = 2;
lmax = 100;
%f = fopen('10k_no_bonding.txt', 'r');
f = fopen('10k_n1.txt', 'r');

    C = textscan(f, '%f %f');
    n1 = 0;
    n2 = 9999;
    %n2 = 4999;
    n = n2 - n1 + 1;

    m = zeros(n, n);
    i = 1;
    while i <= length(C{1})
        k = C{1}(i) - n1 + 1;
        l = C{2}(i) - n1 + 1;
        if l <= n
            m(k, l) = 1;
            m(l, k) = 1;
        end
        i = i + 1;
    end
    fclose(f);
    
    
    
prob = zeros(n, 1);
for i = 0:n-1
for x = 1:n-i
prob(i+1) = prob(i+1) + m(x, x+i)/(n-i);
end
end



a = zeros(n, n);
s = 0;
for i = 1:n
    %a(i, i) = 0;
    %if i ~=n
    %    a(i, i+1) = 1;
    %end
    %for j = i+2:n
    for j = i:n
        if rand() < prob(j-i+1)
            a(i, j) = 1;
            a(j, i) = a(i, j);
        else
            a(i, j) = 0;
            a(j, i) = a(i, j);
        end
    end
end


lmin = 2;
lmax = 100;
ncont = zeros(n-2, lmax);
for len = lmin:lmax
for i = 1:n-len
ncont(i, len) = sum(sum(a(i:i+len, i:i+len)));
end
end

mcont = zeros(lmax-1, 1);
for len = 2:lmax
    mcont(len) = sum(ncont(:, len))/(n-len);
    %mcont(len) = quantile(ncont(:, len), 0.95);
end



k = 1;
dens = 0;
for i = 1:nn-2
for len = 2:lmax
dens(k) = ncont(i, len)/(len^2+2*len+1);
k = k + 1;
end
end

hist(dens, 1000)



