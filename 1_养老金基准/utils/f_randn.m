function rand_n = f_randn(n)

    rand_n = zeros(n,1); 
    for ind1=1:n
        rsq = 100;
        while (rsq>=1 || rsq == 0)
            rand_s = rand();
            v1 = 2.0*rand_s-1.0;
            rand_s = rand();
            v2 = 2.0*rand_s-1.0;
            rsq = v1^2+v2^2;
        end
        fac=sqrt(-2.0*log(rsq)/rsq);
        rand_n(ind1,1) =v1*fac;
    end 
    
end
