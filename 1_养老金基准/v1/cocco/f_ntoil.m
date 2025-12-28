function ind = f_ntoil(value,grid,n)
    if (value >= grid(n,1)) 
            ind=n-1;
    elseif (value < grid(1,1)) 
            ind=1;
    else
            ind = 1+fix((value-grid(1,1))/(grid(2,1)-grid(1,1)));
    end
end