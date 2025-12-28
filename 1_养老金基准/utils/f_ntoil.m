function ind = f_ntoil(value,grid,n)
% 寻找value在grid中最接近的位置，超出的话就返回边界的位置
    if (value >= grid(n,1)) 
            ind=n-1;
    elseif (value < grid(1,1)) 
            ind=1;
    else
            ind = 1+fix((value-grid(1,1))/(grid(2,1)-grid(1,1)));
    end
end