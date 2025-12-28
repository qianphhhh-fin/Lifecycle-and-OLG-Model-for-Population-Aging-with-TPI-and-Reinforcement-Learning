function y = f_sc_splint(xa,ya,y2a,n,x)
    klo = 1;
    khi = n;
    while (khi-klo>1)
        k=floor((khi+klo)/2);
        if (xa(k,1)>x)
            khi=k;
        else
            klo=k;
        end
    end
    h = xa(khi,1)-xa(klo,1);
    a = (xa(khi,1)-x)/h;
    b = (x-xa(klo,1))/h;
    y = a*ya(klo,1)+b*ya(khi,1)+((a^3-a)*y2a(klo,1)+(b^3-b)*y2a(khi,1))*(h^2)/6.0;
    if (y>ya(n,1)) 
        y = ya(n,1);
    elseif (y<ya(1,1)) 
        y = ya(1,1);
    end
    %if (y<ya(1,1)) 
    %   y = ya(1,1)
    %if
end