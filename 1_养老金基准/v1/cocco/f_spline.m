function y2 = f_spline(x,y,n,gam)
    y2 = zeros(n,1);
    u  = zeros(n,1);

    yp1     = x(1,1)^(-gam);
    y2(1,1) = -0.5;
    u(1,1)  = (3.0/(x(2,1)-x(1,1)))*((y(2,1)-y(1,1))/(x(2,1)-x(1,1))-yp1);

    for i=2:n-1
       sig     = (x(i,1)-x(i-1,1))/(x(i+1,1)-x(i-1,1));
       p       = sig*y2(i-1,1)+2.0;
       y2(i,1) = (sig-1.0)/p;
       u(i,1)  = (6.0*((y(i+1,1)-y(i,1))/(x(i+1,1)-x(i,1))-(y(i,1)-y(i-1,1))/(x(i,1)-x(i-1,1)))/(x(i+1,1)-x(i-1,1))-sig*u(i-1,1))/p;
    end
    y2(n,1) = 0.0;
    for k=n-1:-1:1
       y2(k,1) = y2(k,1)*y2(k+1,1)+u(k,1);
    end 
end