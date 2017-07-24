function [G] = first_grad(x)
    e = exp(1);
    th = x(1);r = x(2);
    dth = (4*(cos(e^(-2*r/3))^2)*pi*(6*r-r^2)^2*th)/6561;
    drad = (4*(cos(e^(-2*r/3))^2)*pi*(6-2*r)*(6*r-r^2)*th^2)/6561 - (4*(cos(e^(-2*r/3))^2)*pi*(6*r-r^2)^2*th^2)/19683;
    G = [dth; drad];
end