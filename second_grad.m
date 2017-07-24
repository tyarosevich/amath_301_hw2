function [DG] = second_grad(x)
    e = exp(1);
    th = x(1);r = x(2);
    d2th = (4*(cos(e^(-2*r/3))^2)*pi*(6*r-r^2)^2)^2/6561;
    d2th_r = -(8*cos(e^(-2*r/3))^2*pi*(-6 + r)*r*(18+(-12+r)*r)*th)/19683;
    d2r_th = d2th_r;
    d2r = (8*cos(e^(-2*r/3))^2*pi*(-3+r)*(-54+(-12+r)*(-9+r)*r)*th^2)/59049;
    DG = [d2th d2th_r;
          d2r_th d2r];
end