function [DG] = second_grad(x)
    e = exp(1);
    th = x(1);r = x(2);
    d2th = -((4*e^(-2*r/3)*(-6 + r)^2 * r^2*cos(2*th))/(6561*pi));
    d2th_r = (4*e^(-2*r/3)*(-6 + r)*r*(18-12*r+r^2)*sin(2*th))/(19683*pi);
    d2r_th = d2th_r;
    d2r = (8*e^(-2*r/3)*(162-378*r+171*r^2-24*r^3+r^4)*(cos(th))^2)/(59049*pi);
    DG = [d2th d2th_r;
          d2r_th d2r];
end