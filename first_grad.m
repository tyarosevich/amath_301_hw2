function [G] = first_grad(x)
    e = exp(1);
    th = x(1);r = x(2);
    dth = -((2*e^(-2*r/3)*(-6 + r)^2 * r^2 *sin(2*th))/(6561*pi));
    drad = -((4*e^(-2*r/3)*(-6 + r)*r*(18-12*r+r^2)*(cos(th))^2)/(19683*pi));
    G = [dth; drad];
end