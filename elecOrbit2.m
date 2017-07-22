function Z = elecOrbit2(x0)
    TH = x0(1);
    R = x0(2);
    e = exp(1);
    Z =  - abs(sqrt(2)/81*pi*(6*R - R.^2).*e.^(-R/3).*cos(TH)).^2;
end