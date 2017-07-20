%AMATH 301 Homework 2

%% Exercise 1 SOR

clc; clear all; close all;

A = diag(zeros(1,50)+2) + diag((zeros(1,49)-1), -1) + diag((zeros(1,49) -1), 1);
p = zeros(50,1);
for n = 1:50
    p(n,1) = 2*(1-cos(23*pi/51))*sin(23*pi*n/51);
end

D = diag(diag(A));
L = tril(A) - D;
U = triu(A) - D;
w = 1.5;
S = inv(D+w*L);

c = S*w*p;
M = -S*(w*U+(w-1)*D);
A1 = zeros(50,51);
A1(1:50, 1:50) = M;
A1(1:50, 51) = c;
%% Exercise 1.b.)

A2Holder = zeros(100, 52);
count = 1;
for w_k = 1:0.01:1.99
    S_k = inv(D+w_k*L);
    M_k = -S_k*(w_k*U+(w_k-1)*D);
    A2Holder(count,1:50) = abs(eig(M_k));
    A2Holder(count,51) = w_k;
    count = count + 1;
end

for n = 1:100
    A2Holder(n, 52) = max(A2Holder(n, 1:50));
end
[minEigValue, minEigIndex] = min(A2Holder(1:100,52));

A2 = A2Holder(1:100, 1:50);
A3 = [A2Holder(minEigIndex, 51);A2Holder(minEigIndex, 52)];
save A2.dat A2 -ASCII
save A3.dat A3 -ASCII
%% Exercise 1.c.)


iterations = 0;
A4 = zeros(100,1);

for w_j = 1:0.01:1.99
    S_j = inv(D+w_j*L);
    M_j = -S_j*(w_j*U+(w_j-1)*D);
    c_j = S_j*w_j*p;
    iterations = iterations + 1;
    x_0 = 1 + zeros(50,1);
    x_k = x_0;

    for n = 1:200
        x_0 = x_k;
        x_k= M_j*x_0 + c_j;
    end
    A4(iterations,1) = norm((A*x_k - p), Inf);
    
end
%% Exercise 1.d.)
x0 = zeros(50,1) + 1;
x = x0;
w_opt = 1.89;
iterations = 0;
tol = 1.e-4;
error = 2*tol;
S_optimal = inv(D+w_opt*L);
c_optimal = S_optimal*w_opt*p;
M_optimal = -S_optimal*(w_opt*U+(w_opt-1)*D);
while ((error > tol) && (iterations < 100000))
    iterations = iterations + 1;
    x0 = x;
    x = M_optimal*x0 + c_optimal;
    error = norm((x-x0), Inf);
end
%% Exercise 2
clc;close all;clear all;
load salmon_data.csv
t = (1:length(salmon_data)).';
Q_11 = sum(t.^2); Q_12 = sum(t); Q_21 = Q_12; Q_22 = 77;
R_11 = sum(t.*salmon_data); R_21 = sum(salmon_data);
Q = [Q_11 Q_12;Q_21 Q_22];
R = [R_11; R_21];
A = Q\R;
y_coeffs = A(1,1)*t + A(2,1);
%plot(t, salmon_data)
%hold on
%plot(t, y_coeffs, 'g')
%p = polyfit(t, salmon_data, 1);
%y_coeffs2 =  p(1,1)*t + p(1,2)
%plot(t, y_coeffs2, 'b--o')

%% Exercise 2.b.)

poly2 = polyfit(t, salmon_data, 2);
poly5 = polyfit(t, salmon_data, 5);
poly8 = polyfit(t, salmon_data, 8);
save A10.dat poly2 -ASCII
save A11.dat poly5 -ASCII
save A12.dat poly8 -ASCII



ycoeffs2 = polyval(poly2, t);
ycoeffs5 = polyval(poly5, t);
ycoeffs8 = polyval(poly8, t);

%plot(t, ycoeffs2, 'g')
% hold on
% plot(t, ycoeffs5, 'b')
% plot(t, ycoeffs8, 'k')
% plot(t, salmon_data)

%% Exercise 2. c.)

%Wrote my own polyval function, then remembered how polyval worked. This
%is probably what the guts of polyval look like, albeit in C though!
% test55 = 0
% for n = 1:length(poly5)
%      test55 = test55 + poly5(n)*78^(length(poly5) -n); 
% end

extrap2015_2 = polyval(poly2, 78);
extrap2015_5 = polyval(poly5, 78);
extrap2015_8 = polyval(poly8, 78);

% testextrap = polyval(poly8, t);
% plot(t, testextrap)
% hold on
% plot(t, salmon_data)

A13 = [extrap2015_2;extrap2015_5;extrap2015_8];
save A13.dat A13 -ASCII

%% Exercise 2.d.)
t_coarse = transpose(1:4:77);
salmon_coarse = zeros(length(t_coarse),1);
counter = 1;
for n = 1:4:77
    salmon_coarse(counter) = salmon_data(n);
    counter = counter + 1;
end

save A14.dat salmon_coarse -ASCII

%% Exercise 2.e.)

neighborInterp = interp1(t_coarse, salmon_coarse, t, 'nearest');
linearInterp = interp1(t_coarse, salmon_coarse, t, 'linear');
cubicInterp = interp1(t_coarse, salmon_coarse, t, 'cubic');
splineInterp = interp1(t_coarse, salmon_coarse, t, 'spline');
% plot(t, neighborInterp)
% hold on
% plot(t, linearInterp, 'b')
% plot(t, cubicInterp, 'g')
% plot(t, splineInterp, 'k')
% plot(t, salmon_data, 'm--')
save A15.dat neighborInterp -ASCII
save A16.dat linearInterp -ASCII
save A17.dat cubicInterp -ASCII
save A18.dat splineInterp -ASCII

%% Exercise 2.f.)
interpHolder = [neighborInterp, linearInterp, cubicInterp, splineInterp];
A19 = zeros(4,1);

for n = 1:4
    A19(n) = sqrt(1/77 * sum((salmon_data - interpHolder(:,n)).^2));
end

save A19.dat A19 -ASCII

%% Exercise 3.a.)
clc;clear all;close all;

theta = 0:.05*pi:2*pi;
r = 0:.5:20;
[TH, R] = meshgrid(theta, r);
elecOrbit(TH, R)
% [X, Y] = pol2cart(TH, R);
% surf(X, Y, elecOrbit(TH, R))

function F = elecOrbit(TH, R)
    e = exp(1);
    F = abs(sqrt(2)/81*pi*(6*R - R^2)*e^(-R/3)*cos(TH))^2;
end

