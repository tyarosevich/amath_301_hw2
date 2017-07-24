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
save A1.dat A1 -ASCII

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
save A4.dat A4 -ASCII

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
save A5.dat x -ASCII
save A6.dat iterations -ASCII


%% Exercise 2
clc;close all;clear all;
load salmon_data.csv
t = (1:length(salmon_data)).';
Q_11 = sum(t.^2); Q_12 = sum(t); Q_21 = Q_12; Q_22 = 77;
R_11 = sum(t.*salmon_data); R_21 = sum(salmon_data);
Q = [Q_11 Q_12;Q_21 Q_22];
R = [R_11; R_21];
P = Q\R;
%y_coeffs = P(1,1)*t + P(2,1);
% plot(t, salmon_data)
%hold on
%plot(t, y_coeffs, 'g')
%ptest = polyfit(t, salmon_data, 1);
%y_coeffs2 =  p(1,1)*t + p(1,2)
%plot(t, y_coeffs2, 'b--o')
save A7.dat Q -ASCII
save A8.dat R -ASCII
save A9.dat P -ASCII


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

% plot(t, ycoeffs2, 'g')
% hold on
% plot(t, ycoeffs5, 'b')
% plot(t, ycoeffs8, 'k')
% plot(t, salmon_data)

%% Exercise 2. c.)


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
e = exp(1);

theta = 0:.05*pi:2*pi;
r = 0:.5:20;
[TH, R] = meshgrid(theta, r);
f = [TH R];
F = elecOrbit(f);
[X, Y] = pol2cart(TH, R);
surf(X, Y, F)
xlabel('theta');ylabel('rad');zlabel('height');
hold on
save A20.dat X -ASCII
save A21.dat Y -ASCII
save A22.dat F -ASCII


%% Exercise 3.b.)

A23 = 0; A24 = 0;A25 = 0;A26 = 0;

A23 = transpose(fminsearch('elecOrbit2', [0; 1]));
A24 = transpose(fminsearch('elecOrbit2', [0; 10]));
A25 = transpose(fminsearch('elecOrbit2', [pi; 1]));
A26 = transpose(fminsearch('elecOrbit2', [pi; 10]));
save A23.dat A23 -ASCII
save A24.dat A24 -ASCII
save A25.dat A25 -ASCII
save A26.dat A26 -ASCII


%% Exercise 3.c.)
g1 = [0; 1];g2 = [0; 10];g3 = [pi; 1];g4 = [pi; 10];
A27 = newton_2D(g1);
A28 = newton_2D(g2);
A29 = newton_2D(g3);
A30 = newton_2D(g4);

save A27.dat A27 -ASCII
save A28.dat A28 -ASCII
save A29.dat A29 -ASCII
save A30.dat A30 -ASCII

% tol = 1.e-4;

% delta = [2; 2];
% x = g3;x0 = x;
% iterations = 0;
% while norm(delta, 2) > norm(x0*tol, 2)
%     iterations = iterations + 1;
%     G = first_grad([x0]);
%     DG = second_grad([x0]);
%     delta = DG\G;
%     x0 = x0 - delta;
% end
% 
% iterations
%% Exercise 4.a.)
clc;close all;clear all;

f = [40 50 75];
A31 = transpose(f);
save A31.dat A31 -ASCII

%% Exercise 4.b.)
A = [.1 .23 .31;
    .22 .25 .38;
    -.1 -.23 -.31;
    .22+.15 .25+.17 .38+.27;
    0 0 .31+.38+.27;
    -1 0 0;
    0 -1 0;
    0 0 -1];
b = [5; 7; -4; 15; 2; -3; -4; -2];

save A32.dat A -ASCII
save A33.dat b -ASCII
%% Exercise 4.c.)

x = linprog(f, A, b);
x_real = round(x);

save A34.dat x_real -ASCII




