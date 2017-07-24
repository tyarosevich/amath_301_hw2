function [x] = newton_2D(x0)
    tol = 1.e-4;
    delta = [2; 2];
    iterations = 0;
    x = x0;
    while norm(delta, 2) > norm(x*tol, 2)
        x0=x;
        iterations = iterations + 1;
        G = first_grad([x0]);
        DG = second_grad([x0]);
        delta = DG\G;
        x = x0 - delta;
    end