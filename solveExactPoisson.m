function xsol = solveExactPoisson(X, Y)
    % Compute the exact solution for the Poisson equation
    xsol = sin(pi * X) .* cos(pi * Y);
end