function [L, U] = customPreconditioner(A)
    [L, U] = ichol(A, struct('type', 'ict', 'droptol', 1e-2, 'shape', 'lower'));
end
