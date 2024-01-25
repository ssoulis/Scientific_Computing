    function [L, U] = incompleteCholesky(A)
        % Perform incomplete Cholesky factorization IC(0)
        R = cholinc(A, '0');
        L = R';
        U = R;
    end
