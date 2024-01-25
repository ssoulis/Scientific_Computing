% Initialisation
clear; tol = 1e-8; nd = 3; rng(0); err = zeros(1, nd + 2);
ndim = [2, 3, 4]; Atemp = randi(5, ndim); Btemp = randi(4, ndim);
X = randi([-1, 1], max(ndim), 1); A = tensor(Atemp); B = tensor(Btemp);

% Test ttv_myid
try
    for k = 1:nd
        % Compare with TT command
        err(k) = norm(ttv_1070263(A, X(1:ndim(k), 1), k) - double(ttv(A, X(1:ndim(k), 1), k)));
    end
    assert(max(err) < tol, 'ttv modal multiplication fails');
catch ME1
    ME1.message
end

% Test ttm_myid
try
    % Compare with TT command
    err(nd + 1) = norm(tensor(ttm_1070263(A, X(1:ndim(k), 1), k)) - ttm(A, X(1:ndim(k), 1), k));
    assert(err(nd + 1) < tol, 'ttm modal multiplication fails');
catch ME2
    ME2.message
end

% Test ttt_myid outer product
try
    % Compare with TT command
    C_outer = ttt_1070263(A, B);
    assert(isequal(C_outer, ttt(A, B)), 'ttt outer product fails');
catch ME3
    ME3.message
end

% Test ttt_myid inner product
try
    % Compare with TT command
    t_inner = ttt_1070263(A, B, 'all');
    assert(abs(t_inner - double(ttt(A, B))) < tol, 'ttt inner product fails');
catch ME4
    ME4.message
end
