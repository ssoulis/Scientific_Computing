function Y = ttv_myid_1070263(X, V, N)
    % Check input types and dimensions
    if ~ismatrix(X) || ~iscolumn(V)
        error('Input error: X must be a matrix, and V must be a column vector.');
    end
    
    if numel(size(X)) < N || size(X, N) ~= numel(V)
        error('Input error: Dimensions are not appropriate for the function.');
    end
    
    % Perform tropical multiplication
    Y = X;
    for i = 1:numel(V)
        Y = min(Y, circshift(X, [0 0 0 -i]), [], N);
    end
end
