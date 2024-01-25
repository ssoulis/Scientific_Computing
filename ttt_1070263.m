function C = ttt_myid_1070263(A, B, mode)
    % Check input types and dimensions
    if ~ismatrix(A) || ~ismatrix(B)
        error('Input error: A and B must be matrices.');
    end
    
    if nargin == 2
        % Perform outer product if mode is not specified
        C = tensor(bsxfun(@times, A, permute(B, [2, 1])), [size(A), size(B)]);
    elseif nargin == 3 && ischar(mode) && strcmpi(mode, 'all')
        % Perform inner product if mode is 'all'
        C = sum(A(:) .* B(:));
    else
        error('Input error: Invalid mode for ttt_myid.');
    end
end
