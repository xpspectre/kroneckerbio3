function A = spermute132(A, dim, final)

% Extract dimensions
d1 = dim(1);
d2 = dim(2);
d3 = dim(3);

% Put the swapping columns together
A = reshape(A, d1,d2*d3); % 1_23

% Permute
A = A(:, reshape(1:(d2*d3), d2,d3)'); % 1_32

% Reshape to desired dimensions
A = reshape(A, final);