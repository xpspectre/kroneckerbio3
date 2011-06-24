function [lambda V] = infoeig(F)

% Size
nT = size(F,1);

% Find specials: zeros and infs
zeroInd = find(diag(F) == 0);
n0 = numel(zeroInd);
infInd  = find(diag(F) == inf);
ninf = numel(infInd);

% Eig non-special portion of F
invertableInd = true(nT,1);
invertableInd([zeroInd;infInd]) = false;

if nargout < 2
    lambda = eig(F(invertableInd,invertableInd));
else
    [V lambda] = eig(F(invertableInd,invertableInd));
    lambda = diag(lambda);
end

% Recombine results
lambda = [zeros(n0,1);lambda;inf(ninf,1)];
if nargout >=2
    Vold = V;
    V = zeros(nT,nT);
    V(zeroInd,1:n0) = 1;
    V(invertableInd,n0+1:nT-ninf) = Vold;
    V(zeroInd,nT-ninf+1:nT) = 1;
end