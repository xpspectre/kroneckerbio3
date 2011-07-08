function [lambda V] = infoeig(F)

% Size
nT = size(F,1);
diagF = diag(F);

% Find specials: zeros and infs
zeroInd = find(all(F == 0, 2));
n0 = numel(zeroInd);
negInfInd  = find(diagF == -inf);
posInfInd  = find(diagF == inf);
nNegInf = numel(posInfInd);
nPosInf = numel(negInfInd);

% Eig non-special portion of F
invertableInd = true(nT,1);
invertableInd([zeroInd;negInfInd;posInfInd]) = false;

if nargout <= 1
    lambda = eig(F(invertableInd,invertableInd));
else
    [Vplain lambda] = eig(F(invertableInd,invertableInd));
    lambda = diag(lambda);
end

% Recombine results
positiveLambda = (lambda >= 0);
nNegativeLambda = sum(~positiveLambda);
nPositiveLambda = sum(positiveLambda);
lambda = [inf(nNegInf,1)*-1;
          lambda(~positiveLambda);
          zeros(n0,1);
          lambda(positiveLambda);
          inf(nPosInf,1)];

if nargout >=2
    V = zeros(nT,nT);
    V(sub2ind([nT,nT], negInfInd,(1:nNegInf)')) = 1;% -inf
    V(invertableInd, nNegInf+1:nNegInf+nNegativeLambda) = Vplain(:,~positiveLambda); % negative lambda
    V(sub2ind([nT,nT], zeroInd,(nNegInf+nNegativeLambda+1:nNegInf+nNegativeLambda+n0)')) = 1; % 0
    V(invertableInd,  nNegInf+nNegativeLambda+n0+1:nNegInf+nNegativeLambda+n0+nPositiveLambda) = Vplain(:,positiveLambda); % positive lambda
    V(sub2ind([nT,nT], posInfInd,(nNegInf+nNegativeLambda+n0+nPositiveLambda+1:nT)')) = 1; % +inf
end