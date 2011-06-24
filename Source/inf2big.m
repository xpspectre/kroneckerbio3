function A = inf2big(A)
%inf2big converts infinities to just big numbers
%
%   A = inf2big(A)

A(A == -inf) = -1e6;
A(A == inf)  =  1e6;