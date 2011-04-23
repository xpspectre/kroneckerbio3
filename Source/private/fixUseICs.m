function [UseICs, nTx] = fixUseICs(UseICs, UseModelICs, nx, nCon)
%FIXUSEICS A helper function for functions that support variant initial
%   concentrations in experiments. It converts a variety of ways to
%   specific the active initial conditions and standardizes them into
%   either a vector (use model ICs) or matrix (use experiment ICs) of
%   logicals indicating the active parameters.
%
%   [UseICs, nTx] = fixUseICs(UseICs, UseModelICs, nx, nCon)
%
%   Inputs
%       UseICs - Can be any of the following:
%           If UseModelICs = true
%               1) vector of linear indexes
%               2) vector of logical indexes length of nx
%           If UseModelICs = false
%               1) vector of linear indexes into nx, assumed same for all
%                  conditions
%               2) vector of logical indexes length(nx) to indicate that
%                  all conditions have the same active parameters
%               3) matrix of logical indexes size nx by nCon)
%       UseModelICs - Scalar logical indicating that the model initial
%                     conditions are active and constant for all 
%                     experiments
%       nx - Scalar natural number for how many species are in this model
%       nCon - Scalar natural number for how many experiments are used in
%              this function
%
%   Outputs
%       UseICs - If UseModelICs = true, then UseICs will be a logical
%           column vector length of nx. Otherwise, it will be a logical
%           matrix sizw nx by nCon.
%       nTx - Number of active IC parameters

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

assert(islogical(UseICs) || isempty(UseICs) || (all(isinteger(UseICs)) && all(UseICs >= 1)), 'KroneckerBio:UseICs:Integer', 'Some entries in UseICs are not natural numbers.')
assert(islogical(UseICs) || isempty(UseICs) || (max(UseICs) <= nx), 'KroneckerBio:UseICs:InvalidIndex', 'An entry in UseICs was larger than the number of species in the model.')

if UseModelICs
    if islogical(UseICs)
        UseICs = vec(UseICs);
    else
        temp = UseICs;
        UseICs = zeros(nx, 1);
        UseICs(temp) = 1;
        UseICs = logical(UseICs);
    end
else%if use ICs on conditions
    if islogical(UseICs)
        if numel(UseICs) == nx
            % Parameter ICs are same for all conditions
            temp = UseICs;
            UseICs = zeros(nx, nCon);
            for iCon = 1:nCon
                UseICs(:,iCon) = temp;
            end
            UseICs = logical(UseICs);
        else%if numel(UseICs) == nx*nCon
            % Parameter ICs are provided for all conditions. Format is
            % already correct.
        end
    else%if isnumeric
        temp = UseICs;
        UseICs = false(nx,1);
        UseICs(temp) = true;
        UseICs = repmat(UseICs, 1,nCon);
    end
end
nTx = sum(sum(UseICs));
