function [D All] = FiniteObjectiveGradient(m, con, obj, opts)
%FINITEOBJECTIVEGRADIENT Compute the gradient of an objective function by
%   finite difference steps.
%
%   Given a model and experimental conditions, this function computes the
%   derivative of an objective function with respect to the parameters by
%   taking small steps in the parameters and evaluating the objective
%   function.
%
%   Mathematically: dG/dp = sum(dG/dx(t) * dx/dp, t)
%
%   [D All] = FiniteObjectiveGradient(m, con, obj, opts)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A single experimental condition under which the model
%              will be simulated. (Vector is not supported, yet.)
%       obj  - A vector of objective functions that will be evaluated
%       opts - Scalar structure of options
%           UseModelICs - Logical scalar indicating that the model's
%                         initial conditions should be used instead of
%                         those of the experimental conditions.
%                         Default = true
%           UseParams   - Vector of indexes indicating the rate constants
%                         whose sensitivities will be considered.
%                         Default = 1:m.nP
%           UseICs      - Vector of indexes indicating the initial
%                         concentrations whose sensitivities will be
%                         considered. Default = []
%           Normalized  - Scalar logical indicating that the gradient
%                         should be normalized with respect to the
%                         parameter values. Default = true
%           ObjWeights  - A matrix of size(obj) indicating a post-
%                         evaluation weight on each objective function in
%                         terms of how much it will contribute to the final
%                         value. Default = ones(size(obj))
%           RelTol      - Numerical scalar for the relative tolerance of 
%                         the integration. Default = 1e-6
%           AbsTol      - Cell vector of vectors, numerical scalar, or 
%                         numerical vector for the absolute tolerance of
%                         the integration. Default = 1e-9 
%           Verbose     - Numerical scalar for how much progress 
%                         information should be displayed. Default = 0
%
%   Outputs:
%       D   - The sum of all objective function gradients
%
%       All - A matrix size(obj) for the values for each objective function

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inVuts
assert(nargin >= 3, 'KroneckerBio:FiniteObjectiveGradient:AtLeastThreeInputs', 'FiniteObjectiveSensitivity requires at least 3 input arguments.')
if nargin < 4
	opts = [];
end

assert(isscalar(m), 'KroneckerBio:FiniteObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Options
defaultOpts.UseModelICs = true;
defaultOpts.UseParams	= 1:m.nP;
defaultOpts.UseICs		= [];
defaultOpts.Normalized  = true;
defaultOpts.ObjWeights  = ones(size(obj));
defaultOpts.RelTol      = NaN;
defaultOpts.AbsTol      = NaN;
defaultOpts.Verbose     = 0;

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

opts = mergestruct(defaultOpts, opts);

% Constants
nx = m.nX;
nk = m.nP;
nCon = numel(con);
nObj = size(obj, 2);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);
UseParamsInd = find(opts.UseParams);

% Ensure UseICs is a matrix of logical indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);
UseICsInd = find(opts.UseICs);

nT = nTk + nTx;

% Store starting parameter sets
k = m.p;

if opts.UseModelICs
    x0 = m.ic;
else
    x0 = zeros(nX, nCon);
    for i = 1:nCon
        x0(:,i) = con(i).ic;
    end
end

% Standardize structures
con = pastestruct(Uzero(m), con);
obj = pastestruct(Gzero(m), obj);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);

%% Tolerances
% RelTol
if isempty(opts.RelTol) || isnan(opts.RelTol)
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon, opts.UseAdjoint, opts.UseParams, opts.UseICs, opts.UseModelICs);

%% Loop through conditions
D = zeros(nT,1);

% Initial value
if nargout <= 1
    if verbose; fprintf('Initial round\n'); end
    G = computeObj(m, con, obj, opts);
    
    for iT = 1:nT
        if verbose; fprintf('Step %d of %d\n', iT, nT); end
        
        % Set baseline parameters
        kup = k;
        x0up = x0;
        kdown = k;
        x0down = x0;
        
        % Change current parameter by finite amount
        if iT <= nTk
            Ti = k(UseParamsInd(iT));
            if opts.Normalized
                diff = Ti * 1e-8;
            else
                diff = 1e-8;
            end
            kup(UseParamsInd(iT)) = Ti + diff;
            kdown(UseParamsInd(iT)) = Ti - diff;
        else
            Ti = x0(UseICsInd(iT-nTk));
            if opts.Normalized
                diff = Ti * 1e-8;
            else
                diff = 1e-8;
            end
            x0up(UseICsInd(iT-nVP)) = Ti + diff;
            x0down(UseICsInd(iT-nVP)) = Ti - diff;
        end
        
        % Upper parameter set
        % Update with new parameter set
        newCon = Uzero(nCon);
        if opts.UseModelICs
            m = m.update(kup, x0up);
            for iCon = 1:nCon
                newCon(iCon) = pastestruct(Uzero(m), con(iCon).update(m));
            end
        else
            m = m.update(kup, m.ic);
            for iCon = 1:nCon
                newCon(iCon) = pastestruct(Uzero(m), con(iCon).update(m, x0up(:,iCon)));
            end
        end
        con = newCon;
        
        newObj = Gzero([nCon, nObj]);
        for iCon = 1:nCon
            for iObj = 1:nObj
                newObj(iCon,iObj) = pastestruct(Gzero(m), obj(iCon,iObj).update(m));
            end
        end
        obj = newObj;
        
        % Finitely different goal
        Gup = computeObj(m, con, obj, opts);
        
        % Lower parameter set
        % Update with new parameter set
        newCon = Uzero(nCon);
        if opts.UseModelICs
            m = m.update(kdown, x0down);
            for iCon = 1:nCon
                newCon(iCon) = pastestruct(Uzero(m), con(iCon).update(m));
            end
        else
            m = m.update(kdown, m.ic);
            for iCon = 1:nCon
                newCon(iCon) = pastestruct(Uzero(m), con(iCon).update(m, x0down(:,iCon)));
            end
        end
        con = newCon;
        
        newObj = Gzero([nCon, nObj]);
        for iCon = 1:nCon
            for iObj = 1:nObj
                newObj(iCon,iObj) = pastestruct(Gzero(m), obj(iCon,iObj).update(m));
            end
        end
        obj = newObj;
        
        % Finitely different goal
        Gdown = computeObj(m, con, obj, opts);

        % Compute D
        if opts.Normalized
            D(iT) = Ti * ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
        else
            D(iT) = ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
        end
    end
    
else %nargout == 2
    error('All will not work until computeObjAll is written')
    % Initialize All array
    All = cell(nCon,nObj);
    
    for iCon = 1:nCon
        for iObj = 1:nObj
            if verbose; fprintf('Initial round\n'); end
            G = computeObj(m, con(iCon), obj(iCon,iObj), opts);
            
            % Finite stepping
            curD = zeros(nV,1);
            
            for i = 1:nT
                if verbose; fprintf('Step %d of %d\n', i, nV); end
                
                % Set baseline parameters
                kup = m.p;
                kdown = m.p;
                
                if opts.UseModelICs
                    xup = m.ic;
                    xdown = m.ic;
                else
                    xup = con{iCon}.ic;
                    xdown = con{iCon}.ic;
                end
                
                % Change current parameter by finite amount
                if i <= nVP
                    pi = kup(opts.UseParams(i));
                    diff = kup(opts.UseParams(i)) * 1e-8;
                    kup(opts.UseParams(i)) = pi + diff;
                    kdown(opts.UseParams(i)) = pi - diff;
                else
                    pi = xup(opts.UseICs(i-nVP));
                    diff = xup(opts.UseICs(i)) * 1e-8;
                    xup(opts.UseICs(i-nVP)) = pi + diff;
                    xdown(opts.UseICs(i-nVP)) = pi - diff;
                end
                
                % Run models to get goal function
                if opts.UseModelICs
                    mtemp = m.update(kup, xup);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Gup = computeObj(mtemp, con(iCon), objtemp, opts);
                    mtemp = m.update(kdown, xdown);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Gdown = computeObj(mtemp, con(iCon), objtemp, opts);
                else
                    mtemp = m.update(kup, m.ic);
                    contemp = con(iCon).update(m, xup);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Gup = computeObj(mtemp, contemp, objtemp, opts);
                    mtemp = m.update(kdown, m.ic);
                    contemp = con(iCon).update(m, xdown);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Gdown = computeObj(mtemp, contemp, objtemp, opts);
                end
                
                % Compute D
                if opts.Normalized
                    curD(i) = pi * ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
                else
                    curD(i) = ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
                end
            end
            
            D = D + curD;
            
            All{iCon, iObj} = curD;
        end
    end
end