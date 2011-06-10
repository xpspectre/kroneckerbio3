function D = FiniteObjectiveGradient(m, con, obj, opts)
%FiniteObjectiveGradient evaluates the gradient of a set of objective
%   functions
%
%   D = FiniteObjectiveGradient(m, con, obj, opts)
%
%   This function is identical to ObjectiveGradient except that the
%   gradient is calculated differently. This is useful for testing that
%   custom-made objective functions have been coded correctly and that
%   integration tolerances have been set appropriately.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%       .UseModelICs [ logical scalar {false} ]
%           Indicates that the model's initial conditions should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value
%       .Normalized [ logical scalar {true} ]
%           Indicates if the gradient should be computed in log parameters
%           space
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method.
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%       D: [ real vector nT ]
%           The sum of all objective function gradients

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inVuts
assert(nargin >= 3, 'KroneckerBio:FiniteObjectiveGradient:AtLeastThreeInputs', 'FiniteObjectiveSensitivity requires at least 3 input arguments.')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:ObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Options
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;
defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));
defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.Verbose        = 0;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
nCon = numel(con);
nObj = size(obj, 1);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);
useParamsInd = find(opts.UseParams);

% Ensure UseICs is a matrix of logical indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);
useICsInd = find(opts.UseICs);

nT = nTk + nTx;

% Store starting parameter sets
k = m.k;

if opts.UseModelICs
    x0 = m.x0;
else
    x0 = zeros(nx, nCon);
    for i = 1:nCon
        x0(:,i) = con(i).x0;
    end
end

% Standardize structures
con = pastestruct(Uzero(m), con);
obj = pastestruct(Gzero(m), obj);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon);

%% Loop through conditions
D = zeros(nT,1);

% Initial value
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
        Ti = k(useParamsInd(iT));
        if opts.Normalized
            diff = Ti * 1e-8;
        else
            diff = 1e-8;
        end
        kup(useParamsInd(iT)) = Ti + diff;
        kdown(useParamsInd(iT)) = Ti - diff;
    else
        Ti = x0(useICsInd(iT-nTk));
        if opts.Normalized
            diff = Ti * 1e-8;
        else
            diff = 1e-8;
        end
        x0up(useICsInd(iT-nVP)) = Ti + diff;
        x0down(useICsInd(iT-nVP)) = Ti - diff;
    end
    
    % Upper parameter set
    % Update with new parameter set
    if opts.UseModelICs
        m = m.Update(kup, x0up, m.q);
    else
        m = m.Update(kup, m.x0, m.q);
        for iCon = 1:nCon
            con(iCon) = pastestruct(Uzero(m), con(iCon).Update(x0up(:,iCon), con(iCon).q));
        end
    end
    
    for iCon = 1:nCon
        for iObj = 1:nObj
            obj(iObj,iCon) = pastestruct(Gzero(m), obj(iObj,iCon).Update(m, con(iCon), opts.UseParams, opts.UseICs, opts.UseControls));
        end
    end
    
    % Finitely different goal
    Gup = computeObj(m, con, obj, opts);
    
    % Lower parameter set
    % Update with new parameter set
    if opts.UseModelICs
        m = m.Update(kdown, x0down, m.q);
    else
        m = m.Update(kdown, m.x0, m.q);
        for iCon = 1:nCon
            con(iCon) = pastestruct(Uzero(m), con(iCon).Update(x0down(:,iCon), con(iCon).q));
        end
    end
    
    for iCon = 1:nCon
        for iObj = 1:nObj
            obj(iObj,iCon) = pastestruct(Gzero(m), obj(iObj,iCon).Update(m, con(iCon), opts.UseParams, opts.UseICs, opts.UseControls));
        end
    end
    
    % Finitely different goal
    Gdown = computeObj(m, con, obj, opts);
    
    % Compute D
    if opts.Normalized
        D(iT) = Ti * ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
    else
        D(iT) = ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
    end
end
