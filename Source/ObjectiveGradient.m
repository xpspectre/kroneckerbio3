function D = ObjectiveGradient(m, con, obj, opts)
%ObjectiveGradient evaluates the gradient of a set of objective functions
%
%   D = ObjectiveGradient(m, con, obj, opts)
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
% Clean up inputs
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

% Constants
nx = m.nx;
nk = m.nk;
nCon = numel(con);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a matrix of linear indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Add missing fields to structure
con = pastestruct(Uzero(m), con);
obj = pastestruct(Gzero(m), obj);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, opts.UseAdjoint, opts.UseParams, opts.UseICs, opts.UseModelICs);

%% Run appropriate integration
if opts.UseAdjoint
    [unused D] = integrateObjSensAdj(m, con, obj, opts);
else
    [unused D] = computeObjSens(m, con, obj, opts);
end

%% Normalization
if opts.Normalized
    % Extract parameter set
    kAll = m.k;
    
    if opts.UseModelICs
        x0All = m.x0;
    else
        x0All = zeros(nx, nCon);
        for i = 1:nCon
            x0All(:,i) = con(i).x0;
        end
    end
    
    T = [kAll(opts.UseParams); vec(x0All(opts.UseICs))];
    
    % Normalize
    D = D .* T;
end