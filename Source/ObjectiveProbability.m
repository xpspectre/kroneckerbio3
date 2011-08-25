function p = ObjectiveProbability(m, con, obj, opts)
%ObjectiveProbability evaluates the likelihood of a set of 
%   information-theory-based objective functions
%
%   p = ObjectiveProbability(m, con, obj, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The information-theory-based objective structures defining the
%       probability distribution to be evaluated.
%   opts: [ options struct scalar ]
%       Optional
%       .UseModelICs [ logical scalar {false} ]
%           Indicates that the model's initial conditions should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Which kinetic parameters the gradient will be calculated on.
%           Has an effect only if the objective function is sensitive to
%           this parameter
%       .UseICs [ logical matrix nx by nCon | logical vector nx |
%                 positive integer vector {[]} ]
%           Which initial conditions the gradient will be calculated on Has
%           an effect only if the objective function is sensitive to this
%           parameter
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Which input control parameters the gradient will be calculated
%           on Has an effect only if the objective function is sensitive to
%           this parameter
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
%       p: [ real nonnegative scalar ]
%           The product of all objective function probabilities

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:ObjectiveProbability:MoreThanOneModel', 'The model structure must be scalar')

% Options
defaultOpts.UseModelICs    = true;
defaultOpts.UseModelInputs = false;
defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];

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
nObj = size(obj,1);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a matrix of logical indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTx + nTk;

% Refresh conditions and objectives
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon);

%% Probability
% Initialize variables
p = 1;
intOpts = opts;

for iCon = 1:nCon
    if verbose; disp(['Integrating system for ' con(iCon).Name '...']); end
    % If opts.UseModelICs is false, get the right UseICs
    if ~opts.UseModelICs
        intOpts.UseICs = opts.UseICs(:,iCon);
    end
    
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    tGet = opts.tGet{iCon};
        
    % Integrate
    if opts.continuous(iCon) && opts.complex(iCon)
        sol = integrateObj(m, con(iCon), obj(:,iCon), intOpts);
    elseif opts.complex(iCon)
        sol = integrateSys(m, con(iCon), intOpts);
    elseif opts.continuous(iCon)
        sol = integrateObjSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
    else
        sol = integrateSysSelect(m, con(iCon), tGet, intOpts);
    end
    
    % Probability
    for iObj = 1:nObj
        p = p * obj(iObj,iCon).p(sol);
    end
end
