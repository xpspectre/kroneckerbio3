function [pmy pyTm F m] = TopologyProbability(m, con, obj, opts, F)
% TopologyProbability Compute the relative probability that each member of
%   a set of topologies is true according to a set of
%   information-theory-based objective functions
%
%   [pmy pyTm F m] = TopologyProbability(m, con, obj, opts)

% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:TopologyProbability:TooFewInputs', 'TopologyProbability requires at least 3 input arguments')
if nargin < 5
    F = [];
    if nargin < 4
        opts = [];
    end
end

% Constants
nTop = numel(m);
nCon = size(con,1);
nObj = size(obj,1);

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = NaN;
defaultOpts.UseICs         = NaN;
defaultOpts.UseControls    = NaN;

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;

defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.NeedFit        = true; % Fit the parameters

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = zeros(nTop,1);
nk = zeros(nTop,1);
for iTop = 1:nTop
    nx(iTop) = m(iTop).nx;
    nk(iTop) = m(iTop).nk;
end

%% Standardize structures
% Experimental conditions
assert((size(con,2) == 1 && (opts.UseModelICs || all(nx == nx(1)))) || size(con,2) == nTop, 'KroneckerBio:TopologyProbability:ConSize', 'Second dimension of input "con" must be equal to numel(m) or 1 if opts.UseModelICs is false or every model has the same number of species')
if size(con,2) == 1
    con = repmat(con, 1,nTop);
end

% Objective functions
assert(size(obj,3) == 1 || size(obj,3) == nTop, 'KroneckerBio:TopologyProbability:ObjSize', 'Third dimension of input "obj" must be equal to numel(m) or 1')
if size(obj,3) == 1
    obj = repmat(obj, [1,1,nTop]);
end

%% Active Parameters
% Default UseParams is all kinetic parameters
if isnumeric(opts.UseParams) && (isempty(opts.UseParams) || any(isnan(opts.UseParams)))
    opts.UseParams = cell(nTop,1);
    for iTop = 1:nTop
        opts.UseParams{iTop} = 1:m(iTop).nk;
    end
end

% Default UseICs is no initial condition parameters
if isnumeric(opts.UseICs) && (isempty(opts.UseICs) || any(isnan(opts.UseICs)))
    opts.UseICs = cell(nTop,1);
    for iTop = 1:nTop
        opts.UseICs{iTop} = [];
    end
end

% Default UseControls is no controls
if isnumeric(opts.UseControls) && (isempty(opts.UseControls) || any(isnan(opts.UseControls)))
    opts.UseControls = cell(nCon,nTop);
    for iTop = 1:nTop
        opts.UseControls{iTop} = [];
    end
end

% Ensure UseParams is vector of logical indexes within a cell array
nTk = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseParams{iTop} nTk(iTop)] = fixUseParams(opts.UseParams{iTop}, nk(iTop));
end

% Ensure UseICs is a matrix of logical indexes within a cell array
nTx = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseICs{iTop} nTx(iTop)] = fixUseICs(opts.UseICs{iTop}, opts.UseModelICs, nx(iTop), nCon);
end

% Ensure UseControls is a cell array of logical vectors
nTq = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseControls{iTop} nTq(iTop)] = fixUseControls(opts.UseControls{iTop}, opts.UseModelInputs, nCon, m(iTop).nq, cat(1,con(:,iTop).nq));
end

nT = nTk + nTx + nTq;

%% Refresh structures
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
if isnumeric(opts.AbsTol)
    opts.AbsTol = repmat({vec(opts.AbsTol)}, 1,nTop);
end
if isstruct(opts.AbsTol)
    opts.AbsTol = vec(opts.AbsTol).';
end
tempAbsTol = cell(nTop,1);
for iTop = 1:nTop
    % Distribute AbsTol because TopologyProbability never uses it directly
    tempAbsTol{iTop} = opts.AbsTol(:,iTop);
end
opts.AbsTol = tempAbsTol;

%% Fix bounds
if isnumeric(opts.LowerBound)
    opts.LowerBound = repmat({opts.LowerBound}, nTop,1);
end
if isnumeric(opts.UpperBound)
    opts.UpperBound = repmat({opts.UpperBound}, nTop,1);
end
for iTop = 1:nTop
    opts.LowerBound{iTop} = fixBounds(opts.LowerBound{iTop}, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
    opts.UpperBound{iTop} = fixBounds(opts.UpperBound{iTop}, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
end

%% Distribute options for topology specific functions
optsTop = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsTop(iTop).UseParams   = opts.UseParams{iTop};
    optsTop(iTop).UseICs      = opts.UseICs{iTop};
    if ~opts.UseModelICs
        optsTop(iTop).UseICs  = optsTop(iTop).UseICs(:,1:nCon);
    end
    optsTop(iTop).UseControls = opts.UseControls{iTop};
    optsTop(iTop).AbsTol      = opts.AbsTol{iTop};
    optsTop(iTop).LowerBound  = opts.LowerBound{iTop};
    optsTop(iTop).UpperBound  = opts.UpperBound{iTop};
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).Name ' to objectives...\n']); end
        m(iTop) = FitObjective(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    end
end

%% Information
% Compute the information if not provided
if isempty(F)
    F = cell(nTop,1);
    for iTop = 1:nTop
        F{iTop} = ObjectiveInformation(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    end
end

% Extract eigenvalues
lambda = cell(nTop,1);
for iTop = 1:nTop
    % Eigendecompose
    lambda{iTop} = infoeig(F{iTop});
    
    % Bend flat directions
    lambda{iTop}(lambda{iTop} < 1e-16) = 1e-16;
end

%% Compute p_y|m for each topology
% p_y|m will be computed in log space in order to give more digits to the
% exponent for computing this number, which is often too small to be
% represented by float64.
pym = zeros(nTop,1); % In log space
pyTm = ones(nTop,1);
for iTop = 1:nTop
    if verbose; fprintf(['Computing probability of ' m(iTop).Name '...\n']); end
    % Likelihood
    pyTm(iTop) = ObjectiveProbability(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    
    % Equation for linearization about maximum a posteriori
    %pym(iTop) = pyTm(iTop) * (2*pi).^(nT(iTop)/2) * det(F{iTop}).^(-1/2); % Without log space
    pym(iTop) = log(pyTm(iTop)) + (nT(iTop)/2) * log(2*pi) + -1/2 * sum(log(lambda{iTop}));
end

%% Compute p_m|y for the set
% Apply topology prior
%pmy = pym .* opts.PriorTopology; % Without log space
pmy = pym + log(opts.PriorTopology); % In log space

% Rescale p_y|m
% Since the entire distribution is normalized at the end, this operation
% has no mathematical effect on the answer. However, it ensures that all
% the probabilities are representable by float64 when numbers are
% exponentiated back into regular space.
pmy = pmy - max(pmy);

% Return to regular space
pmy = exp(pmy);

% Normalize 
pmy = pmy ./ sum(pmy);
