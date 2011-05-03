function [pmy pyTm F m] = TopologyProbability(m, con, obj, opts)
% Clean-up inputs
if nargin < 4
    opts = [];
end

% Constants
nTop = numel(m);
nCon = size(con,1);
nObj = size(obj,1);

% Options
defaultOpts.RelTol         = NaN; % 1e-6
defaultOpts.AbsTol         = NaN; % in fixAbsTol
defaultOpts.UseModelICs    = true;
defaultOpts.UseAdjoint     = false;

defaultOpts.UseParams      = NaN; % All kinetic parameters
defaultOpts.UseICs         = NaN; % No initial conditions

defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.NeedFit        = true; % Fit the parameters

defaultOpts.Verbose        = 0;

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
con = refreshCon(m, con);

% Objective functions
assert(size(obj,3) == 1 || size(obj,3) == nTop, 'KroneckerBio:TopologyProbability:ObjSize', 'Third dimension of input "obj" must be equal to numel(m) or 1')
if size(obj,3) == 1
    obj = repmat(obj, [1,1,nTop]);
end
obj = refreshObj(m, con, obj);

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

nT = nTk + nTx;

%% Tolerances
if isnan(opts.RelTol) || isempty(opts.RelTol) 
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
if isnumeric(opts.AbsTol)
    opts.AbsTol = repmat(vec(opts.AbsTol), 1,nTop);
end
tempAbsTol = cell(nCon,nTop);
for iTop = 1:nTop
    tempAbsTol(:,iTop) = fixAbsTol(opts.AbsTol(:,iTop), 2, false(nCon,1), nx(iTop), nCon, opts.UseAdjoint, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseModelICs);
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
    opts.LowerBound{iTop} = fixBounds(opts.LowerBound{iTop}, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseModelICs);
    opts.UpperBound{iTop} = fixBounds(opts.UpperBound{iTop}, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseModelICs);
end

%% Distribute options for topology specific functions
optsTop = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsTop(iTop).UseParams  = opts.UseParams{iTop};
    optsTop(iTop).UseICs     = opts.UseICs{iTop};
    if ~opts.UseModelICs
        optsTop(iTop).UseICs = optsTop(iTop).UseICs(:,1:nCon);
    end
    optsTop(iTop).AbsTol     = opts.AbsTol(1:nCon,iTop);
    optsTop(iTop).LowerBound = opts.LowerBound{iTop};
    optsTop(iTop).UpperBound = opts.UpperBound{iTop};
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).Name ' to objectives...\n']); end
        m(iTop) = FitObjective(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    end
end

%% Active parameters
T = cell(nTop,1);
for iTop = 1:nTop
    T{iTop} = collectActiveParameters(m(iTon), con(iCon,iTop), opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseModelICs);
end

%% Compute p_y|m for each topology
pym = zeros(nTop,1);
pyTm = ones(nTop,1);
F = cell(nTop,1);
for iTop = 1:nTop
    if verbose; fprintf(['Computing probability of ' m(iTop).Name '...\n']); end
    % Likelihood
    pyTm(iTop) = ObjectiveProbability(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    
    % Information
    F{iTop} = ObjectiveInformation(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    F{iTop} = posdef(F{iTop});
    
    % Equation for linearization about maximum a posteriori
    pym(iTop) = pyTm(iTop) * (2*pi).^(nT(iTop)/2) * det(F{iTop}).^(-1/2);
end

%% Compute p_m|y for the set
% Apply topology prior
pmy = pym .* opts.PriorTopology;

% Normalize 
pmy = pmy ./ sum(pmy);
