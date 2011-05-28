function p = ObjectiveProbability(m, con, obj, opts)
% p = ObjectiveProbability(m, con, obj, opts)

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:ObjectiveProbability:MoreThanOneModel', 'The model structure must be scalar')

% Options
defaultOpts.UseModelICs    = true;
defaultOpts.UseModelInputs = false;
defaultOpts.ObjWeights     = ones(size(obj));
defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.Verbose        = 0;

opts = mergestruct(defaultOpts, opts);
verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nCon = numel(con);
nObj = size(obj,1);

% Add missing fields to structure
con = pastestruct(Uzero(m), con);
obj = pastestruct(Gzero(m), obj);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

%% Tolerances
% RelTol
if isempty(opts.RelTol) || isnan(opts.RelTol)
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon);

%% Probability
% Initialize variables
p = 1;
intOpts = opts;

for iCon = 1:nCon
    if verbose; disp(['Integrating system for ' con(iCon).name '...']); end
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

