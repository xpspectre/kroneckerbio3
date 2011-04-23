function [GTolRatio DTolRatio] = OptimalAbsTol(m, con, obj, opts)

%% Options
defaultOpts.Coverage                = 1;
defaultOpts.UseParams               = 1:m.nP;
defaultOpts.UseICs                  = [];
defaultOpts.UseModelICs             = true;
defaultOpts.ObjWeights              = ones(size(obj));
defaultOpts.UseAdjoint              = true;
defaultOpts.Verbose                 = 1;

defaultOpts.AbsTol                  = NaN; % in fixAbsTol
defaultOpts.RelTol                  = NaN; % 1e-6

opts = mergestruct(defaultOpts, opts);

% Constants
nX = m.nX;
nP = m.nP;
nCon = numel(con);
nObj = size(obj,2);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nVP] = fixUseParams(opts.UseParams, nP);

% Ensure UseICs is a matrix of logical indexes
[opts.UseICs, nVX] = fixUseICs(opts.UseICs, opts.UseModelICs, nX, nCon);

nV = nVX + nVP;

%% Integration type: simple, continuous, complex, or both
% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);

%% Tolerances
if isempty(opts.RelTol) || isnan(opts.RelTol)
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nX, nCon, opts.UseAdjoint, opts.UseParams, opts.UseICs, opts.UseModelICs, true);

%% Determine the order of the integration
if nargout <= 1
    order = 0;
elseif nargout <= 2
    order = 1;
end

intOpts = opts;
intOpts.Verbose = max(opts.Verbose-1,0);

GTolRatio = cell(nCon,1);
DTolRatio = cell(nCon,1);

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.tGet = opts.tGet{iCon};
    intOpts.ObjWeights = opts.ObjWeights(iCon,:);
    
    if opts.Verbose; fprintf(['Computing optimal tolerance ratio for ' con(iCon).name '...\n']); end
    % Integrate the basic system
    if order == 0
        if opts.continuous(iCon)
            sol = integrateObjFwd(m, con, obj, intOpts);
        else
            sol = integrateObjFwdComplex(m, con, obj, intOpts);
        end
    elseif order == 1
        if opts.continuous(iCon)
            sol = integrateObjSensFwd(m, con, obj, intOpts);
        else
            [GTolRatio{iCon} DTolRatio{iCon}] = integrateOptimalAbsTolSensSimple(m, con, obj, intOpts);
        end
    end
    
end
