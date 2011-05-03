function [GTolRatio DTolRatio] = OptimalAbsTol(m, con, obj, opts)

%% Options
defaultOpts.Coverage                = 1;

defaultOpts.UseModelICs             = true;
defaultOpts.UseModelInputs          = true;
defaultOpts.UseParams               = 1:m.nk;
defaultOpts.UseICs                  = [];
defaultOpts.UseControls             = [];%TODO

defaultOpts.ObjWeights              = ones(size(obj));
defaultOpts.Normalized              = true;
defaultOpts.UseAdjoint              = false;

defaultOpts.AbsTol                  = NaN; % in fixAbsTol
defaultOpts.RelTol                  = NaN; % in fixRelTol
defaultOpts.Verbose                 = 1;

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

nT = nTk + nTx;

%% Integration type: simple, continuous, complex, or both
% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, opts.UseAdjoint, opts.UseParams, opts.UseICs, opts.UseModelICs, true);

%% Determine the order of the integration
if nargout <= 1
    order = 0;
elseif nargout <= 2
    order = 1;
end

intOpts = opts;

GTolRatio = cell(nCon,1);
DTolRatio = cell(nCon,1);

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.tGet = opts.tGet{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    
    if verbose; fprintf(['Computing optimal tolerance ratio for ' con(iCon).Name '...\n']); end
    % Integrate the basic system
    if order == 0
        if opts.continuous(iCon)
            %sol = integrateObjFwd(m, con, obj, intOpts);
        else
            %sol = integrateObjFwdComplex(m, con, obj, intOpts);
        end
    elseif order == 1
        if opts.continuous(iCon)
            %sol = integrateObjSensFwd(m, con, obj, intOpts);
        else
            [GTolRatio{iCon} DTolRatio{iCon}] = integrateOptimalAbsTolSensSimple(m, con, obj, intOpts);
        end
    end
    
end
