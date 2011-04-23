function [best data] = BestTopologyExperiment(m, con, obj, posCon, posObj, target, opts, pmyStart, pmyNew)

% Clean-up inputs
if nargin < 9
    pmyNew = [];
    if nargin < 8
        pmyStart = [];
        if nargin < 7
            opts = [];
        end
    end
end

% Constants
nTop = numel(m);
nCon = size(con,1);
nObj = size(obj,2);
nPosCon = size(posCon,1);
nPosObj = size(posObj,2);

% Options
defaultOpts.RelTol      = NaN; % 1e-6
defaultOpts.AbsTol      = NaN; % in fixAbsTol
defaultOpts.UseModelICs = true;

defaultOpts.UseAdjoint  = false;
defaultOpts.UseParams   = NaN; % All kinetic parameters
defaultOpts.UseICs      = NaN; % No initial conditions

defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.NeedFit        = true; % Fit the parameters

defaultOpts.BestTopologyMethod = 'MonteCarlo';
defaultOpts.TargetTol = 0.05;
defaultOpts.MinTargetIter = 0;
defaultOpts.MaxTargetIter = 200;

defaultOpts.Verbose     = 0;

opts = mergestruct(defaultOpts, opts);
verbose = logical(max(opts.Verbose, 0));
verbose2 = logical(max(opts.Verbose-1, 0));
opts.Verbose = max(opts.Verbose-2, 0);

% Constants
nx = zeros(nTop,1);
nk = zeros(nTop,1);
for iTop = 1:nTop
    nx(iTop) = m(iTop).nX;
    nk(iTop) = m(iTop).nP;
end

%% Standardize structures
% Experimental conditions
assert((size(con,2) == 1 && (opts.UseModelICs || all(nx == nx(1)))) || size(con,2) == nTop, 'KroneckerBio:TopologyProbability:ConSize', 'Second dimension of input "con" must be equal to numel(m) or 1 if opts.UseModelICs is false or every model has the same number of species')
if size(con,2) == 1
    con = repmat(con, 1,nTop);
end
con = refreshCon(m, con);

assert((size(posCon,2) == 1 && (opts.UseModelICs || all(nx == nx(1)))) || size(posCon,2) == nTop, 'KroneckerBio:TopologyProbability:PosConSize', 'Second dimension of input "posCon" must be equal to numel(m) or 1 if opts.UseModelICs is false or every model has the same number of species')
if size(posCon,2) == 1
    posCon = repmat(posCon, 1,nTop);
end
posCon = refreshCon(m, posCon);

% Objective functions
assert(size(obj,3) == 1 || size(obj,3) == nTop, 'KroneckerBio:TopologyProbability:ObjSize', 'Third dimension of input "obj" must be equal to numel(m) or 1')
if size(obj,3) == 1
    obj = repmat(obj, [1,1,nTop]);
end
obj = refreshObj(m, obj);

% Possible objective functions must be blind to the topology
assert(size(posObj,3) == 1, 'KroneckerBio:TopologyProbability:PosObjSize', 'Third dimension of input "posObj" must be equal to 1')
posObj = refreshObj(m(1), posObj);

%% Active Parameters
% Default UseParams is all kinetic parameters
if isempty(opts.UseParams) || any(isnan(opts.UseParams))
    opts.UseParams = cell(nTop,1);
    for iTop = 1:nTop
        opts.UseParams{iTop} = 1:m(iTop).nP;
    end
end

% Default UseICs is no initial condition parameters
if isempty(opts.UseICs) || any(isnan(opts.UseICs))
    opts.UseICs = cell(nTop,1);
    for iTop = 1:nTop
        opts.UseICs{iTop} = [];
    end
end

% Ensure UseParams is column vector of logical indexes
nTk = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseParams{iTop} nTk(iTop)] = fixUseParams(opts.UseParams{iTop}, nk(iTop));
end

% Ensure UseICs is a matrix of logical indexes
nTx = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseICs{iTop} nTx(iTop)] = fixUseICs(opts.UseICs{iTop}, opts.UseModelICs, nx(iTop), nCon + nPosCon);
end

nT = nTk + nTx;

%% Integration type: simple, continuous, complex, or both
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);
[posContinuous, posComplex, postGet] = fixIntegrationType(posObj);

%% Tolerances
if isnan(opts.RelTol) || isempty(opts.RelTol) 
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
if isnumeric(opts.AbsTol)
    opts.AbsTol = repmat(vec(opts.AbsTol), 1,nTop);
end
tempAbsTol = cell(nCon+nPosCon,nTop);
for iTop = 1:nTop
    tempAbsTol(:,iTop) = fixAbsTol(opts.AbsTol(:,iTop), 2, false(nCon+nPosCon,1), nx(iTop), nCon + nPosCon, opts.UseAdjoint, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseModelICs);
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
    optsTop(iTop).continuous = opts.continuous(:,iTop);
    optsTop(iTop).complex    = opts.complex(:,iTop);
    optsTop(iTop).tGet       = opts.tGet(:,iTop);
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).name ' to existing data...\n']); end
        m(iTop) = FitObjective(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    end
    opts.NeedFit = false;
end

%% Information and starting p_m|y
% Information
V = cell(nTop,1);
for iTop = 1:nTop
    if verbose; fprintf(['Computing posterior variance of ' m(iTop).name '...\n']); end
    V{iTop} = ObjectiveInformation(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    V{iTop} = posdef(V{iTop}); % Posterior information
    V{iTop} = inv(V{iTop});    % Posterior variance
end

% Starting p_m|y
if isempty(pmyStart)
    if verbose; fprintf('Computing starting probability...\n'); end
    pmyStart = TopologyProbability(m, con, obj, opts);
end

%% Monte Carlo parameters
nCarlo = opts.MaxTargetIter; % Maximum number of Monte Carlo iterations 
pdrawmax = zeros(nTop,1); % Maximum value coming from mvnpdf on Monte Carlo parameters
prealmax = zeros(nTop,1); % Maximum value coming from likelihood for fit to data

for iTop = 1:nTop
    pdrawmax(iTop) = mvnpdf(zeros(nT(iTop),1), zeros(nT(iTop),1), V{iTop});
    prealmax(iTop) = ObjectiveProbability(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
end

%% Expected target function value for each possible experiment
Etarget = zeros(nPosCon,nPosObj);
pmyAll    = cell(nPosCon,nPosObj,nTop);
targetAll = cell(nPosCon,nPosObj,nTop);
weightAll = cell(nPosCon,nPosObj,nTop);
wmeanAll  = cell(nPosCon,nPosObj);
wvarAll   = cell(nPosCon,nPosObj);
errAll    = cell(nPosCon,nPosObj);
for iPosCon = 1:nPosCon
    for iPosObj = 1:nPosObj
        wmeanAll{iPosCon,iPosObj}  = zeros(nCarlo,1);
        wvarAll{iPosCon,iPosObj}   = zeros(nCarlo,1);
        errAll{iPosCon,iPosObj}    = zeros(nCarlo,1);
        for iTop = 1:nTop
            pmyAll{iPosCon,iPosObj,iTop} = zeros(nTop,nCarlo);
        end
    end
end

% Initialize Gibbs sampler
iSam = inf(nTop,1); % Start with replenishment
sam = cell(nTop,1);
mu  = cell(nTop,1);
lb  = cell(nTop,1);
ub  = cell(nTop,1);
for iTop = 1:nTop
    sam{iTop} = collectActiveParameters(m(iTop), con(:,iTop), optsTop(iTop).UseParams, optsTop(iTop).UseICs, opts.UseModelICs);
    mu{iTop}  = m(iTop).p(optsTop(iTop).UseParams);
    lb{iTop}  = optsTop(iTop).LowerBound;
    ub{iTop}  = optsTop(iTop).UpperBound;
    if opts.Normalized
        sam{iTop} = log(sam{iTop});
        mu{iTop}  = log(mu{iTop});
        lb{iTop}  = log(lb{iTop});
        ub{iTop}  = log(ub{iTop});
    end
end

for iPosCon = 1:nPosCon
    for iPosObj = 1:nPosObj
        if verbose; fprintf(['Monte Carlo sampling for ' posObj(iPosCon,iPosObj,1).name '...\n']); end
        
        index = 0;
        indexes = zeros(nTop,1);
        
        % Reset growing variables
        targets = cell(iTop,1);
        weights = cell(iTop,1);
        for iTop = 1:nTop
            targets{iTop} = zeros(nCarlo,1);
            weights{iTop} = zeros(nCarlo,1);
        end
        
        while true %dowhile
            % Draw random topology
            if verbose; fprintf('Drawing random model...\n'); end
            trueTop = randp(pmyStart);
            index = index + 1; % Total count for this experiment
            indexes(trueTop) = indexes(trueTop) + 1; % Count for a particular topology
            if verbose2; fprintf('Model %d (%s) was chosen\n', trueTop, m(trueTop).name); end
            
            % Draw parameter set
            while true %dowhile
                % Replenish sample supply when empty
                if iSam(trueTop) > size(sam{trueTop}, 2)
                    sam{trueTop} = mvnbndrndgibbs(mu{trueTop}, V{trueTop}, lb{trueTop}, ub{trueTop}, sam{trueTop}(:,end), nCarlo, 0, 0);
                    sam{trueTop} = sam{trueTop}(:,randperm(nCarlo)); % Scramble to reduce serial correlation
                    iSam(trueTop) = 1;
                end
                
                % Fetch random parameter set to try
                Trand = sam{trueTop}(:,iSam(trueTop));
                chi2Trand = (mu{trueTop} - Trand).' * (V{trueTop} \ (mu{trueTop} - Trand)); % Needed to test if this draw is safe
                if opts.Normalized
                    Trand = exp(Trand);
                end
                iSam(trueTop) = iSam(trueTop) + 1;
                
                % Probability of drawing this set
                if opts.Normalized
                    pdraw = mvnpdf(log(Trand), log(m(trueTop).p(optsTop(trueTop).UseParams)), V{trueTop}) / pdrawmax(trueTop); % Chance of drawing in monte carlo
                else
                    pdraw = mvnpdf(Trand, m(trueTop).p(optsTop(trueTop).UseParams), V{trueTop}) / pdrawmax(trueTop); % Chance of drawing in monte carlo
                end
                if verbose2; fprintf('pdraw = %-6.4g\n', pdraw); end
                
                if chi2pvalue(chi2Trand, nT(trueTop)) > opts.TargetTol % pdraw should be safe
                    % Update with new parameter set
                    [mRand conRand objRand] = updateAll(m(trueTop), con(:,trueTop), obj(:,:,trueTop), Trand, optsTop(trueTop).UseParams, optsTop(trueTop).UseICs, opts.UseModelICs);
                    
                    % Probability that this set is real
                    preal = ObjectiveProbability(mRand, conRand, objRand, optsTop(trueTop)) / prealmax(trueTop);
                    if verbose2; fprintf('preal = %-6.4g\n', preal); end
                    
                    % Weight of this draw
                    weight = preal / pdraw;
                    
                    if weight > opts.TargetTol % preal should be meaningful
                        % Contributable parameter set found
                        break
                    end
                end
            end % Gibbs
            
            % Generate data
            if verbose; fprintf('Generating data from Monte Carlo model...\n'); end
            optsSub = opts;
            optsSub.AbsTol = opts.AbsTol{nCon+iPosCon,trueTop};
            optsSub.tGet   = postGet{iPosCon};
            if posContinuous(iPosCon) && posComplex(iPosCon)
                sol = integrateObj(mRand, posCon(iPosCon), posObj(iPosCon,iPosObj), optsSub);
            elseif posComplex(iPosCon)
                sol = integrateSys(mRand, posCon(iPosCon,trueTop), optsSub);
            elseif posContinuous(iPosCon)
                sol = integrateObjDisc(mRand, posCon(iPosCon,trueTop), posObj(iPosCon,iPosObj), optsSub);
            else
                sol = integrateSysDisc(mRand, posCon(iPosCon,trueTop), optsSub);
            end
            
            % Create objective function appropriate to each topology
            sol.c = mRand.c;
            newObj = pastestruct(Gzero(mRand), posObj(iPosCon,iPosObj).AddData(sol, posCon(iPosCon,iTop).u));
            
            % Fit the new data
            mfit = m;
            for iTop = 1:nTop
                if verbose; fprintf(['Fitting ' m(iTop).name ' to existing data plus Monte Carlo data...\n']); end
                optsSub = optsTop(iTop);
                if ~opts.UseModelICs
                    optsSub.UseICs = [optsSub.UseICs, opts.UseICs{iTop}(:,nCon+iPosCon)];
                end
                optsSub.AbsTol = [optsSub.AbsTol; opts.AbsTol(nCon+iPosCon,iTop)]; % Final experiment has no tolerance needs
                mfit(iTop) = FitObjective(mfit(iTop), [con(:,iTop); posCon(iPosCon,iTop)], [obj(:,:,iTop); [newObj, repmat(Gzero(m(iTop)), 1,nObj-1)]], optsSub);
            end
            
            % Predicted p_m|y
            if verbose; fprintf('Computing topology probability after adding Monte Carlo data...\n'); end
            optsSub = opts;
            optsSub.AbsTol = optsSub.AbsTol([1:nCon,iPosCon],:);
            if ~opts.UseModelICs
                for iTop = 1:nTop
                    optsSub.UseICs{iTop} = optsSub.UseICs{iTop}(:,[1:nCon,iPosCon]);
                end
            end
            pmy = TopologyProbability(mfit, [con; posCon(iPosCon,:)], [obj; [repmat(newObj, [1,1,nTop]), repmat(Gzero(m(1)), [1,nObj-1,nTop])]], optsSub);
            
            % Evaluate target function
            targets{trueTop}(indexes(trueTop)) = target(pmy);
            if verbose; fprintf('Target Value = %-6.4g\n', targets{trueTop}(indexes(trueTop))); end
            
            % Importance sampling weight
            weights{trueTop}(indexes(trueTop)) = weight; % Keep all the weights
            
            % Reweight the importance weights to balance the topologies
            targetStacked     = zeros(index,1); % Stack targets into a vector
            reweightedWeights = zeros(index,1);
            rwwInd = 0;
            for iTop = 1:nTop
                % Compute mean topology weight
                topoWeight = mean(weights{iTop}(1:indexes(iTop)));
                if isnan(topoWeight)
                    % This topology has no draws yet, ignore it
                    continue
                end
                
                % Reweight accordingly
                nextrwwInd = rwwInd+indexes(iTop);
                targetStacked(rwwInd+1:nextrwwInd) = targets{iTop}(1:indexes(iTop));
                reweightedWeights(rwwInd+1:nextrwwInd) = weights{iTop}(1:indexes(iTop)) ./ topoWeight;
                rwwInd = nextrwwInd;
            end
            
            wmean  = weightedmean(targetStacked, reweightedWeights, 1);
            wvar   = weightedvar(targetStacked, reweightedWeights, 1);
            err    = sqrt(wvar * sum((reweightedWeights./sum(reweightedWeights)).^2));
            if verbose; fprintf('Mean = %-6.4g Error = %-6.4g\n', wmean, err); end
            
            % Store iteration results
            pmyAll{iPosCon,iPosObj,trueTop}(:,indexes(trueTop)) = pmy;
            wmeanAll{iPosCon,iPosObj}(index) = wmean;
            wvarAll{iPosCon,iPosObj}(index)  = wvar;
            errAll{iPosCon,iPosObj}(index)   = err;
            
            if ((err <= opts.TargetTol) && (index >= opts.MinTargetIter)) || (index >= opts.MaxTargetIter)
                if verbose; fprintf('Monte Carlo terminated\n'); end
                Etarget(iPosCon,iPosObj) = wmeanAll{iPosCon,iPosObj}(index);
                for iTop = 1:nTop
                    pmyAll{iPosCon,iPosObj,iTop} = pmyAll{iPosCon,iPosObj,iTop}(:,1:indexes(iTop));
                    targetAll{iPosCon,iPosObj,iTop} = targets{iTop}(1:indexes(iTop));
                    weightAll{iPosCon,iPosObj,iTop} = weights{iTop}(1:indexes(iTop));
                end
                wmeanAll{iPosCon,iPosObj}  = wmeanAll{iPosCon,iPosObj}(1:index);
                wvarAll{iPosCon,iPosObj}   = wvarAll{iPosCon,iPosObj}(1:index);
                errAll{iPosCon,iPosObj}    = errAll{iPosCon,iPosObj}(1:index);
                break
            end
        end % Monte Carlo
    end % iPosObj
end % iPosCon

[unused best] = min(Etarget); %#ok<ASGLU>

if nargout >= 2
    data.Targets               = Etarget;
    data.Startingpmy           = pmyStart;
    data.Allpmy                = pmyAll;
    data.AllTargets            = targetAll;
    data.AllWeights            = weightAll;
    data.AllWeightedMeans      = wmeanAll;
    data.AllWeightedVars       = wvarAll;
    data.AllWeightedMeanErrors = errAll;
end
