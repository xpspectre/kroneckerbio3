function useExperiments = ParameterExperimentElimination(m, con, obj, conPos, objPos, opts, F, EFs)
% Clean up inputs
assert(nargin >= 6, 'KroneckerBio:BestParameterExperiment:AtLeastSixInputs', 'FindBestExperiment requies at least 6 input arguments.')
if nargin < 8
    EFs = [];
    if nargin < 7
        F = [];
        if nargin < 6
            opts = [];
        end
    end
end

% Options
defaultOpts.UseExperiments      = true(size(obj));
defaultOpts.Cost                = zeros(size(obj));
defaultOpts.Budget              = inf;
defaultOpts.Verbose             = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nPosCon = numel(conPos);
nPosObj = size(objPos,1);

%% Fix UseExperiments
useExperiments = fixUseExperiments(opts.UseExperiments, nPosObj, nPosCon);
remainingCons = find(opts.UseExperiments); % The indexes of the conditions allowed to be chosen
nPos = numel(remainingCons);

%% Compute existing information
if isempty(F)
    F = ObjectiveInformation(m, con, obj, opts);
end
nT = size(F,1);

%% Fetch expected hessians
if isempty(EFs)
    [unused EFs] = ObjectiveInformation(m, conPos, objPos, opts);
end

%% Generate all information singles
%EFsingles = cellfun(@plus, repmat({F}, size(EFs)), EFs, 'UniformOutput', false);
%EFsingles = cellfun(@plus, repmat({eye(nT)}, size(EFs)), EFs, 'UniformOutput', false);
EFsingles = EFs;

%% Single elimination
rCons = numel(remainingCons);
conList = remainingCons(combinator(rCons, 2, 'c', 'n'));   % Generate list of pairs of experiments, no repeats
for iPos = 1:size(conList,1)
    firstCon = conList(iPos,1);
    secondCon = conList(iPos,2);
    [iObj1 iCon1] = ind2sub([nPosObj nPosCon], firstCon);
    [iObj2 iCon2] = ind2sub([nPosObj nPosCon], secondCon);
    
    firstZeros  = all(EFsingles{firstCon} == 0);
    secondZeros = all(EFsingles{secondCon} == 0);
    firstInf    = any(EFsingles{firstCon} == inf);
    secondInf   = any(EFsingles{secondCon} == inf);
    
    % Domination: 0 tied, 1 first dominates, 2 second dominates, -1 complimentary
    if all(firstZeros == secondZeros)
        zeroDom = 0;
    elseif all(firstZeros >= secondZeros)
        zeroDom = 2;
    elseif all(firstZeros <= secondZeros)
        zeroDom = 1;
    else
        zeroDom = -1;
    end
    if all(firstInf == secondInf)
        infDom = 0;
    elseif all(firstInf >= secondInf)
        infDom = 1;
    elseif all(firstInf <= secondInf)
        infDom = 2;
    else
        infDom = -1;
    end
    
    % Only continue if they are not already mutually excllusive
%     if ((zeroDom == 1 || zeroDom == 0) && (infDom == 1 || infDom == 0) && (opts.Cost(firstCon) <= opts.Cost(secondCon))) || ...
%             (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
%         % Compute dominance in sensible parameter space
%         % Invert
%         Vfirst  = infoinv(EFsingles{firstCon});
%         Vsecond = infoinv(EFsingles{secondCon});
%         
%         % Remove directions already under consideration
%         keepDir = ~(firstZeros | secondZeros | firstInf | secondInf);
%         keepSize = sum(keepDir);
%         Vfirst  = Vfirst(keepDir,keepDir);
%         Vsecond = Vsecond(keepDir,keepDir);
%         
%         % Eigendecompose
%         [sigmafirst Qfirst]   = infoeig(Vfirst);
%         [sigmasecond Qsecond] = infoeig(Vsecond);
%         
%         % Rotate second to basis of first
%         Ftransformed = (spdiags(sqrt(sigmasecond), 0, keepSize, keepSize) * Qsecond.') * (Qfirst * spdiags(sqrt(1./sigmafirst), 0, keepSize, keepSize));
%         Ftransformed = Ftransformed.' * Ftransformed;
%         
%         % Extract eigenvalues
%         sigmatransformed = infoeig(Ftransformed);
%         
%         % Compare for dominance
%         if all(sigmatransformed >= 1 - 1e-6) && (zeroDom == 1 || zeroDom == 0) && (infDom == 1 || infDom == 0) && (opts.Cost(firstCon) <= opts.Cost(secondCon))
%             if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', conPos(iCon1).Name, objPos(firstCon).Name, firstCon, conPos(iCon2).Name, objPos(secondCon).Name, secondCon); end
%             useExperiments(secondCon) = false;
%         elseif all(sigmatransformed <= 1 + 1e-6) && (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
%             if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', conPos(iCon2).Name, objPos(secondCon).Name, secondCon, conPos(iCon1).Name, objPos(firstCon).Name, firstCon); end
%             useExperiments(firstCon) = false;
%         end
%     end

    % Only continue if they are not already mutually exclusive
    if ((zeroDom == 1 || zeroDom == 0) && (infDom == 1 || infDom == 0) && (opts.Cost(firstCon) <= opts.Cost(secondCon))) || ...
            (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
        % Find difference in information between experiments
        Fdiff = EFsingles{firstCon} - EFsingles{secondCon};
        
        % Eigendecompose the difference
        lambdadiff = infoeig(Fdiff);
        
        % Compare for dominace
        if all(lambdadiff >= -sqrt(opts.RelTol*32)) && (zeroDom == 1 || zeroDom == 0) && (infDom == 1 || infDom == 0) && (opts.Cost(firstCon) <= opts.Cost(secondCon))
            if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', conPos(iCon1).Name, objPos(firstCon).Name, firstCon, conPos(iCon2).Name, objPos(secondCon).Name, secondCon); end
            useExperiments(secondCon) = false;
        elseif all(lambdadiff <= sqrt(opts.RelTol*32)) && (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
            if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', conPos(iCon2).Name, objPos(secondCon).Name, secondCon, conPos(iCon1).Name, objPos(firstCon).Name, firstCon); end
            useExperiments(firstCon) = false;
        end
    end

end