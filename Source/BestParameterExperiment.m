function [bestCons data] = BestParameterExperiment(m, con, obj, posCon, posObj, goal, opts, H, EHs)
%BESTPARAMETEREXPERIMENT Determine which experiments will most efficiently
%   minimize a goal function of the fisher information matrix (FIM) of the
%   objective function. Traditionally, this algorithm is used to minimize
%   some uncertainty in the parameters directions. (The variance-covariance
%   matrix of the parameters is approximated by the inverse of the
%   FIM if the objective function is some form of least squares.) 
%
%   This algorithm computes the current FIM according to the current
%   experiments defined by con and obj. It then adds the expected FIM from
%   each of the possible experiments defined by posCon and posObj and asks
%   the goal function to evaluate the resulting FIM. Finally, the function
%   returns the index of the experiment that has the lowest expected goal
%   function value.
%
%   [...] = BestParameterExperiment(m, con, obj, posCon, posObj, goal, opts)
%
%   Inputs:
%       m   - The KroneckerBio model of interest
%       con - A vector of structures, with each item corresponding to the
%             experimental conditions of the experiments that have already
%             been performed
%       obj - A vector of structures, with each item corresponding to the
%             objective function containing the measurements for the
%             experiment
%       posCon - A vector of experimental conditions that will be examined
%       posObj - A vector of objective functions corresponding to the
%                measurement techique that will be applied under the
%                possible experimental conditions
%       goal   - A function handle @(H) that returns a scalar evaluation of
%                a hessian. This function chooses the best experiments
%                based on which set has the lowest goal function.
%       opts - Optional function options
%           UseParams - Vector of indexes indicating the rate constants
%                       that will be considered in the hessian.
%           UseICs    - Vector of indexes indicating the initial
%                       concentrations whose sensitivities will be
%                       considered
%           UseModelICs   - Boolean that determines whether to use the
%                           initial conditions of the model or the
%                           conditions. This is used both for simulating
%                           the model and calculating the sensitivities of
%                           the IC parameters.
%           ReturnCount   - Integer indicating the number of experiments
%                           to return
%           MaxGreedySize - Integer indicating how many experiments should
%                           be considered simultaneously. If it set to 1,
%                           the algorithm is perfectly greedy, choosing
%                           the best first experiment and then choosing
%                           the second best, given that the first was
%                           already chosen. If it is greater than or equal
%                           to ReturnCount, then the algorithm is
%                           perfectly optimal, choosing the combination
%                           that all together minimizes the goal function.
%           TerminalGoal  - Scalar that terminates the search if the
%                           expected goal is less than or equal to this
%                           value. Default = -inf
%           AllowRepeats  - Boolean that determines whether or to allow
%                           experiments to be repeated
%           Verbose       - Print progress to command window
%
%   Outputs:
%       bestCons = FindBestExperiment(m, con, obj, posCon, posObj, goal, opts)
%       	bestCons - A vector of the indexes to the best experiments
%
%       [bestCons data] = FindBestExperiment(m, con, obj, posCon, posObj, goal, opts)
%        	data - Includes a structure that holds the hessians and
%                  goal function values at each iteration. This data
%                  may be useful for analysis. It also can require a lot of
%                  memory.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargin >= 6, 'KroneckerBio:BestParameterExperiment:AtLeastSixInputs', 'FindBestExperiment requies at least 6 input arguments.')
if nargin < 9
    EHs = [];
    if nargin < 8
        H = [];
        if nargin < 7
            opts = [];
        end
    end
end

% Options
defaultOpts.UseParams           = 1:m.nP;
defaultOpts.UseICs              = [];
defaultOpts.UseModelICs         = true;
defaultOpts.ReturnCount         = 1;        % Number of experiments to return
defaultOpts.MaxGreedySize       = 1;        % Number of experiments to consider at once for the greedy search. Inf = not greedy
defaultOpts.TerminalGoal        = -inf;
defaultOpts.AllowRepeats        = true;
defaultOpts.Verbose             = false;

opts = mergestruct(defaultOpts, opts);

% Constants
nPK = length(opts.UseParams);
nPX = length(opts.UseICs);
nP = nPK + nPX;
nCon = length(con);
nPos = numel(posObj);
blockSize = min(opts.MaxGreedySize, opts.ReturnCount);  % The number of experiments over which to maximize cannot be greater than the number of experiments to be returned
maxSearchSize = nchoosek(nPos, blockSize);
nIterations = ceil(opts.ReturnCount / blockSize);

%% Compute old hessian
if isempty(H)
    H = ObjectiveInformation(m, con, obj, opts);
end

%% Fetch expected hessians
if isempty(EHs)
    [unused EHs] = ObjectiveInformation(m, posCon, posObj, opts);
end

%% Find best experiment
if opts.Verbose; fprintf('Combining information to find best experiment...\n'); end
    
bestCons = zeros(opts.ReturnCount, 1);
bestGoal = goal(H);

if opts.Verbose
    fprintf('Starting goal is %d\n', bestGoal)
end

% Allocate algorithm storage if it is requested
if nargout >= 2
    allHessians = cell(nIterations, maxSearchSize);
    allGoals = zeros(nIterations, maxSearchSize);
    bestHessians = cell(nIterations, 1);
    bestGoals = zeros(nIterations, 1);
    currentIteration = 1; % Keeps track of which iteration the loop is on
end

remainingCount = opts.ReturnCount; % The number of best experiments remaining to be found
remainingCons = (1:numel(EHs)).'; % The indexes of the conditions allowed to be chosen
currentH = H; % Updates after each block to hold the current expected value of the hessian

while (remainingCount > 0 && bestGoal > opts.TerminalGoal)
    if opts.Verbose; fprintf('%d experiments remaining\n', remainingCount); end
    
    % CurrentBlockSize is the set size of the experiments to optimize over.
    % It is equal to blockSize until the number of experiments remaining to
    % be found is less than blockSize.
    currentBlockSize = min(blockSize, remainingCount);
    
    % Combinatorics
    rCons = length(remainingCons);
    conList = combinator(rCons, currentBlockSize, 'c', 'r');   % Generate list of combinations of experiments
    nSearch = size(conList, 1);                               % Total number of possible sets
    
    if nargout >= 2
        % Compute expected hessians
        newEHs = cell(nSearch,1);
        for iSearch = 1:nSearch
            newEHs{iSearch} = currentH;                                                            % Expected hessian is old hessian...
            for iBlock = 1:currentBlockSize
                newEHs{iSearch} = newEHs{iSearch} + EHs{remainingCons(conList(iSearch, iBlock))};  % ...plus contributions from each hypothetical experiment
            end
        end
        
        % Compute goal values
        goalValues = zeros(nSearch, 1);
        for iSearch = 1:nSearch
            goalValues(iSearch) = goal(newEHs{iSearch});
        end
        
        % Pick best goal
        bestGoal = min(goalValues);
        bestSetInd = find(bestGoal == goalValues, 1);                       % Row of conList corresponding to the minimum goal function
    else
        % This way takes less memory by not storing data
        bestGoal = inf;
        bestSetInd = 1;
        for iSearch = 1:nSearch
            % Compute expected hessian
            newEH = currentH;
            for iBlock = 1:currentBlockSize
                newEH = newEH + EHs{remainingCons(conList(iSearch, iBlock))};
            end
            
            % Compute goal and keep it if it is better
            goalValue = goal(newEH);
            if goalValue < bestGoal
                bestGoal = goalValue;
                bestSetInd = iSearch;
            end
        end
    end
    
    % Find corresponding best set
    bestSet = conList(bestSetInd, :);                                   % Experiments that are part of the best set (indexes of remainingCons)
    bestTrueSet = remainingCons(bestSet);                               % Those experiments according to the indexes of obj
    blockInd1 = opts.ReturnCount - remainingCount + 1;                  % Where in bestCons to starting storing these experiment indexes
    blockInd2 = opts.ReturnCount - remainingCount + currentBlockSize;   % Where in bestCons is the last index to store
    bestCons(blockInd1:blockInd2) = bestTrueSet;                        % Store values of the best set
    
    if opts.Verbose
        for i = 1:currentBlockSize
            fprintf([posCon(bestTrueSet(i)).name ' was chosen\n']);
        end
        fprintf('Expected goal is %d\n', bestGoal);
    end
    
    % Apply information for best set
    for i = 1:currentBlockSize
        currentH = currentH + EHs{bestTrueSet(i)};
    end
    
    % Remove conditions if repeats are not allowed
    if ~opts.AllowRepeats
        remainingCons(bestSet) = [];
    end
    
    % Store iteration data only if it is requested
    if nargout >= 2
        allHessians(currentIteration, 1:nSearch) = newEHs;
        allGoals(currentIteration, 1:nSearch) = goalValues;
        bestHessians{currentIteration} = currentH;
        bestGoals(currentIteration) = bestGoal;
        currentIteration = currentIteration + 1;
    end
    
    % Decrement
    remainingCount = remainingCount - currentBlockSize;
    
end

%% Work-down
% Clean up if terminated on TerminalGoal rather than ReturnCount
bestCons = bestCons(bestCons ~= 0);

if nargout >= 2
    data.StartingFIM = H;            % FIM of system before possible experiments are applied
    data.AllFIMs     = EHs;          % FIMs for each possible experiment
    data.SearchFIMs  = allHessians;  % Sum of FIMs in each search
    data.AllGoals    = allGoals;     % All the goal values of the possible iterations
    data.BestFIMs    = bestHessians; % FIM after each iteration
    data.BestGoals   = bestGoals;    % Goal after each iteration
end

end
