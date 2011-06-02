function [F All] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%OBJECTIVEINFORMATION Compute the fisher information matrix of an objective
%   function.
%
%   Given a model and experimental conditions, this function computes the
%   fisher information matrix provided by the objective function, given the
%   model and conditions and evaluated at p.
%
%   Mathematically: F = dy/dp' * V^-1 * dy/dp
%
%   [...] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A structure vector of the experimental conditions under
%              which the hessian will be evaluated
%       obj  - A structure array of the objective functions under which the
%              hessian will be evaluated. Each row of obj is matched to the
%              corresponding entry in con.
%       opts - Optional function options
%           .UseParams   - Vector of indexes indicating the rate constants
%                          whose sensitivities will be considered
%           .UseICs      - Vector of indexes indicating the initial
%                          concentrations whose sensitivities will be
%                          considered
%           .UseModelICs - Boolean that determines whether to use the
%                          initial concentrations of the model or the
%                          conditions. Default = true
%           .Normalized  - Logical that determines if the simple
%                          information or the normalized information will
%                          be computed. The normalized information is
%                          normalized with respect to the values of the
%                          parameters. Default = true
%           .Verbose     - Print progress to command window
%       dxdTSol   - A structure vector containing the solution to the model
%                   sensitivities under each condition. Optional, but
%                   speeds up the calculation.
%
%   Outputs:
%       F = ObjectiveInformation(m, con, obj, ...)
%           F - A matrix
%
%       [F, All] = ObjectiveInformation(m, con, obj, ...)
%           All - The fisher information matrix (FIM) is the sum of all
%                 fisher information matrices assuming there is no
%                 covariance between errors in seperate experiments. All
%                 provides the individual FIMs for each experiment.
%
%   Additional info:
%   - The experimental condition vector can also be a cell vector
%   - The objective function array can also be a cell array. Empty entries 
%   in the cell array and entries in the structure array with empty
%   values are ignored. This way, conditions can have different numbers of
%   objective functions associated with them.
%   - The optional solutions, dxdTSol, can also be cell vectors. Empty
%   entries in the cell vector and entries in the structure vector with
%   empty values are ignored.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:ObjectiveInformation:AtLeastThreeInputs', 'ObjectiveExpectedHessian requires at least 3 input arguments.')
if nargin < 5
    dxdTSol = [];
    if nargin < 4
        opts = [];
    end
end

assert(isscalar(m), 'KroneckerBio:ObjectiveInformation:MoreThanOneModel', 'The model structure must be scalar')

% Options
defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];%TODO
defaultOpts.UseModelICs    = true;
defaultOpts.UseModelInputs = false;
defaultOpts.Normalized     = true;
defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.Verbose        = 0;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
nCon = length(con);
nObj = size(obj, 1);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a matrix of linear indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Standardize structures
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, {opts.UseParams}, {opts.UseICs}, {opts.UseControls});

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

nT = nTx + nTk;

% Make solutions consistent as cell vectors
if isstruct(dxdTSol)
    dxdTSol = num2cell(dxdTSol);
end

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, false, opts.UseParams, opts.UseICs, opts.UseModelICs);

%% Loop through conditions
F = zeros(nT,nT);
Txind = nTk+1; % Stores the position in F where the first variable IC goes for each iCon
intOpts = opts;

% Initialize All array if requested
if nargout >= 2
    All = cell(nObj,nCon);
end

for iCon = 1:nCon
    % If opts.UseModelICs is false, the number of variables can change
    if opts.UseModelICs
        inTx = nTx;
        inT = nT;
    else
        intOpts.UseICs = opts.UseICs(:,iCon);
        inTx = sum(intOpts.UseICs);
        inT = nTk + inTx;
    end
    
    % Sensitivity integration if not provided
    if isempty(dxdTSol) || isempty(dxdTSol{iCon})
        if ~opts.UseModelICs
            intOpts.UseICs = opts.UseICs(:,iCon);
        end
        
        
        % Modify opts structure
        intOpts.AbsTol = opts.AbsTol{iCon};
        tGet = opts.tGet{iCon};
        
        % Integrate
        if verbose; fprintf(['Integrating sensitivities for ' con(iCon).Name '...']); end
        if opts.continuous(iCon) && opts.complex(iCon)
            sol = integrateObjSens(m, con(iCon), obj(:,iCon), intOpts);
        elseif opts.complex(iCon)
            sol = integrateSens(m, con(iCon), intOpts);
        elseif opts.continuous(iCon)
            sol = integrateObjSensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
        else
            sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
        end
        if verbose; fprintf('done.\n'); end
    else
        sol = dxdTSol{iCon};
    end
    
    % Sum all FIMs as computed by each objective function
    if opts.Normalized
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            if opts.UseModelICs
                T = [m.k(opts.UseParams); m.x0(opts.UseICs)]; % model initial conditions
                Fi = obj(iObj,iCon).Fn(sol, T);
                F = F + Fi;
            else
                T = [m.k(opts.UseParams); con(iCon).x0(opts.UseICs(:,iCon))]; % con initial conditions
                Fi = obj(iObj,iCon).Fn(sol, T);
                % Add to correct place in F
                linInd = [1:nTk, Txind:Txind+inTx-1]; % linear indexes of parameters
                F(linInd,linInd) = F(linInd,linInd) + Fi;
            end
            
            % Store FIM if requested
            if nargout >= 2
                All{iObj,iCon} = Fi;
            end
        end
    else
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            if opts.UseModelICs
                Fi = obj(iObj,iCon).F(sol);
                F = F + Fi;
            else
                Fi = obj(iObj,iCon).F(sol);
                % Add to correct place in F
                linInd = [1:nTk, Txind:Txind+inTx-1]; % linear indexes of parameters
                F(linInd,linInd) = F(linInd,linInd) + Fi;
            end
            
            % Store FIM if requested
            if nargout >= 2
                All{iObj,iCon} = Fi;
            end
        end
    end
    
    % Update condition x0 position
    if ~opts.UseModelICs
        Txind = Txind + inTx;
    end
end