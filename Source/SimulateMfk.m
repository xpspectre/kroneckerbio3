function [varargout] = SimulateMfk(m, con, opts)
%SIMULATE Integrate the concentration of every species over time
%
%   Mathematically: x = Integral(f)
%   
%   [...] = Simulate(m, con, opts)
%   
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - An array of experimental conditions under which the model
%              will be simulated
%       opts - Optional function options
%           UseModelICs - Logical scalar indicating that the model's
%                         initial conditions should be used instead of
%                         those of the experimental conditions.
%                         Default = true
%           RelTol      - Numerical scalar for the relative tolerance of 
%                         the integration. Default = 1e-6
%           AbsTol      - Cell vector of vectors, numerical scalar, or 
%                         numerical vector for the absolute tolerance of
%                         the integration. Default = 1e-9 
%           V0          - Initial Conditions of the Covariance matrix.
%                         Default = zeros(nx) 
%           Verbose     - Numerical scalar for how much progress 
%                         information should be displayed. Default = 0
%
%   Outputs:
%       Simulate(m, con, opts)
%        - Plots the concentrations under each condition
%
%       sim = SimulateSensitivity(m, con, opts)
%        - An array of structures with each entry being the sensitivities
%          under one of the conditions.
%           sim.t   - If opts.tGet is empty: A vector of timepoints 
%                     chosen by the ode solver. 
%                     If opts.tGet is specified: This is opts.tGet.
%           sim.y   - If opts.tGet is empty: A function handle @(t, y) that
%                     evaluates the system at some particular timepoints
%                     (t) with respect to some particular output indexes
%                     (y).
%                     If opts.tGet is specified: A matrix m.nY by length of
%                     tGet. Each column is the concentrations at the
%                     corresponding tGet.
%           sim.x   - If opts.tGet is empty: A function handle @(t, x) that
%                     evaluates the system at some particular timspoints
%                     (t) with respect to some particular species indexes
%                     (x).
%                     If opts.tGet is specified: A matrix m.nx by length of
%                     tGet. Each column is the concentrations at the
%                     corresponding tGet.
%           sim.sol - The integrator solution to the sensitivies.
%       
%       [t, y] = Simulate(m, con, opts)
%        - sim.t, sim.y
%
%       [t, y, x] = Simulate(m, con, opts)
%        - sim.t, sim.y, sim.x
%
% (c) 2010 David R Hagen, Joshua F Apgar, Jared E Toettcher, & Bruce Tidor
% This work is released under the MIT license.
%       [t, y, x, sol] = Simulate(m, con, opts)
%        - sim.t, sim.y, sim.x, sim.sol
%       
%	Special:
%       This function can also be requested to return an empty array of
%       structures with the same fields as simulation. This may be
%       necessary to initialize an array that the user will later fill.
%
%       empty = Simulate(m)
%           For integer m, create an empty simulation array m by 1
%
%       empty = Simulate(m,n)
%           For integers m and n, create an empty simulation array m by n
%
%       empty = Simulate([m,n,p...])
%           For integers array, create an empty simulation array size
%           [m,n,p...]

% (c) 2010 David R Hagen, Joshua F Apgar, Jared E Toettcher, & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargout <= 4, 'KroneckerBio:Simulate:FourOrFewerOutputs', 'Simulate must have between 0 and 4 outputs.')
if nargin < 3
    opts = [];
    if nargin < 2
        con = [];
        if nargin < 1
            m = [];
        end
    end
end

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    assert(nargout <= 1, 'KroneckerBio:Simulate:OneOutputOnEmpty', 'Simulation must have only one output when an empty simulation structure is requested.')
    temp = struct('Type', [], 'Name', [], 't', [], 'y', [], 'x', [], 'sol', []);
    if isempty(m)
        m = 1;
    end
    if ~isscalar(m)
        varargout{1}(prod(m),1) = temp;
        varargout{1} = reshape(varargout{1}, m);
        return
    end
    if isempty(con)
        con = 1;
    end
    varargout{1}(m, con) = temp;
    return
end

assert(isscalar(m), 'KroneckerBio:Simulate:MoreThanOneModel', 'The model structure must be scalar')

% Constants
nx = m.nx;
nCon = numel(con);

% Options
defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;
defaultOpts.Verbose        = 0;
defaultOpts.V0             = zeros(nx);
opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);


% RelTol
if isempty(opts.RelTol) || isnan(opts.RelTol)
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, false(nCon,1), nx + nx*nx, nCon);

%% Run integration for each experiment
sim(nCon) = struct('Type', [], 'Name', [], 't', [], 'y', [], 'x', [], 'sol', []);
intOpts = opts;

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};

    % Integrate system
    if verbose; fprintf(['Integrating system for ' con(iCon).Name '...']); end
    sol = integrateMfk(m, con(iCon), intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type   = 'Simulation.MassActionKinectics';
    sim(iCon).Name   = [m.Name ' in ' con.Name];
    sim(iCon).t      = sol.x;
    sim(iCon).y      = @(t, varargin)evaluateOutputs(sol, t, varargin{:});
    sim(iCon).x      = @(t, varargin)evaluateStates(sol, t, varargin{:});
    sim(iCon).sol    = sol;
end

%% Work-down
switch (nargout)
    case 0
        % Save the hold state of the figure
        holdState = ishold;
        % Draw each result
        for iCon = 1:nCon
            plotExperiment(m, sim(iCon));
            hold on;
        end
        % Reset the hold state
        if ~holdState
            hold off;
        end
    case 1
        varargout{1} = sim;
    case 2
        varargout{1} = sim.t;
        varargout{2} = sim.y;
    case 3
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.x;
    case 4
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.x;
        varargout{4} = sim.sol;
    otherwise
        error(nargoutchk(0, 4, nargout, 'struct'));
end

end
% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluation functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = evaluateOutputs(sol, t, ind)
if nargin < 3
    val = sol.C1 * deval(sol, t) + sol.C2 * sol.u(t) + repmat(sol.c, 1,numel(t));
else
    val = sol.C1(ind,:) * deval(sol,t) + sol.C2(ind,:) * sol.u(t) + repmat(sol.c(ind,:), 1,numel(t));
end
end

function val = evaluateStates(sol, t, ind)
if nargin < 3
    val = deval(sol, t);
else
    val = deval(sol, t, ind);
end
end
