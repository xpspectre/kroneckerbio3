function [varargout] = SimulateSensitivitySelect(m, con, tGet, opts)
%SimulateSensitivitySelect integrates the sensitivities of every species with
%   respect to every parameter over all time in the mass action kinetics
%   framework and returns the values at select time points
%
%   Mathematically: dx/dT = Integral(df/dx * dx/dT + df/dT, t=0:tF)
%   
%   [...] = SimulateSensitivitySelect(m, con, tGet, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   tGet: [ nonegative vector ]
%       Indicates which time points will be returned. This does not need to
%       be sorted. Times larger than con.tF will return NaN for all values.
%   opts: [ options struct scalar ]
%       Optional
%       .UseModelICs [ logical scalar {false} ]
%           Indicates that the model's initial conditions should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters whose sensitivities are
%           desired
%       .UseICs [ logical matrix nx by nCon | logical vector nx |
%                 positive integer vector {[]} ]
%           Indicates the initial conditions of the state species whose
%           sensitivites are desired
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters whose sensitivites are
%           desired
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%   SimulateSensitivitySelect(m, con, opts)
%   	Plots the sensitivities under each condition
%
%   sim = SimulateSensitivitySelect(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t tGet
%       .dydT [ matrix ny*nT by numel(tGet) ]
%           The value of the sensitivites of the outputs at each selected
%           time point
%       .dxdT [ matrix nx by numel(tGet) ]
%           The value of the sensitivities of the states at each selected
%           time point
%       .sol [ struct scalar ]
%           The discrete integrator solution to the system
%       
%   The following are permitted if nCon = 1
%   [t, dydT] = SimulateSensitivitySelect(m, con, opts)
%     	sim.t, sim.dydT
%
%   [t, dydT, dxdT] = SimulateSensitivitySelect(m, con, opts)
%     	sim.t, sim.dydT, sim.dxdT
%
%   [t, dydT, dxdT, sol] = SimulateSensitivitySelect(m, con, opts)
%     	sim.t, sim.dydT, sim.dxdT, sim.sol
%       
%	Special
%   This function can also be requested to return an empty array of
%   structures with the same fields as simulation. This may be necessary to
%   initialize an array that the user will later fill.
%
% 	empty = SimulateSensitivitySelect(m)
%   	For integer m, create an empty simulation array m by 1
%
%   empty = SimulateSensitivitySelect([m,n,p...])
%       For integers array, create an empty simulation array size
%       [m,n,p...]

% (c) 2010 David R Hagen, Joshua F Apgar, Jared E Toettcher, & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargout <= 4, 'KroneckerBio:SimulateSensitivity:FourOrFewerOutputs', 'SimulateSensitivity must have between 0 and 4 outputs.')
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
    varargout{1} = emptystruct(m, 'Type', 'Name', 't', 'dydT', 'dxdT', 'sol');
    return
end

assert(isscalar(m), 'KroneckerBio:SimulateSensitivity:MoreThanOneModel', 'The model structure must be scalar')

% Options
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;
defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];%TODO
defaultOpts.AbsTol         = NaN;
defaultOpts.RelTol         = NaN;
defaultOpts.Verbose        = 0;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ny = m.ny;
nk = m.nk;
nCon = numel(con);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a matrix of linear indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

nT = nTx + nTk;

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, false(nCon,1), nx, nCon, false, opts.UseParams, opts.UseICs, opts.UseModelICs);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'dydT', 'dxdT', 'sol');
intOpts = opts;

for iCon = 1:nCon
    % If opts.UseModelICs is false, the number of variables can change
    if opts.UseModelICs
        inTx = nTx;
        inT = nT;
    else
        inTx = sum(opts.UseICs(:,iCon));
        inT = nTk + inTx;
    end
	
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    
    % Integrate dx/dp over time
    if verbose; fprintf(['Integrating sensitivities for ' con(iCon).Name '...']); end
    sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type = 'Simulation.MassActionKineticsSensitivity.SelectPoints';
    sim(iCon).Name = [m.Name ' in ' con.Name];
    sim(iCon).t    = sol.x;
    sim(iCon).dydT = reshape(sol.C1*reshape(sol.y, nx,(inT+1)*length(sol.x)), ny*(inT+1),length(sol.x));
    sim(iCon).dxdT = sol.y;
    sim(iCon).sol  = sol;
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
        varargout{2} = sim.dydT;
    case 3
        varargout{1} = sim.t;
        varargout{2} = sim.dydT;
        varargout{3} = sim.dxdT;
    case 4
        varargout{1} = sim.t;
        varargout{2} = sim.dydT;
        varargout{3} = sim.dxdT;
        varargout{4} = sim.sol;
    otherwise
        error(nargoutchk(0, 4, nargout, 'struct'));
end
end