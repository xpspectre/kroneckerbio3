function absTolRatio = GoodAbsTol(m, con, sd, opts)
%GoodAbsTol Make a reasonable estimate as what the absolute tolerance on
%   the species and sensitivities for objective calculations
%
%   absTolRatio = GoodAbsTol(m, con, sd, opts)
%
%   This function uses the lowest uncertainty in the outputs, as given by
%   sd at an output value of 0, as the floor on the relavent output value.
%   This output value is then propogated to a floor on the relavent species
%   values. This floor is the AbsTolRatio. AbsTolRatio times RelTol is a
%   good guess as to what should be used as AbsTol for ODE integration of
%   the states and sensitivities.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model
%   con: [ experiment struct vector ]
%       The experimental conditions 
%   sd: [ handle @(t,yInd,yVal) returns positive scalar ]
%       The standard deviation as a function of the time, output index, and
%       output value.
%   opts: [ options struct scalar ]
%       Optional
%       .UseModelICs [ logical scalar {false} ]
%           Indicates that the model's initial conditions should be used
%           instead of those of the experimental conditions. This will
%           determine both which parameters are used for simulation as well
%           as what parameters will be varied in the optimization.
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions. This will determine both
%           which parameters are used for simulation as well as what
%           parameters will be varied in the optimization.
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters that will be allowed to vary
%           during the optimization
%       .UseICs [ logical matrix nx by nCon | logical vector nx |
%                 positive integer vector {[]} ]
%           Indicates the initial conditions of the state species that will
%           be allowed to vary during the optimzation. If UseModelICs is
%           true then UseICs can be a vector of linear indexes or a vector
%           of logicals length of nx. If UseModelICs is false then UseICs
%           can be a matrix of logicals size nx by nCon. It can also be a
%           vector of length nx, and every experiment will be considered to
%           have the same active IC parameters. It can also be a vector of
%           linear indexes into the nx vector and assumed the same for all
%           conditions.
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%
%   Outputs
%   absTolRatio: [ postive vector nx+nx*nT ]
%       Multiply this vector by RelTol to get the good AbsTol

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Constants
nx = m.nx;
ny = m.ny;
nk = m.nk;
nCon = numel(con);

% Default options
defaultOpts.UseParams               = 1:m.nk;
defaultOpts.UseICs                  = [];
defaultOpts.UseControls             = [];%TODO
defaultOpts.UseAdjoint              = false;

defaultOpts.TolOptim                = 1e-5;

opts = mergestruct(defaultOpts, opts);

% Ensure UseRates is column vector of logical indexes
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a matrix of logical indexes
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

nT = nTx + nTk;

%% xGoodRatio
% Get the floor on the uncertainty of the outputs
yFloor = zeros(ny,1);
for iy = 1:ny
    yFloor(iy) = sd(0, iy, 0);
end

% Transfer the floor back to the species
% Out of all the outputs that a species affects, we should take the most
% conservative floor.
xFloor = bsxfun(@rdivide, yFloor, m.C1); % Propogation of floors
xFloor = min(xFloor, [], 1); % Keep most conservative
xFloor = vec(xFloor); % Straighten it

% Assume that the effect between species is one-to-one and can travel an
% infinite distance along the network.
% net = NetworkDistance(m); % Distance between species
% net = spones(net); % Convert to one-to-one mapping
% xGoodRatio = bsxfun(@rdivide, xFloor, net); % Propogation
xGoodRatio = repmat(min(xFloor, [], 1), nx,1); % Keep most conservative and assign it to all
xGoodRatio = vec(xGoodRatio); % Straighten it

%% Assemble the vectors for each experiment
if ~opts.UseAdjoint
    absTolRatio = cell(nCon,1);
    for iCon = 1:nCon
        % Experiment specific parameters
        if opts.UseModelICs
            T = [m.k(opts.UseParams); m.x0(opts.UseICs)]; % model initial conditions
        else
            T = [m.k(opts.UseParams); con(iCon).x0(opts.UseICs(:,iCon))]; % con initial conditions
        end
        
        % dxdTGoodRatio
        % If the model is linear, then the maximum value of dx/dT is x/T. A good
        % ratio would therefore be x/T
        dxdTGoodRatio = bsxfun(@rdivide, xGoodRatio, T.'); % x_T
        dxdTGoodRatio = vec(dxdTGoodRatio); % xT_
        
        % Combine into one AbsTol ratio
        absTolRatio{iCon} = [xGoodRatio; dxdTGoodRatio];
    end
else%UseAdjoint
    absTolRatio = cell(nCon,1);
    for iCon = 1:nCon
        % Experiment specific parameters
        if opts.UseModelICs
            UseICsi = opts.UseICs;
            T = [m.k(opts.UseParams); m.x0(UseICsi)]; % model initial conditions
        else
            UseICsi = opts.UseICs(:,iCon);
            T = [m.k(opts.UseParams); con(iCon).x0(UseICsi)]; % con initial conditions
        end
        
        % D ratio
        % These elements are not computed in log space so their absolute
        % value will be inversely proportational to theta. Let the floor we
        % care about be TolOptim.
        DGoodRatio = opts.TolOptim ./ T;
        
        % lambda ratio
        % The effect that lambda has on D is through df/dT. Take the most
        % conservative effect. If the model is linear, then the maximum
        % value of df/dT is x/T. A good ratio would therefore be x/T
        dxdTGoodRatio = bsxfun(@rdivide, xGoodRatio, T.'); % x_T
        effect = bsxfun(@rdivide, DGoodRatio, dxdTGoodRatio.'); % T / T_x
        lAbsTol = max(effect, [], 1);
        
        % Zero entries are for those with no effect
        lAbsTol(lAbsTol == 0) = inf;
        
        % Combine into one AbsTol ratio
        absTolRatio{iCon} = [xGoodRatio; vec(lAbsTol); DGoodRatio];
    end
end