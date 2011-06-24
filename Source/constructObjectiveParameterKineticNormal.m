function obj = constructObjectiveParameterKineticNormal(m, kbar, Vkbar, UseParams, normalized, name)

% Clean up inputs
if nargin < 6
    name = [];
    if nargin < 5
        normalized = [];
        if nargin < 4
            UseParams = [];
            if nargin < 3
                Vkbar = [];
                if nargin < 2
                    kbar = [];
                    if nargin < 1
                        m = [];
                    end
                end
            end
        end
    end
end

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    obj = emptystruct(m, 'Type', 'Name', 'Continuous', 'Complex', 'G', 'dGdk', 'd2Gdk2', 'F', 'Fn', 'p', 'pvalue', 'n', 'Update');
    return
end

% Check inputs
n = numel(kbar);
assert(ndims(Vkbar) == 2 && all(size(Vkbar) == [n,n]), 'KroneckerBio:constructObjectiveParameterNormal:Vsize', 'Input "Vkbar" must be a square matrix of numel(kbar)')

% Find unique times
discreteTimes = 0;

% Constants
nx = m.nx;
nk = m.nk;

% Fix UseParams
[UseParams, nTk] = fixUseParams(UseParams, nk);

% Pre-select kinetic parameters
k = m.k;
Tk = k(UseParams);
Tkbar = kbar(UseParams);
VTbark = Vkbar(UseParams,UseParams);

% Normalization
if normalized
    Tk = log(Tk);
    VTbark = spdiags(Tkbar.^(-1),0,nTk,nTk) * VTbark * spdiags(Tkbar.^(-1),0,nTk,nTk);
    Tkbar = log(Tkbar);
end

FTbark = infoinv(VTbark);

% Gut m
temp = m;
clear m
m.nx = temp.nx;
m.nk = temp.nk;
m.k = temp.k;
clear temp

if isempty(name)
	name = 'PriorKineticParametersNormal';
end

% Objective structure
obj.Type = 'Objective.Information';
obj.Name = name;

obj.Continuous    = false;
obj.Complex       = false;
obj.Linked        = 0;
obj.DiscreteTimes = discreteTimes;

obj.G      = @G;
obj.dGdk   = @dGdk;
obj.d2Gdk2 = @d2Gdk2;

obj.F      = @F;
obj.Fn     = @Fn;
obj.p      = @p;
obj.pvalue = @pvalue;
obj.n      = nTk;

obj.Update = @update;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val stopTimes] = G(sol)
        stopTimes = discreteTimes;
        diff = Tk - Tkbar;
        val = diff.' * (VTbark \ diff);
    end

    function val = dGdk(t,sol)
        val = zeros(nk,1);
        if normalized
            val(UseParams) = 2 * diag(k(UseParams).^(-1)) * (VTbark \ (Tk - Tkbar));
        else
            val(UseParams) = 2 * (VTbark \ (Tk - kbar));
        end
    end

    function val = d2Gdk2(t,sol)
        val = zeros(nk,nk);
        if normalized
            val(UseParams,UseParams) = 2 * diag(k(UseParams).^(-1)) * (VTbark \ diag(k(UseParams).^(-1))) - 2 * diag(k(UseParams).^(-2)) * diag(VTbark \ (Tk - kbar));
        else
            val(UseParams,UseParams) = 2 * inv(VTbark);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Information theory %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Likelihood function
    function p = p(sol)
        diff = Tk - Tkbar;
        p = (2*pi).^(-nTk/2) * det(VTbark).^(-1/2) * exp(-1/2 * diff.' * (VTbark \ diff));
    end

%% Fisher information
    function val = F(dxdTSol)
        % This objective only provides information on the first nTk parameters
        nT = (size(dxdTSol.y, 1) - nx) / nx;
        val = zeros(nT,nT);
        
        if ~normalized
            val(1:nTk,1:nTk) = FTbark;
        else
            % Unnormalize
            TkbarTemp = exp(Tkbar);
            val(1:nTk,1:nTk) = diag(TkbarTemp) * FTbark * diag(TkbarTemp);
        end
    end

    function val = Fn(dxdTSol, T)
        % This objective only provides information on the first nTk parameters
        nT = (size(dxdTSol.y, 1) - nx) / nx;
        val = zeros(nT,nT);
        
        if ~normalized
            % Normalize
            val(1:nTk,1:nTk) = diag(Tkbar.^(-1)) * FTbark * diag(Tkbar.^(-1));
        else
            val(1:nTk,1:nTk) = FTbark;
        end
    end

%% P-Value
    function pval = pvalue(sol, obj)
        chi2 = 0;
        ntot = 0;
        
        % Sum chi-square and n across all data
        for i = 1:length(sol)
            chi2 = chi2 + obj(i).G(sol(i));
            ntot = ntot + obj(i).n;
        end
        
        pval = chi2pvalue(chi2, ntot);
    end

%% Update
    function objNew = update(m, con, UseParams, UseICs, UseControls)
        objNew = constructObjectiveParameterKineticNormal(m, kbar, Vkbar, UseParams, normalized);
    end
end