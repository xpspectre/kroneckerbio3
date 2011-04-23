function [optimalGRatio optimalDRatio] = integrateOptimalAbsTolSensSimple(m, con, obj, opts)
% There can only be one con

%% Constants
nX = m.nX;
nVP = sum(opts.UseParams);
nVX = sum(sum(opts.UseICs));
nV  = nVP + nVX;
nObj = numel(obj);

verbose = opts.Verbose;
opts.Verbose = max(opts.Verbose-1,0);

absTolStart = nX+nX*nV+1;
absTolEnd   = nX+nX*nV+nX*nX+nX*nV*nX;
dxdvStart   = nX+1;
dxdvEnd     = nX+nX*nV;
dxdx0Start  = 1;
dxdx0End    = nX*nX;
dsdx0Start  = nX*nX+1;
dsdx0End    = nX*nX+nX*nV*nX;

%% Integrate the basic system
% Complete integration required even if objective function
if verbose; fprintf('Integrating initial sensitivities...'); end
sol = integrateSens(m, con, opts);
if verbose; fprintf('done.\n'); end

nt = numel(sol.x) - 1; % First point has no error

%% Chop down coverage of integration
if opts.Coverage <= 1
    nPick = ceil(opts.Coverage * nt);
elseif opts.Coverage > 1
    nPick = ceil(opts.Coverage);
end
selectIndex = ceil((0:nPick-1) * nt / nPick) + 2; % Add 1 to undo 0 index, add 1 more to skip first point

if verbose; fprintf('Covering %d out of %d times...\n', nPick, nt); end

%% Objective sensitivity to species at all times
% Discrete Times
discreteTimes = [];
for iObj = 1:nObj
    discreteTimes = [discreteTimes; vec(obj(iObj).discreteTimes)];
end
discreteTimes = unique(discreteTimes);
nDisc = numel(discreteTimes);

% Precompute dG/dx at all times
dGdx = zeros(nX,nDisc);
for iDisc = 1:nDisc
    for iObj = 1:nObj
        dGdx(:,iDisc) = dGdx(:,iDisc) + opts.ObjWeights(iObj) * vec(obj(iObj).dGdx(discreteTimes(iDisc), sol, con, 1));
    end
end

%% Sensitivities at all timepoints
effect = zeros(nX+nV*nX,nPick);

% Initial conditions are the same for all runs
ic = sparse(nX*nX+nX*nX*nV,1);
ic(dxdx0Start:dxdx0End) = speye(nX);

for it = 1:nPick
    % Construct system
    [jDer, jJac] = constructSystem();
    
    % Set the integration intervals
    t0 = sol.x(selectIndex(it));
    Disc0 = find(discreteTimes >= t0, 1); % First index in discreteTimes that matters
    stopTimes = unique([t0; con.tF; con.tStop((con.tStop > t0) & (con.tStop < con.tF))]);
    
    % Integrate [dx/dx0; d2x/dx0dv]
    if verbose; fprintf('Integrating from time %d...', t0); end
    if numel(stopTimes) > 1
        % There is actual integration to be done
        soli = accumulateSol(jDer, jJac, [], [], ic, con, opts.AbsTol(absTolStart:absTolEnd), opts.RelTol, stopTimes, 1, [], discreteTimes(Disc0:nDisc));
    else
        % This is the final timepoint, no need to integrate
        soli.x = stopTimes;
        soli.y = ic;
    end
    if verbose; fprintf('done.\n'); end
    
    % Compute the effect over remaining time
    for iu = Disc0:nDisc % index into discreteTimes
        effect(:,it) = effect(:,it) + vec(dGdx(:,iu).' * reshape(soli.y(:,iu-Disc0+1), nX,nX+nV*nX)); % multiply by dG/dx, and sum over time
    end
end

% The effect is absolute
effect = abs(effect);

clear soli % very large matrix

%% Maximum impact
maximpact = zeros(1+nV,1);
for it = 1:nPick
    % The impact is absolute
    x = abs(sol.y(1:nX,selectIndex(it))); % x_
    xbox = spdiags(x,0,nX,nX); % x_x diagonal
    dxdv = abs(reshape(sol.y(nX+1:nX+nX*nV,selectIndex(it)),nX,nV)); % x_v
    
    xeffect = effect(1:nX,it);
    xeffectbox = spdiags(xeffect,0,nX,nX);
    dxdpeffect = reshape(effect(nX+1:nX+nX*nV),nV,nX);
    
    % Impact that x has on G
    maximpact(1) = max([maximpact(1); xeffect .* x]); % x_ .* x_ --> x_
    
    % Impact that x has on dG/v
    maximpact(2:end) = max(maximpact(2:end), max(dxdpeffect*xbox, [],2)); % (vx_ --> v_x) * x_x --> v_x --> v_
    
    % Impact that dx/dv has on dG/dv
    maximpact(2:end) = max(maximpact(2:end), vec(max(xeffectbox*dxdv, [],1))); % x_x * x_v --> x_v --> _v --> v_
end

%% Optimal impact from ratio
% Take the maximum over time beforehand
effect = max(effect, [], 2);

% Impact ratio for x's impact on G
optimalGRatio = maximpact(1) ./ effect(1:nX);

% Initialize D
optimalDRatio = zeros(nX+nX*nV,1);

% Impact ratio for x's impact on dG/dv
% There can only be one abstol, take the most conservative
optimalDRatio(1:nX) = min(optimalGRatio, min(spdiags(maximpact(2:end), 0,nV,nV) * reshape(effect(nX+1:nX+nV*nX), nV,nX).^-1, [],1).');

% Impact ratio for dx/dp's impact on dG/dp
% There is no need to take the maximum over v because the entries are
% unique along v. dx/dv only affects dG/dv of the same v
optimalDRatio(nX+1:nX+nX*nV) = (maximpact(2:end) * (effect(1:nX).^(-1)).').'; % v_ * (x_ --> _x) --> v_x --> x_v --> xv_

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [jDer, jJac] = constructSystem()
        
        Ix = speye(nX);
        Ipx = speye(nV*nX);
    
        dfdx    = m.dfdx;
        d2fdx2  = m.d2fdx2;
        d2fdvdx = @d2fdvdxSub;
        
        jDer = @der;
        jJac = @jac;
        
        % Derivative of [dx/dx0; d2x/dx0dv] with respect to time
        function val = der(t, joint, u)
            u = u(t);
            jointold = deval(sol, t, 1:dxdvEnd); % x+xv_
            x = jointold(1:nX); % x_
            dxdv = reshape(jointold(dxdvStart:dxdvEnd), nX,nV); % x_v
            dxdx0 = reshape(joint(dxdx0Start:dxdx0End), nX,nX); % x_x0
            dsdx0 = reshape(joint(dsdx0Start:dsdx0End), nX,nV*nX); % x_vx0
            
            % Derivative of dx/dx0
            dxdx0dot = vec(dfdx(t,x,u) * dxdx0); % x_x * x_x0 --> x_x0 --> xx0_
            
            % Derivative of ds/dx0
            dsdx0dot = d2fdx2(t,x,u) * dxdv + d2fdvdx(t,x,u); % fx_x * x_v + fx_v --> fx_v
            dsdx0dot = full(dsdx0dot); % fx_v
            dsdx0dot = reshape(dsdx0dot, nX,nX,nV); % fx_v --> f_x_v
            dsdx0dot = permute(dsdx0dot, [1,3,2]); % f_x_v --> f_v_x
            dsdx0dot = reshape(dsdx0dot, nX*nV,nX); % f_v_x --> fv_x
            dsdx0dot = sparse(dsdx0dot); % fv_x
            dsdx0dot = vec(dsdx0dot * dxdx0) + vec(dfdx(t,x,u) * dsdx0); % (fv_x * x_x0 --> fv_x0 --> fvx0_) + (f_x * x_vx0 --> f_vx0 --> fvx0_) --> fvx0_
            
            % Combine all three
            val = [dxdx0dot; dsdx0dot];
        end
        
        % Jacobian of [dx/dx0; d2x/dx0dv] derivative
        function val = jac(t, joint, u)
            u = u(t);
            jointold = deval(sol, t, 1:nX+nX*nV); % x+xv_
            x = jointold(1:nX); % x_
            dxdv = reshape(jointold(dxdvStart:dxdvEnd), nX,nV); % x_v
            
            % Jacobian of ds/dx0 wrt dx/dx0
            d2sdxdx0dx0 = d2fdx2(t,x,u) * dxdv + d2fdvdx(t,x,u); % fx_x * x_v + fx_v --> fx_v
            d2sdxdx0dx0 = full(d2sdxdx0dx0); % fx_v
            d2sdxdx0dx0 = reshape(d2sdxdx0dx0, nX,nX,nV); % fx_v --> f_x_v
            d2sdxdx0dx0 = permute(d2sdxdx0dx0, [1,3,2]); % f_x_v --> f_v_x
            d2sdxdx0dx0 = reshape(d2sdxdx0dx0, nX*nV,nX); % f_v_x --> fv_x
            d2sdxdx0dx0 = sparse(d2sdxdx0dx0); % fv_x
            
            % Combine all
            val = [kron(Ix, dfdx(t,x,u)), sparse(nX*nX,nX*nV*nX);
                   kron(Ix, d2sdxdx0dx0), kron(Ipx, dfdx(t,x,u))];
        end
        
        % Modifies d2fdpdx to relate only to the parameters of interest
        function val = d2fdvdxSub(t, x, u)
            val = m.d2fdpdx(t,x,u);
            val = [val(:, opts.UseParams) zeros(nX*nX, nVX)];
        end
        
    end
end
