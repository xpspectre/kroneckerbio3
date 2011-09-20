function sol = integrateLna(m, con, opts)

% Constants
nr = m.nr;
nx = m.nx;

upperVInd = upperInd(nx);
lowerVInd = lowerInd(nx);

% Construct system
[der, jac] = constructSystem();

% Initial conditions
% V is symmetric. Only integrate the upper half of the matrix
if ~con.SteadyState
    if opts.UseModelICs
        ic = [m.x0; opts.V0(upperVInd)];
    else
        ic = [con.x0; opts.V0(upperVInd)];
    end
else
    error('Steady-state not yet implemented for LNA')
    ic = steadystateSys(m, con, opts);
end

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate x over time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+numel(upperVInd)));
sol.u = u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating the linear noise approximation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        S    = m.S;
        r    = m.r;
        drdx = m.drdx;
        f    = m.f;
        dfdx = m.dfdx;
        d2fdx2 = m.d2fdx2;
        
        Ix = speye(nx);
        
        VStart = nx+1;
        VEnd   = nx+numel(upperVInd);
        diagVInd = sub2ind([nx,nx], vec(1:nx), vec(1:nx));
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [x;V] with respect to time
        function val = derivative(t, joint, u)
            x = joint(1:nx); % x_
            V = zeros(nx,nx);
            V(upperVInd) = joint(VStart:VEnd); % x_x
            V(lowerVInd) = joint(VStart:VEnd); % x_x
            u   = u(t); % u_
            
            % xdot = f
            % Vdot = dfdx *{x.x} V +{f+V.x;V.x+f} dfdx *{x.x} V +{f+x;f+x} S *{rxr} r *{r.r} S
            Vdot = dfdx(t,x,u) * V; % f_x * v1_V.x -> f_V.x
            Vdot = Vdot + Vdot.' + bsxfun(@times, S, r(t,x,u).') * S.'; % f_V.x + V.x_f + x_r .* r * r_x -> x_x
            
            val = [f(t,x,u); Vdot(upperVInd)];
        end
        
        % Jacobian of x derivative
        function val = jacobian(t, joint, u)
            x = joint(1:nx); % x_
            V = zeros(nx,nx);
            V(upperVInd) = joint(VStart:VEnd); % x_x
            V(lowerVInd) = joint(VStart:VEnd); % x_x
            u   = u(t); % u_
            
            % dxdot/dx = dfdx
            % dxdot/dV = 0
            % dVdot/dx = d2fdx2 *{x1.x} V +{f+V.x;V.x+f;x2+x2} d2fdx2 *{x1.x} V +{f+S.x;f+S.x;x2+drdx.x} S *{rxr} drdx *{r.r} S
            % dVdot/dV = dfdx *{x.x} Iv +{f+V1.x;V1.x+f} dfdx *{x.x} dfdx
            dVdotdx = d2fdx2(t,x,u) * V; % fx_x * V.x_V.x -> fx_V.x
            dVdotdx = spermute132(dVdotdx, [nx,nx,nx], [nx*nx,nx]) + reshape(dVdotdx.', [nx*nx,nx]); % (fx_V.x -> f,V.x_x) + (fx_V.x -> V.x_fx -> V.x,f_x) -> ff_x
            dVdotdx = dVdotdx + reshape(S * reshape(bsxfun(@times, vec(S.'), repmat(drdx(t,x,u), nx,1)), nr,nx*nx), nx*nx,nx); % ff_x + (f_r * ((f_r -> r_f -> rf_) .* (r_x -> rf_x) -> rf_x -> r_fx) -> f_fx -> ff_x) -> ff_x
            dVdotdx = dVdotdx(upperVInd,:);
            
            dVdordV = kron(Ix, dfdx(t,x,u)) + kron(dfdx(t,x,u), Ix) + sparse(repmat(vec(1:nx*nx), nx,1), vec(repmat(1:nx*nx, nx,1)), vec(repmat(dfdx(t,x,u), nx,1)), nx*nx,nx*nx) + sparse(vec(repmat(1:nx*nx, nx,1)), repmat(vec(1:nx*nx), nx,1), vec(repmat(dfdx(t,x,u).', nx,1)), nx*nx,nx*nx);
            dVdordV(:,diagVInd) = dVdordV(:,diagVInd) ./ 2; % Diagonal does not get copied, so it has half the normal effect
            dVdordV = dVdordV(upperVInd,upperVInd);
            
            val = [dfdx(t,x,u), zeros(nx, numel(upperVInd));
                   dVdotdx,     dVdordV];
        end
    end
end