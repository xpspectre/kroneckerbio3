function sol = integrateSens(m, con, opts)

% Constants
nx = m.nx;
nTk = sum(sum(opts.UseParams));
nTx = sum(sum(opts.UseICs));
nT  = nTk + nTx;

% Construct system
[der, jac] = constructSystem();

%TODO: allow for starting from steady-state

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Initial conditions [x0; vec(dxdT0)]
if opts.UseModelICs
    x0 = m.x0;
else
    x0 = con.x0;
end

% Initial effect of rates on sensitivities is 0
dxdT0 = zeros(nx, nTk); % Active rate parameters

% Initial effect of ics on sensitivities is 1 for that state
dxdx0                = zeros(nx,nTx);
dxdx0(opts.UseICs,:) = eye(nTx);

% Combine them into a vector
ic = [x0; vec([dxdT0, dxdx0])];

% Integrate [x; G; dxdT; dGdv] with respect to time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT));
sol.u = u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and dxdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        
        IT = speye(nT);

        dxdTStart = nx+1;
        dxdTEnd   = nx+nx*nT;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdTdx = @d2fdTdxSub;
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ --> x_T
            
            % Compute derivative of dxdT
            dxdTdot = dfdx(t, x, u) * dxdT + dfdT(t, x, u); % f_x * x_T + f_T --> f_T
            
            % Combine
            val = [f(t, x, u); vec(dxdTdot)];
        end
        
        % Jacobian of [x; dxdT] derivative
        function val = jacobian(t, joint, u)
            u = u(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Compute d/dx(dfdT)
            d2xdxdT = d2fdx2(t,x,u) * dxdT + d2fdTdx(t,x,u); % fx_T
            d2xdxdT = spermute132(d2xdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(t, x, u), sparse(nx,nx*nT);
                         d2xdxdT, kron(IT, dfdx(t,x,u))];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t, x, u);
            val = [val(:, opts.UseParams) zeros(nx, nTx)];
        end
        
        % Modifies d2fdkdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u);
            val = [val(:, opts.UseParams) zeros(nx*nx, nTx)];
        end
        
    end

end