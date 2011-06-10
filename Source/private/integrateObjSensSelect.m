function sol = integrateObjSensSelect(m, con, obj, tGet, opts)

% Constants
nx = m.nx;
nTk = sum(sum(opts.UseParams));
nTx = sum(sum(opts.UseICs));
nT  = nTk + nTx;
nObj = numel(obj);

% Construct system
[der, jac, events] = constructSystem();

%TODO: events

if ~con.SteadyState
    % Initial conditions [x0; G0; vec(dxdT0); vec(dGdv0)]
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
    ic = [x0; 0; vec([dxdT0, dxdx0]); zeros(nT,1)];
else
    % Run to steady-state first
    ic = steadystateSens(m, con, opts);
    
    % Include objective and gradient
    ic = [ic(1:nx); 0; ic(nx+1:nx+1+nx*nT); zeros(nT,1)];
end

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate [x; G; dxdT; dGdv] with respect to time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+nx*nT), [], [], [], [], [], tGet);
sol.u = u(tGet);
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x, G, dxdp, and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, jEvents] = constructSystem()
        
        IT = speye(nT);

        dxdTStart = nx+1+1;
        dxdTEnd   = nx+1+nx*nT;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdTdx = @d2fdTdxSub;
        
        der = @derivative;
        jac = @jacobian;
        if false%~isempty(con.events)
            %TODO: events
            events  = obj.Events;
            jEvents = @eventFun;
        else
            jEvents = [];
        end
        
        % Derivative of [x; G; dxdT; dGdT] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % xT_ --> x_T
            
            % Compute derivative of dxdT
            dxdTdot = dfdx(t, x, u) * dxdT + dfdT(t, x, u); % f_x * x_T + f_T --> f_T
            
            % Sum continuous objective functions
            g = 0;
            dgdT = zeros(nT,1);
            for i = 1:nObj
                g = g + opts.ObjWeights(i) * obj(i).g(t, x, u);
                temp = vec(obj(i).dgdk(t, x, u)); % k_
                temp = [temp(:,opts.UseParams); sparse(nTx,1)]; % k_ --> T_   % remove rows for inactive parameters
                dgdT = dgdT + opts.ObjWeights(i) * (vec(vec(obj(i).dgdx(t, x, u)).' * dxdT) + temp); % T_ + (_x * x_T --> _T --> T_) + T_ --> T_
            end
                        
            val = [f(t, x, u); g; vec(dxdTdot); dgdT];
        end
        
        % Jacobian of [x; G; dxdT; dGdv] derivative
        function val = jacobian(t, joint, u)
            u = u(t);            
            x = joint(1:nx); % x_
            dxdT = reshape(joint(dxdTStart:dxdTEnd), nx,nT); % x_T
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx); % _x
            d2gdxdT = sparse(nT,nx); % T_x
            for i = 1:nObj
                dgdx = dgdx + opts.ObjWeights(i) * vec(obj(i).dgdx(t, x, u)).'; % _x
                
                % Compute d/dx(dgdT)
                temp = obj(i).d2gdxdk(t,x,u); % Partial d2gdpdx k_x
                temp = [temp(opts.UseParams, :); sparse(nTx,nx)]; % k_x --> T_x
                d2gdxdT = d2gdxdT + opts.ObjWeights(i) * (dxdT.' * obj(i).d2gdx2(t,x,u) + temp); % T_x + T_x * x_x + T_x --> T_x
            end
            
            % Compute d/dx(dfdT)
            d2xdxdT = d2fdx2(t,x,u) * dxdT + d2fdTdx(t,x,u); % fx_T
            d2xdxdT = spermute132(d2xdxdT, [nx,nx,nT], [nx*nT,nx]);
            
            % Combine
            val = [dfdx(t, x, u), sparse(nx,1+nx*nT+nT);
                            dgdx, sparse(1,1+nx*nT+nT);
                         d2xdxdT, sparse(nx*nT,1), kron(IT, dfdx(t,x,u)), sparse(nx*nT,nT);
                         d2gdxdT, sparse(nT,1), kron(IT, dgdx), sparse(nT,nT)];
        end
        
        % Modifies dfdp to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t, x, u);
            val = [val(:, opts.UseParams) zeros(nx, nTx)];
        end
        
        % Modifies d2fdpdx to relate only to the parameters of interest
        function val = d2fdTdxSub(t, x, u)
            val = m.d2fdkdx(t,x,u);
            val = [val(:, opts.UseParams) zeros(nx*nx, nTx)];
        end
        
        % Event function for system
%         function [val, dir, term] = eventFun(t, joint, u)
%             u = u(t);
%             x = joint(1:nx);
%             [val, dir, term] = events(t, x, u);
%         end
    end

end