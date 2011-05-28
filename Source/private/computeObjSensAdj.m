function [G, D] = computeObjSensAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);
if verboseAll; tic; end

% Constants
nx = m.nx;
nTk = sum(opts.UseParams);
nTx = sum(sum(opts.UseICs));
nT  = nTk + nTx;
nCon = length(con);
nObj = size(obj,1);

% Initialize variables
G = 0;
D = zeros(nT,1);
intOpts = opts;

if opts.Verbose; disp('Integrating adjoint:'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % * Integrate system *
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    
    % Integrate
    % Do not use select methods since the solution is needed at all time
    if opts.continuous(iCon)
        xSol = integrateObj(m, con(iCon), obj(:,iCon), intOpts);
    else
        xSol = integrateSys(m, con(iCon), intOpts);
    end
    
    % Extract continuous term
    if opts.continuous(iCon)
        contG = xSol.y(nx+1,end);
    else
        contG = 0;
    end
    
    % Compute discrete term
    discG = 0;
    discreteTimes = [];
    for iObj = 1:nObj
        [iDiscG, temp] = obj(iObj,iCon).G(xSol);
        discreteTimes = [discreteTimes; vec(temp)];
        discG = discG + opts.ObjWeights(iObj,iCon) * iDiscG;
    end
    
    % Remove repetitive discreteTimes
    discreteTimes = unique(discreteTimes);
    nDisc = numel(discreteTimes);

    % Add to cumulative goal value
    G = G + contG + discG;
    
    % * Integrate Adjoint *
    [der, jac, del] = constructSystem();
    
    % Set initial conditions
    ic = zeros(nx+nTk,1);
    
    % Input
    if opts.UseModelInputs
        u = m.u;
    else
        u = con(iCon).u;
    end
    
    % Integrate [lambda; D] backward in time
    sol = accumulateOde(der, jac, 0, con(iCon).tF, ic, u, [con(iCon).Discontinuities; discreteTimes], [], opts.RelTol, opts.AbsTol{iCon}(nx+opts.continuous(iCon)+1:nx+opts.continuous(iCon)+nx+nT), del, -1, [], [], [], 0);
    
    % *Add contributions to derivative*
    % (Subtract contribution because gradient was integrated backward)
    
    % Rate parameters
    curD = zeros(nT,1);
    curD(1:nTk) = -sol.y(nx+1:end,end);
    
    % Initial conditions
    if opts.UseModelICs
        lambda = sol.y(1:nx,end);
        dx0dx0 = speye(nx,nx);
        curD(nTk+1:nT) = dx0dx0(opts.UseICs,:) * lambda;
    else
        lambda = sol.y(1:nx,end).';
        dx0dx0 = speye(nx,nx);
        curD(nTk+1:nT) = dx0dx0(opts.UseICs(:,iCon),:) * lambda;
    end
    
    % Add to cumulative goal value
    D = D + curD;
    
    if verboseAll; fprintf('iCon = %d\t|dGdp| = %g\tTime = %0.2f\n', iCon, norm(curD), toc); end    
end

if opts.Verbose; fprintf('Summary: |dGdp| = %g\n', norm(D)); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        dfdx = m.dfdx;
        dfdT = @dfdTSub;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [lambda; D] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);
            x = deval(xSol, t, 1:nx);
            l = joint(1:nx);
            
            % Sum continuous objective functions
            dgdx = zeros(nx,1);
            dgdT = zeros(nT,1);
            for iObj = 1:nObj
                dgdx = dgdx + opts.ObjWeights(iObj,iCon)*vec(obj(iObj,iCon).dgdx(t, x, u));
                dgdk = obj(iObj,iCon).dgdk(t, x, u);
                dgdT = dgdT + opts.ObjWeights(iObj,iCon)*dgdk(opts.UseParams);
            end
            
            val = [dgdx; dgdT] - [dfdx(t,x,u).'; dfdT(t,x,u).'] * l;
        end
        
        % Jacobian of [lambda; D] derivative
        function val = jacobian(t, joint, u)
            u = u(t);
            x = deval(xSol, t, 1:nx);
            
            val = [-dfdx(t, x, u).', sparse(nx, nT);
                   -dfdT(t, x, u).', sparse(nT, nT)];
        end
        
        % Discrete effects of the objective function
        function val = delta(t)
            dGdx = zeros(nx,1);
            dGdT = zeros(nT,1);
            for iObj = 1:nObj
                dGdx = dGdx + opts.ObjWeights(iObj,iCon)*vec(obj(iObj,iCon).dGdx(t, xSol));
                dGdk = obj(iObj,iCon).dGdk(t, xSol);
                dGdT = dGdT + opts.ObjWeights(iObj,iCon)*dGdk(opts.UseParams);
            end
            
            val = [dGdx; dGdT];
        end
        
        % Modifies dfdk to relate only to the parameters of interest
        function val = dfdTSub(t, x, u)
            val = m.dfdk(t, x, u);
            val = val(:,opts.UseParams);
        end
    end

end