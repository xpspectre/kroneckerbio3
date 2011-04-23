function [G, D] = integrateObjSensAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);
if verboseAll; tic; end

% Run forward method
[G, xSol, discreteTimes] = integrateObjFwdPreAdj(m, con, obj, opts);

% Constants
nX = m.nX;
nVP = sum(opts.UseParams);
nVX = sum(sum(opts.UseICs));
nV  = nVP + nVX;
nCon = length(con);
nObj = size(obj,2);

% Initialize variables
D = zeros(1,nV);

if opts.Verbose; disp('Integrating adjoint:'); end
for iCon = 1:nCon
    % *Run adjoint method*
    [jDer, jJac, jDelta] = constructSystem(m, iCon);
    
    % Set initial conditions
    j0 = zeros(1, nX+nVP);
    
    % Get the times that the solver will need to stop
    stopTimes = unique([0; con(iCon).tF; discreteTimes{iCon}; con(iCon).tStop((con(iCon).tStop < con(iCon).tF))]);
    
    % Integrate [lambda; D] backward in time
    sol = accumulateSol(jDer, jJac, [], jDelta, j0, con(iCon), opts.AbsTol{iCon}(nX+2:nX+1+nX+nV), opts.RelTol, stopTimes, -1);
    
    % *Add contributions to derivative*
    % (Subtract contributions because sensitivities were integrated
    % backward)
    
    % Rate parameters (Subtract contributions because sensitivities were
    % integrated backward)
    curD = zeros(1,nV);
    curD(1:nVP) = curD(1:nVP) - sol.y(nX+1:end,end).';
    
    % Initial conditions
    if opts.UseModelICs
        lambda = sol.y(1:nX,end).';
        dx0dx0 = eye(nX,nX);
        curD(nVP+1:nV) = curD(nVP+1:nV) + lambda * dx0dx0(:,opts.UseICs);
    else
        lambda = sol.y(1:nX,end).';
        dx0dx0 = eye(nX,nX);
        curD(nVP+1:nV) = curD(nVP+1:nV) + lambda * dx0dx0(:,opts.UseICs(:,iCon));
    end
    
    % Add to cumulative goal value
    D = D + curD;
    
    if verboseAll; fprintf('iCon = %d\t||dGdp|| = %g\tTime = %0.2f\n', iCon, norm(curD), toc); end    
end

if opts.Verbose;fprintf('Summary: ||dGdp|| = %g\n', norm(D));end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating lambda and D %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [jDer, jJac, jDelta] = constructSystem(m, iCon)
        dfdx    = m.dfdx;
        dfdp    = @dfdpSub;
        
        jDer    = @der;
        jJac    = @jac;
        jDelta  = @delta;
        
        % Derivative of [lambda; -D] with respect to time
        function val = der(t, joint, u)
            u = u(t);
            x = deval(xSol(iCon), t, 1:nX);
            l  = joint(1:nX);
            
            % Sum continuous objective functions
            dG = zeros(nX+nVP,1);
            for i = 1:nObj
                dgdp = obj(iCon,i).dgdp(t, x, u);
                dgdp = dgdp(:,opts.UseParams);
                dG = dG + opts.ObjWeights(iCon,i)*[-obj(iCon,i).dgdx(t, x, u), dgdp].';
            end
            
            dF = [-dfdx(t, x, u),  dfdp(t, x, u)].';
            val = dF*l + dG;
        end
        
        % Jacobian of [lambda; -D] derivative
        function val = jac(t, joint, u)
            u = u(t);
            x = deval(xSol(iCon), t, 1:nX);
            
            val = [-dfdx(t, x, u).', sparse(nX, nV);
                    dfdp(t, x, u).', sparse(nV, nV)];
        end
        
        % Discrete effects of the objective function
        function val = delta(t, u)
            u = u(t);
            
            val = zeros(1,nX+nVP);
            for i = 1:nObj
                dGdp = obj(iCon,i).dGdp(t, xSol, u, iCon);
                dGdp = dGdp(:,opts.UseParams);
                val = val + opts.ObjWeights(iCon,i)*[-obj(iCon,i).dGdx(t, xSol, u, iCon), dGdp];
            end
        end
        
        % Modifies dfdp to relate only to the parameters of interest
        function val = dfdpSub(t, x, u)
            val = m.dfdp(t, x, u);
            val = val(:,opts.UseParams);
        end
    end

end