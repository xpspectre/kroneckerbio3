function [G, sol, discreteTimes] = integrateObjFwdPreAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nCon = numel(con);
nObj = size(obj,1);

% Initialize variables
ic = zeros(nx+1,1);
G = 0;
discreteTimes = cell(nCon,1);
sol = emptystruct(nCon, 'solver', 'extdata', 'x', 'y', 'xe', 'ye', 'ie', 'stats', 'idata', 'u', 'C1', 'C2', 'c');

if opts.Verbose; disp('Integrating forward...'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % *Run simulation for each experiment*
    [der, jac] = constructSystem();
    
    %TODO: allow for starting from steady-state
    
    % Determine initial conditions
    if opts.UseModelICs
        ic(1:nx) = m.x0;
    else
        ic(1:nx) = con(iCon).x0;
    end
    
    % Reset to zero initial condition of continuous objective function
    ic(nx+1) = 0;
    
    % Input
    if opts.UseModelInputs
        u = m.u;
    else
        u = con.u;
    end
    
    % Integrate [x; g] over time
    sol(iCon) = accumulateOde(der, jac, 0, con(iCon).tF, ic, u, con(iCon).Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx));
    sol(iCon).u = u;
    sol(iCon).C1 = m.C1;
    sol(iCon).C2 = m.C2;
    sol(iCon).c  = m.c;
                
    % Extract continuous term
    contG = sol(iCon).y(nx+1,end);
    
    % Compute discrete term
    discG = 0;
    for iObj = 1:nObj
        [iDiscG, temp] = obj(iCon,iObj).G(sol(iCon));
        discreteTimes{iCon} = [discreteTimes{iCon}; vec(temp)];
        discG = discG + opts.ObjWeights(iCon,iObj) * iDiscG;
    end
    
    % Remove repetitive discreteTimes
    discreteTimes{iCon} = unique(discreteTimes{iCon});
    
    % Add to cumulative goal value
    G = G + contG + discG;
    
    if verboseAll; fprintf('iCon = %d\tG = %g\tTime = %0.2f\n', iCon, contG + discG, toc); end
end

if opts.Verbose; fprintf('Summary: G = %g\n', G); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and g %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        
        f      = m.f;
        dfdx   = m.dfdx;
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [x; g] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            
            % Sum continuous objective functions
            g = 0;
            for iObj = 1:nObj
                g = g + opts.ObjWeights(iObj,iCon) * obj(iObj,iCon).g(t,x,u);
            end
            
            val = [f(t,x,u); g];
        end
        
        % Jacobian of [x; g] derivative
        function val = jacobian(t, joint, u)
            u = u(t);
            x = joint(1:nx);
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx);
            for iObj = 1:nObj
                dgdx = dgdx + opts.ObjWeights(iObj,iCon) * obj(iObj,iCon).dgdx(t,x,u);
            end
            
            val = [dfdx(t,x,u), sparse(nx,1);
                          dgdx,            0];
        end
    end

end