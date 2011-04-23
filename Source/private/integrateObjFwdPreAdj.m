function [G, gSol, discreteTimes] = integrateObjFwdPreAdj(m, con, obj, opts)
verboseAll = max(opts.Verbose-1,0);

% Constants
nX = m.nX;
nCon = length(con);
nObj = size(obj,2);

% Initialize variables
ic = zeros(nX+1,1);
G = 0;
discreteTimes = cell(nCon,1);
gSol(nCon,1) = struct('solver', [], 'extdata', [], 'x', [], 'y', [], 'xe', [], 'ye', [], 'ie', [], 'stats', [], 'idata', []);

if opts.Verbose; disp('Integrating forward...'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % *Run simulation for each experiment*
    [gDer, gJac, gEvents] = constructSystem(m, iCon);
    
    %TODO: allow for starting from steady-state
    
    % Determine initial conditions
    if opts.UseModelICs
        ic(1:nX) = m.ic;
    else
        ic(1:nX) = con(iCon).ic;
    end
    
    % Reset to zero initial condition of continuous objective function
    ic(nX+1) = 0;
    
    % Get the times that the solver will need to stop
    stopTimes = unique([0; con(iCon).tF; con(iCon).tStop((con(iCon).tStop < con(iCon).tF))]);
    
    % Integrate [x; g] over time
    gSol(iCon) = accumulateSol(gDer, gJac, gEvents, [], ic, con(iCon), opts.AbsTol{iCon}(1:nX+1), opts.RelTol, stopTimes, 1, 1:m.nX);
                
    % Extract continuous term
    contG = gSol(iCon).y(nX+1,end);
    
    % Compute discrete term
    discG = 0;
    for iObj = 1:nObj
        [iDiscG, temp] = obj(iCon,iObj).G([], gSol, con(iCon).u, iCon);
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
    function [gDer, gJac, gEvents] = constructSystem(m, iCon)
        
        f      = m.f;
        dfdx   = m.dfdx;
        
        gDer = @der;
        gJac = @jac;
        if ~isempty(con(iCon).events)
            events  = con.events.events;
            gEvents = @eventFun;
        else
            gEvents = [];
        end
        
        % Derivative of [x; g] with respect to time
        function val = der(t, joint, u)
            u = u(t);            
            x = joint(1:nX);
            
            % Sum continuous objective functions
            g = 0;
            for i = 1:nObj
                g = g + opts.ObjWeights(iCon,i) * obj(iCon,i).g(t, x, u);
            end
            
            val = [f(t, x, u); g];
        end
        
        % Jacobian of [x; g] derivative
        function val = jac(t, joint, u)
            u = u(t);
            x = joint(1:nX);
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nX);
            for i = 1:nObj
                dgdx = dgdx + opts.ObjWeights(iCon,i) * obj(iCon,i).dgdx(t, x, u);
            end
            
            val = [dfdx(t, x, u), sparse(nX,1);
                            dgdx,            0];
        end
        
        % Event function for system
        function [val, dir, term] = eventFun(t, joint, u)
            u = u(t);
            x = joint(1:nX);
            [val, dir, term] = events(t, x, u);
        end
    end

end