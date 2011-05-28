function [G, D] = computeObjSens(m, con, obj, opts)
% Switch to Adjoint method if requested
if opts.UseAdjoint
    [G, D] = computeObjSensAdj(m, con, obj, opts);
    return
end

% Continue with forward method
verboseAll = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nTk = sum(sum(opts.UseParams));
nTx = sum(sum(opts.UseICs));
nT  = nTk + nTx;
nCon = length(con);
nObj = size(obj,1);

% Initialize variables
G = 0;
D = zeros(nT,1);
Txind = 1; % Stores the position in D where the first x0 parameter goes for each iCon
intOpts = opts;

if opts.Verbose; disp('Integrating sensitivities:'); end
for iCon = 1:nCon
    if verboseAll; tic; end
    
    % If opts.UseModelICs is false, the number of variables can change
    if opts.UseModelICs
        inTx = nTx;
        inT = nT;
    else
        intOpts.UseICs = opts.UseICs(:,iCon);
        inTx = sum(intOpts.UseICs);
        inT = nTk + inTx;
    end

    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    tGet = opts.tGet{iCon};
    
    % Integrate
    if opts.continuous(iCon) && opts.complex(iCon)
        sol = integrateObjSens(m, con(iCon), obj(:,iCon), intOpts);
    elseif opts.complex(iCon)
        sol = integrateSens(m, con(iCon), intOpts);
    elseif opts.continuous(iCon)
        sol = integrateObjSensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
    else
        sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
    end
    
    % *Compute G*
    % Extract continuous term
    if opts.continuous(iCon)
        contG = sol.y(nx+1,end);
    else
        contG = 0;
    end
    
    % Compute discrete term
    discG = 0;
    discreteTimes = [];
    for iObj = 1:nObj
        [iDiscG, temp] = obj(iObj,iCon).G(sol);
        discreteTimes = [discreteTimes; vec(temp)];
        discG = discG + opts.ObjWeights(iObj,iCon) * iDiscG;
    end

    % Remove repetitive discreteTimes
    discreteTimes = unique(discreteTimes);
    nDisc = numel(discreteTimes);

    % Add to cumulative goal value
    G = G + contG + discG;
    
    % *Compute D*
    % Extract continuous term
    if opts.continuous(iCon)
        dGdTStart = nx+1+nx*inT+1;
        dGdTEnd   = nx+1+nx*inT+inT;
        contD = sol.y(dGdTStart:dGdTEnd,end);
    else
        contD = sparse(nT,1);
    end
    
    % Compute discrete term
    dxdTStart = nx+opts.continuous(iCon)+1;
    dxdTEnd   = nx+opts.continuous(iCon)+nx*inT;
    dxdT = deval(sol, discreteTimes, dxdTStart:dxdTEnd); % xv_t
    discD = zeros(inT,nObj*nDisc); % v_ot
    for iObj = 1:nObj
        objDiscD = zeros(inT,1); % v_
        for iDisc = 1:nDisc
            objDiscD = objDiscD + vec(vec(obj(iObj,iCon).dGdx(discreteTimes(iDisc), sol)).' * reshape(dxdT(:,iDisc), nx, inT)); % v_ + (_x * x_v --> _v --> v_) --> v_
            temp = vec(obj(iObj,iCon).dGdk(discreteTimes(iDisc), sol)); % p_ % partial dGdp(i)
            temp = [temp(opts.UseParams); sparse(inTx,1)]; % p_ --> v_
            objDiscD = objDiscD + temp; % v_
        end
        discD(:,(iObj-1)*nDisc + iDisc) = opts.ObjWeights(iObj,iCon) * objDiscD; % v_ as a row in v_ot
    end
    
    % Sorting by abs value to minimize numerical error
    [unused I] = sort(abs(discD),2);
    for i = 1:inT
        discD(i,:) = discD(i,I(i,:));
    end
    discD = sum(discD,2);
    
    % Sum discrete and continuous terms
    if opts.UseModelICs
        % All conditions have the same variables
        D = D + contD + discD;
    else
        % Rate parameters are the same
        D(1:nTk) = D(1:nTk) + contD(1:nTk) + discD(1:nTk);
        % IC parameters are different
        D(Txind:Txind+inTx-1) = contD(nTk+1:nTk+inTx) + discD(nTk+1:nTk+inTx);
        % Increment index of variable ICs
        Txind = Txind + inTx;
    end
    
    % Return the adaptive abstol, if requested
%     if nargout >= 3
%         abstol{iCon} = opts.RelTol * abstolObjSimple(m, con(iCon), obj(:,iCon), sol);
%     end
    
    if verboseAll; fprintf('iCon = %d\t|dGdT| = %g\tTime = %0.2f\n', iCon, norm(contD + discD), toc); end
end

if opts.Verbose;fprintf('Summary: |dGdp| = %g\n', norm(D));end
