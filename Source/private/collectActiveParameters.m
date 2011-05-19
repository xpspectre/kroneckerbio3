function T = collectActiveParameters(m, con, useParams, useICs, useModelICs)

% Constants
nCon = size(con, 1);
nx   = m.nx;

% Store complete parameter sets
k = m.k;

if useModelICs
    x0 = m.x0;
else
    x0 = zeros(nx, nCon);
    for iCon = 1:nCon
        x0(:,iCon) = con(iCon).x0;
    end
end

% Construct starting variable parameter set
T = [k(useParams); vec(x0(useICs))];
