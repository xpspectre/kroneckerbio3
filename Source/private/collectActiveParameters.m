function T = collectActiveParameters(m, con, useParams, useICs, useModelICs)

% Constants
nCon = size(con, 1);
nx   = m.nX;

% Store complete parameter sets
k = m.p;

if useModelICs
    x0 = m.ic;
else
    x0 = zeros(nx, nCon);
    for iCon = 1:nCon
        x0(:,iCon) = con(iCon).ic;
    end
end

% Construct starting variable parameter set
T = [k(useParams); vec(x0(useICs))];
