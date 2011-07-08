function bounds = fixBounds(bounds, useParams, useICs, useModelICs)

% Constants
nx = size(useICs, 1);
nk = numel(useParams);
nTk = nnz(useParams);
nTx = nnz(useICs);
nT = nTk + nTx;
nCon = size(useICs, 2);

if useModelICs
    % Bounds can be (nk+nx), nT, or nk long
    l = numel(bounds);
    
    if l == nk+nx
        bounds = bounds([useParams; useICs]);
    elseif l == nT
        %bounds = bounds;
    elseif nTx == 0 && l == nk
        bounds = bounds(useParams);
    elseif l == 1
        bounds = zeros(nT,1) + bounds;
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.nx, number of varible parameters, m.nP if there are no variable ICs, or scalar')
    end
else
    % Bounds can be nk+(nx*nCon), nk+nx, nT, or nk long
    l = numel(bounds);

    if l == nk+(nx*nCon)
        bounds = bounds([useParams; vec(useICs)]);
    elseif l == nk+nx
        temp = bounds(nk+1:nk+nx); % Extract nx vector
        temp = vec(kron(temp,ones(1,nCon))); % Expand it
        temp = temp(vec(useICs)); % Extract active IC parameters
        bounds = [bounds(useParams); temp]; % Concatenate active rate parameters
    elseif l == nT
        %bounds = bounds;
    elseif nTx == 0 && l == nk
        bounds = bounds(useParams);
    elseif l == 1
        bounds = zeros(nT,1) + bounds;
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+(m.nx*length(con)), m.nk+m.nx, number of varible parameters, m.nk if there are no variable ICs, or scalar')
    end
end
