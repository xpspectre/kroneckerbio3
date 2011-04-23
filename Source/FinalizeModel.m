function m = FinalizeModel(m)
%FinalizeModel updates the mathematical components of the model to reflect
%   changes made to the model
%
%   m = FinalizeModel(m)
%
%   Inputs
%   m: [ model struct scalar ]
%
%   Outputs
%   m: [ model struct scalar ]
%
%   Some components of Kronecker models are cyclically dependent
%   (compartment volumes depend on species and species are placed into
%   compartments). Because of this, it is not possible to design Kronecker
%   in such a way that the model is ready to use after the addition of any
%   component. Furthermore, rebuilding the mathematical components of
%   Kronecker takes a non-trivial amount of time. Since most components are
%   added in groups, simply storing the components until they are all added
%   saves time. Unfortunately, this requires the user to remember to call
%   this function on the model when he is done.
%
%   Forgetting to call this function and using a model which has Ready set
%   to false will result in undefined behavior.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Place compartments
nvNew = m.add.nv;

% Check if item by this name already exists
handled = false(nvNew,1);
for iv = 1:nvNew
    % Check existing compartments
    matchPosition = find(strcmp(m.add.Compartments(iv).Name, {m.Compartments.Name}));
    if ~isempty(matchPosition)
        handled(iv) = true;
        m.Compartments(matchPosition) = m.add.Compartments(iv);
    end
    
    % Check newly added compartments
    matchPosition = find(strcmp(m.add.Compartments(iv).Name, {m.add.Compartments(1:iv-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new item
m.Compartments = [m.Compartments; m.add.Compartments(~handled)];

% Update count
m.nv = numel(m.Compartments);

%% Place species
nxuNew = m.add.nxu;

% Check if item by this name already exists
handled = false(nxuNew,1);
for ixu = 1:nxuNew
    % Check existing states
    matchPosition = find(strcmp(m.add.Species(ixu).Name, {m.Species.Name}) & strcmp(m.add.Species(ixu).Compartment, {m.Species.Compartment}));
    if ~isempty(matchPosition)
        handled(ixu) = true;
        m.Species(matchPosition) = rmfield(m.add.Species(ixu), 'Units');
    end
    
    % Check newly added states
    matchPosition = find(strcmp(m.add.Species(ixu).Name, {m.add.Species(1:ixu-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new item
m.Species = [m.Species; rmfield(m.add.Species(~handled), 'Units')];

% Update count
isu = cat(1, m.Species.IsInput, false(0,1));
m.nu = nnz(isu);
m.nx = numel(m.Species) - m.nu;
xuNamesFull = strcat({m.Species.Compartment}, '.', {m.Species.Name});
xNamesFull = strcat({m.Species(~isu).Compartment}, '.', {m.Species(~isu).Name});
uNamesFull = strcat({m.Species(isu).Compartment}, '.', {m.Species(isu).Name});

%% Place outputs
nyNew = m.add.ny;

% Check if item by this name already exists
handled = false(nyNew,1);
for iy = 1:nyNew
    % Check existing outputs
    matchPosition = find(strcmp(m.add.Outputs(iy).Name, {m.Outputs.Name}));
    if ~isempty(matchPosition)
        handled(iy) = true;
        m.Outputs(matchPosition) = m.add.Outputs(iy);
    end
    
    % Check newly added outputs
    matchPosition = find(strcmp(m.add.Outputs(iy).Name, {m.add.Outputs(1:iy-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new items
m.Outputs = [m.Outputs; m.add.Outputs(~handled)];

% Update count
m.ny = numel(m.Outputs);

%% Place parameters
nkNew = m.add.nk;

% Check if item by this name already exists
handled = false(nkNew,1);
for ik = 1:nkNew
    % Check existing parameters
    matchPosition = find(strcmp(m.add.Parameters(ik).Name, {m.Parameters.Name}));
    if ~isempty(matchPosition)
        handled(ik) = true;
        m.Parameters(matchPosition) = rmfield(m.add.Parameters(ik), 'Units');
    end
    
    % Check newly added parameters
    matchPosition = find(strcmp(m.add.Parameters(ik).Name, {m.add.Parameters(1:ik-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new items
m.Parameters = [m.Parameters; rmfield(m.add.Parameters(~handled), 'Units')];

% Update count
m.nk = numel(m.Parameters);

%% Expand reactions
nrNew = m.add.nr; % Number of reaction specifications

% Loop over added reaction specifications
rNew = cell(nrNew,1); % Containers for all the spawned reactions
rCount = zeros(nrNew,1); % Number of reactions that have spawned
for ir = 1:nrNew
    % Determine possible compartments as a cell array of strings
    if isempty(m.add.Reactions(ir).Compartment)
        possibleComp = {m.Compartments.Name};
    elseif ischar(m.add.Reactions(ir).Compartment)
        possibleComp = {m.add.Reactions(ir).Compartment};
    elseif iscell(m.add.Reactions(ir).Compartment)
        possibleComp = m.add.Reactions(ir).Compartment;
    end
    
    % Container for all the reactions that come from this specification
    nvPossible = numel(possibleComp);
    rNew{ir} = emptyReactions(nvPossible);
    
    % Match reactants and products
    if any(m.add.Reactions(ir).Reactants{1} == '.')
        reactant1Complete = true;
        possibleReactant1 = strcmp(m.add.Reactions(ir).Reactants{1}, xuNamesFull);
    else
        reactant1Complete = false;
        possibleReactant1 = strcmp(m.add.Reactions(ir).Reactants{1}, {m.Species.Name});
    end
    if any(m.add.Reactions(ir).Reactants{2} == '.')
        reactant2Complete = true;
        possibleReactant2 = strcmp(m.add.Reactions(ir).Reactants{2}, xuNamesFull);
    else
        reactant2Complete = false;
        possibleReactant2 = strcmp(m.add.Reactions(ir).Reactants{2}, {m.Species.Name});
    end
    if any(m.add.Reactions(ir).Products{1} == '.')
        product1Complete = true;
        possibleProduct1 = strcmp(m.add.Reactions(ir).Products{1}, xuNamesFull);
    else
        product1Complete = false;
        possibleProduct1 = strcmp(m.add.Reactions(ir).Products{1}, {m.Species.Name});
    end
    if any(m.add.Reactions(ir).Products{2} == '.')
        product2Complete = true;
        possibleProduct2 = strcmp(m.add.Reactions(ir).Products{2}, xuNamesFull);
    else
        product2Complete = false;
        possibleProduct2 = strcmp(m.add.Reactions(ir).Products{2}, {m.Species.Name});
    end
    
    % Loop over possible compartments
    aReactionWasFound = false; % Used to throw a warning if this remains false
    for iv = 1:nvPossible
        % Match compartments
        possibleCompartment = strcmp(possibleComp{iv}, {m.Species.Compartment});
        
        % Find reactants and products for this compartment
        reactant1Exists = ~isempty(m.add.Reactions(ir).Reactants{1});
        reactant2Exists = ~isempty(m.add.Reactions(ir).Reactants{2});
        product1Exists  = ~isempty(m.add.Reactions(ir).Products{1});
        product2Exists  = ~isempty(m.add.Reactions(ir).Products{2});
        reactant1Found  = ~reactant1Exists || (any(possibleReactant1) && reactant1Complete) || any(possibleReactant1 & possibleCompartment);
        reactant2Found  = ~reactant2Exists || (any(possibleReactant2) && reactant2Complete) || any(possibleReactant2 & possibleCompartment);
        product1Found   = ~product1Exists  || (any(possibleProduct1)  && product1Complete)  || any(possibleProduct1 & possibleCompartment);
        product2Found   = ~product2Exists  || (any(possibleProduct2)  && product2Complete)  || any(possibleProduct2 & possibleCompartment);
        parameter1Exists = ~isempty(m.add.Reactions(ir).Parameters{1});
        parameter2Exists = ~isempty(m.add.Reactions(ir).Parameters{2});
        
        % Check for errors
        if ((reactant1Exists || reactant2Exists) && (reactant1Found && reactant2Found) && (~product1Found || ~product2Found) && parameter1Exists)
            error('KroneckerBio:FinalizeModel:MissingProduct', 'Reaction %s (#%i) has the necessary reactants in compartment %s but not the necessary products', m.add.Reactions(ir).Names{1}, ir, possibleComp{iv})
        end
        if ((product1Exists || product2Exists) && (product1Found && product2Found) && (~reactant1Found || ~reactant2Found) && parameter2Exists)
            error('KroneckerBio:FinalizeModel:MissingProduct', 'Reverse reaction %s (#%i) has the necessary reactants in compartment %s but not the necessary products', m.add.Reactions(ir).Names{2}, ir, possibleComp{iv})
        end
        
        % Check if this reaction takes place
        if ~(reactant1Found && reactant2Found && product1Found && product2Found)
            continue
        end
        
        % Forward reaction
        if parameter1Exists
            % Forward reaction exists in this compartment
            rCount(ir) = rCount(ir) + 1;
            
            % Add reaction to list
            rNew{ir}(rCount(ir)).Name         = m.add.Reactions(ir).Names{1};
            rNew{ir}(rCount(ir)).Reactants{1} = fixReactionSpecies(m.add.Reactions(ir).Reactants{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Reactants{2} = fixReactionSpecies(m.add.Reactions(ir).Reactants{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{1}  = fixReactionSpecies(m.add.Reactions(ir).Products{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{2}  = fixReactionSpecies(m.add.Reactions(ir).Products{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Parameter    = m.add.Reactions(ir).Parameters{1};
            
            aReactionWasFound = true;
        end
        
        % Reverse reaction
        if parameter2Exists
            % Forward reaction exists in this compartment
            rCount(ir) = rCount(ir) + 1;
            
            % Add reaction to list
            rNew{ir}(rCount(ir)).Name         = m.add.Reactions(ir).Names{2};
            rNew{ir}(rCount(ir)).Reactants{1} = fixReactionSpecies(m.add.Reactions(ir).Products{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Reactants{2} = fixReactionSpecies(m.add.Reactions(ir).Products{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{1}  = fixReactionSpecies(m.add.Reactions(ir).Reactants{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{2}  = fixReactionSpecies(m.add.Reactions(ir).Reactants{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Parameter    = m.add.Reactions(ir).Parameters{2};

            aReactionWasFound = true;
        end
    end
    
    % Warn if no compartment had the species for this reaction specification
    if ~aReactionWasFound
        warning('KroneckerBio:FinalizeModel:UnusedReaction', 'Reaction %s (#%i) was not found to take place in any compartment', m.add.Reactions(ir).Names{1}, ir)
    end
end

%% Place reactions
% Total number of new elementary reactions
nrTotal = sum(rCount);

% Add room for new reactions
m.Reactions = [m.Reactions; emptyReactions(nrTotal)];

% Insert items
rEndIndex = m.nr;
for ir = 1:nrNew
    rStartIndex = rEndIndex + 1;
    rEndIndex   = rEndIndex + rCount(ir);
    m.Reactions(rStartIndex:rEndIndex) = rNew{ir}(1:rCount(ir));
end

% Update count
m.nr = numel(m.Reactions);

%% Useful information
% Constants
nv = m.nv;
nx = m.nx;
nu = m.nu;
nxu = nx + nu;
ny = m.ny;
nk = m.nk;
nr = m.nr;

% State compartments
m.vxInd = zeros(nx,1);
xInd = find(~isu);
for ix = 1:nx
    m.vxInd(ix) = find(strcmp(m.Species(xInd(ix)).Compartment, {m.Compartments.Name}));
end

% Input compartments
m.vuInd = zeros(nu,1);
uInd = find(isu);
for iu = 1:nu
    m.vuInd(iu) = find(strcmp(m.Species(uInd(iu)).Compartment, {m.Compartments.Name}));
end

% Map from species to x0 and u
x0uInd = zeros(nxu,1);
x0uInd(~isu) = 1:nx;
x0uInd(isu) = 1:nu;

%% Process compartments
% Dimensions
m.d = [m.Compartments.Dimension].';

% Entries in each sparse matrix for compartment conversion
nB1Entries = 0;
nB2Entries = 0;
nbEntries = 0;

B1Entries = zeros(0,2);
B1Values  = zeros(0,1);
B2Entries = zeros(0,2);
B2Values  = zeros(0,1);
bEntries  = zeros(0,2);
bValues   = zeros(0,1);

for iv = 1:nv
    nExpr = numel(m.Compartments(iv).Expressions);
    for iExpr = 1:nExpr
        % Find states that match the expression
        match = find(~cellfun(@isempty, regexp(xNamesFull, m.Compartments(iv).Expressions{iExpr}, 'once')));
        nAdd = numel(match);
        nB1Entries = nB1Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(B1Entries,1);
        if nB1Entries > currentLength
            addlength = max(currentLength, nAdd);
            B1Entries = [B1Entries; zeros(addlength,2)];
        end
        
        % Add entries
        B1Entries(nB1Entries-nAdd+1:nB1Entries,1) = iv;
        B1Entries(nB1Entries-nAdd+1:nB1Entries,2) = match;
        B1Values(nB1Entries-nAdd+1:nB1Entries) = m.Compartments(iv).Values(iExpr);

        % Find inputs that match the expression
        match = find(~cellfun(@isempty, regexp(uNamesFull, m.Compartments(iv).Expressions{iExpr}, 'once')));
        nAdd = numel(match);
        nB2Entries = nB2Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(B2Entries,1);
        if nB2Entries > currentLength
            addlength = max(currentLength, nAdd);
            B2Entries = [B2Entries; zeros(addlength,2)];
        end
        
        % Add entries
        B2Entries(nB2Entries-nAdd+1:nB2Entries,1) = iv;
        B2Entries(nB2Entries-nAdd+1:nB2Entries,2) = match;
        B2Values(nB2Entries-nAdd+1:nB2Entries) = m.Compartments(iv).Values(iExpr);

        % Find empty expressions, which are constants
        if isempty(m.Compartments(iv).Expressions{iExpr})
            nbEntries = nbEntries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(bEntries,1);
            if nbEntries > currentLength
                addlength = max(currentLength, 1);
                bEntries = [bEntries; zeros(addlength,2)];
            end
            
            % Add entries
            bEntries(nbEntries,1) = iv;
            bEntries(nbEntries,2) = 1;
            bValues(nbEntries) = m.Compartments(iv).Values(iExpr);
        end
    end
end

% Remove duplicate entries
[B1Entries ind] = unique(B1Entries(1:nB1Entries,:), 'rows');
B1Values = B1Values(ind);

[B2Entries ind] = unique(B2Entries(1:nB2Entries,:), 'rows');
B2Values = B2Values(ind);

[bEntries ind] = unique(bEntries(1:nbEntries,:), 'rows');
bValues = bValues(ind);

% Construct matrices
m.B1 = sparse(B1Entries(:,1), B1Entries(:,2), B1Values, m.nv, m.nx);
m.B2 = sparse(B2Entries(:,1), B2Entries(:,2), B2Values, m.nv, m.nu);
m.b  = sparse(bEntries(:,1),  bEntries(:,2),  bValues,  m.nv, 1);

%% Compute non-concentration species values
x0Handled = true(nx,1);
uHandled = true(nu,1);
xuNewInd = zeros(nxuNew,1);
for ixuNew = 1:nxuNew
    % Index in Species.Names that this new Species applies
    xuIndi = find(strcmp(m.add.Species(ixuNew).Name, {m.Species.Name}));
    xuNewInd(ixuNew) = xuIndi;
    
    if m.add.Species(ixuNew).IsInput
        [m.Species(xuIndi).Value.Function uHandled(x0uInd(xuIndi))] = fixInputValueFunction(m.add.Species(ixuNew).Value.Function, m.add.Species(ixuNew).Units);
    else
        [m.Species(xuIndi).Value x0Handled(x0uInd(xuIndi))] = fixStateInitialValue(m.add.Species(ixuNew).Value, m.add.Species(ixuNew).Units);
    end
end

%% Compute starting compartment volumes
% Assert that no unhandled species play a part in the compartment volume
if nnz(m.B1(:,~x0Handled)) > 0
    errorInd = find(m.B1(:,~x0Handled));
    errorInd = ind2sub(errorInd, m.nv, nnz(~x0Handled));
    error('KroneckerBio:FinalizeModel:ConcentrationAffectsVolume', 'Initial condition for species %s was given in concentration, but its value affects compartment %s; the amount in the concentration cannot be resolved', m.States(errorInd(2)).Name, m.Compartments(errorInd(1)).Name)
end
if nnz(m.B2(:,~uHandled)) > 0
    errorInd = find(m.B1(:,~x0Handled));
    errorInd = ind2sub(errorInd, m.nv, nnz(~x0Handled));
    error('KroneckerBio:FinalizeModel:ConcentrationAffectsVolume', 'Input %s was given in concentration, but its value affects a compartment; the amount in the concentration cannot be resolved', m.Inputs(errorInd(2)).Name, m.Compartments(errorInd(1)).Name)
end

x0 = cat(1, m.Species(~isu).Value, zeros(0,1));
u = completeInputFunction(m.Species(isu));
v0 = m.B1 * x0 + m.B2 * u(0) + m.b;

%% Compute remaining species values
for ixuNew = 1:nxuNew
    % Index in m.Species matching this new Species
    xuIndi = xuNewInd(ixuNew);
    
    % Only do something if it was not already handled
    if m.Species(xuIndi).IsInput && ~uHandled(x0uInd(xuIndi))
        % Index in m.Compartments
        vInd = strcmp(m.Species(xuIndi).Compartment, {m.Compartments.Name});
        
        m.Species(xuIndi).Value.Function = fixInputValueFunction(m.add.Species(ixuNew).Value.Function, m.add.Species(ixuNew).Units, v0(vInd));
    elseif ~m.Species(xuIndi).IsInput && ~x0Handled(x0uInd(xuIndi))
        % Index in m.Compartments
        vInd = strcmp(m.Species(xuIndi).Compartment, {m.Compartments.Name});
        
        m.Species(xuIndi).Value = fixStateInitialValue(m.add.Species(ixuNew).Value, m.add.Species(ixuNew).Units, v0(vInd));
    end
end

%% Process states
% Put initial conditions into x0 vector
m.x0 = cat(1, m.Species(~isu).Value, zeros(0,1));

%% Process inputs
% Condense time varying inputs into u(t) function
m.u = completeInputFunction(m.Species(isu));

% Put input parameters into q vector
m.nqu = zeros(nu,1);
uInd = find(isu);
for iu = 1:nu
    m.nqu(iu) = numel(m.Species(uInd(iu)).Value.Parameters);
end
uValue = cat(1, m.Species(isu).Value, struct('Function', cell(0,1), 'Parameters', cell(0,1)));
m.q = cat(1, uValue.Parameters, zeros(0,1));
nq = numel(m.q);
m.nq = nq;

%% Process outputs
% Entries in each sparse matrix for compartment conversion
nC1Entries = 0;
nC2Entries = 0;
ncEntries = 0;

C1Entries = zeros(0,2);
C1Values  = zeros(0,1);
C2Entries = zeros(0,2);
C2Values  = zeros(0,1);
cEntries  = zeros(0,2);
cValues   = zeros(0,1);

for iy = 1:ny
    nExpr = numel(m.Outputs(iy).Expressions);
    for iExpr = 1:nExpr
        % Find states that match the expression
        match = find(~cellfun(@isempty, regexp(xNamesFull, m.Outputs(iy).Expressions{iExpr}, 'once')));
        nAdd = numel(match);
        nC1Entries = nC1Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(C1Entries,1);
        if nC1Entries > currentLength
            addlength = max(currentLength, nAdd);
            C1Entries = [C1Entries; zeros(addlength,2)];
            C1Values  = [C1Values;  zeros(addlength,1)];
        end
        
        % Add entries
        C1Entries(nC1Entries-nAdd+1:nC1Entries,1) = iy;
        C1Entries(nC1Entries-nAdd+1:nC1Entries,2) = match;
        C1Values(nC1Entries-nAdd+1:nC1Entries) = m.Outputs(iy).Values(iExpr);

        % Find inputs that match the expression
        match = find(~cellfun(@isempty, regexp(uNamesFull, m.Outputs(iy).Expressions{iExpr}, 'once')));
        nAdd = numel(match);
        nC2Entries = nC2Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(C2Entries,1);
        if nC2Entries > currentLength
            addlength = max(currentLength, nAdd);
            C2Entries = [C2Entries; zeros(addlength,2)];
            C2Values  = [C2Values; zeros(addlength,1)];
        end
        
        % Add entries
        C2Entries(nC2Entries-nAdd+1:nC2Entries,1) = iy;
        C2Entries(nC2Entries-nAdd+1:nC2Entries,2) = match;
        C2Values(nC2Entries-nAdd+1:nC2Entries) = m.Outputs(iy).Values(iExpr);

        % Find empty expressions, which are constants
        if isempty(m.Outputs(iy).Expressions{iExpr})
            ncEntries = ncEntries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(cEntries,1);
            if nbEntries > currentLength
                addlength = max(currentLength, 1);
                cEntries = [cEntries; zeros(addlength,2)];
                cValues = [cValues; zeros(addlength,1)];
            end
            
            % Add entries
            cEntries(ncEntries,1) = iy;
            cEntries(ncEntries,2) = 1;
            cValues(ncEntries) = m.Outputs(iy).Values(iExpr);
        end
    end
end

% Remove duplicate entries
[C1Entries ind] = unique(C1Entries(1:nC1Entries,:), 'rows');
C1Values = C1Values(ind);

[C2Entries ind] = unique(C2Entries(1:nC2Entries,:), 'rows');
C2Values = C2Values(ind);

[cEntries ind] = unique(cEntries(1:ncEntries,:), 'rows');
cValues = cValues(ind);

% Construct matrices
m.C1 = sparse(C1Entries(:,1), C1Entries(:,2), C1Values, m.ny, m.nx);
m.C2 = sparse(C2Entries(:,1), C2Entries(:,2), C2Values, m.ny, m.nu);
m.c  = sparse(cEntries(:,1),  cEntries(:,2),  cValues,  m.ny, 1);

%% Process parameters
% Put rate parameter into k vector
m.k = cat(1, m.Parameters.Value, zeros(0,1));

%% Process Reactions
nA1Entries = 0;
nA2Entries = 0;
nA3Entries = 0;
nA4Entries = 0;
nA5Entries = 0;
nA6Entries = 0;
naEntries = 0;

dA1dkEntries = zeros(0,2);
dA1dkValues  = zeros(0,1);
dA2dkEntries = zeros(0,2);
dA2dkValues  = zeros(0,1);
dA3dkEntries = zeros(0,2);
dA3dkValues  = zeros(0,1);
dA4dkEntries = zeros(0,2);
dA4dkValues  = zeros(0,1);
dA5dkEntries = zeros(0,2);
dA5dkValues  = zeros(0,1);
dA6dkEntries = zeros(0,2);
dA6dkValues  = zeros(0,1);
dadkEntries = zeros(0,2);
dadkValues  = zeros(0,1);

for ir =1:nr
    % Determine reaction type and species indexes
    if isempty(m.Reactions(ir).Reactants{1})
        reactant1Exists = false;
    else
        reactant1Exists = true;
        reactant1IsState = true;
        reactant1 = find(strcmp(m.Reactions(ir).Reactants{1}, xNamesFull));
        if isempty(reactant1)
            reactant1IsState = false;
            reactant1 = find(strcmp(m.Reactions(ir).Reactants{1}, uNamesFull));
        end
    end
    
    if isempty(m.Reactions(ir).Reactants{2})
        reactant2Exists = false;
    else
        reactant2Exists = true;
        reactant2IsState = true;
        reactant2 = find(strcmp(m.Reactions(ir).Reactants{2}, xNamesFull));
        if isempty(reactant2)
            reactant2IsState = false;
            reactant2 = find(strcmp(m.Reactions(ir).Reactants{2}, uNamesFull));
        end
    end
    
    product1 = find(strcmp(m.Reactions(ir).Products{1}, xNamesFull));
    product1IsState = ~isempty(product1);
    
    product2 = find(strcmp(m.Reactions(ir).Products{2}, xNamesFull));
    product2IsState = ~isempty(product2);
    
    parameter = find(strcmp(m.Reactions(ir).Parameter, {m.Parameters.Name}));
    assert(~isempty(parameter), 'KroneckerBio:FinalizeModel:MissingReactionParameter', 'Reaction %s (#%i) requires parameter %s, but no parameter by that name was found', m.Reactions(ir).Name, ir, m.Reactions(ir).Parameter)
    
    % Switch on reactant state
    if reactant1Exists && ~reactant2Exists && reactant1IsState
        % A1 reaction
        nAdd = 1 + product1IsState + product2IsState;
        nA1Entries = nA1Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(dA1dkEntries,1);
        if nA1Entries > currentLength
            addlength = max(currentLength, nAdd);
            dA1dkEntries = [dA1dkEntries; zeros(addlength,2)];
            dA1dkValues  = [dA1dkValues;  zeros(addlength,1)];
        end
        
        % Add entries
        % Substract reactant
        dA1dkEntries(nA1Entries-nAdd+1,1) = sub2ind([nx,nx], reactant1, reactant1);
        dA1dkEntries(nA1Entries-nAdd+1,2) = parameter;
        dA1dkValues(nA1Entries-nAdd+1)    = -1;
        
        % Add product 1
        if product1IsState
            dA1dkEntries(nA1Entries-nAdd+2,1) = sub2ind([nx,nx], product1, reactant1);
            dA1dkEntries(nA1Entries-nAdd+2,2) = parameter;
            dA1dkValues(nA1Entries-nAdd+2)    = 1;
        end
        
        % Add product 2
        if product2IsState
            dA1dkEntries(nA1Entries-nAdd+product1IsState+2,1) = sub2ind([nx,nx], product2, reactant1);
            dA1dkEntries(nA1Entries-nAdd+product1IsState+2,2) = parameter;
            dA1dkValues(nA1Entries-nAdd+product1IsState+2)    = 1;
        end
    elseif reactant1Exists && reactant2Exists && reactant1IsState && reactant2IsState
        % A2 reaction
        % Order reactants so that freest species is second
        if m.d(m.vxInd(reactant1)) > m.d(m.vxInd(reactant2)) || (m.d(m.vxInd(reactant1)) == m.d(m.vxInd(reactant2)) && reactant1 > reactant2)
            [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
        end
        
        % Increment entries
        nAdd = 2 + product1IsState + product2IsState;
        nA2Entries = nA2Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(dA2dkEntries,1);
        if nA2Entries > currentLength
            addlength = max(currentLength, nAdd);
            dA2dkEntries = [dA2dkEntries; zeros(addlength,2)];
            dA2dkValues  = [dA2dkValues;  zeros(addlength,1)];
        end
        
        % Add entries
        % Substract reactant 1
        dA2dkEntries(nA2Entries-nAdd+1,1) = sub2ind([nx,nx,nx], reactant1, reactant2, reactant1);
        dA2dkEntries(nA2Entries-nAdd+1,2) = parameter;
        dA2dkValues(nA2Entries-nAdd+1)    = -1;
        
        % Substract reactant 2
        dA2dkEntries(nA2Entries-nAdd+2,1) = sub2ind([nx,nx,nx], reactant2, reactant2, reactant1);
        dA2dkEntries(nA2Entries-nAdd+2,2) = parameter;
        dA2dkValues(nA2Entries-nAdd+2)    = -1;

        % Add product 1
        if product1IsState
            dA2dkEntries(nA2Entries-nAdd+3,1) = sub2ind([nx,nx,nx], product1, reactant2, reactant1);
            dA2dkEntries(nA2Entries-nAdd+3,2) = parameter;
            dA2dkValues(nA2Entries-nAdd+3)    = 1;
        end
        
        % Add product 2
        if product2IsState
            dA2dkEntries(nA2Entries-nAdd+product1IsState+3,1) = sub2ind([nx,nx,nx], product2, reactant2, reactant1);
            dA2dkEntries(nA2Entries-nAdd+product1IsState+3,2) = parameter;
            dA2dkValues(nA2Entries-nAdd+product1IsState+3)    = 1;
        end
    elseif reactant1Exists && reactant2Exists && xor(reactant2IsState, reactant1IsState)
        % Order reactants so that freest species is second
        if reactant1IsState
            if m.d(m.vxInd(reactant1)) > m.d(m.vuInd(reactant2))
                [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
                [reactant1IsState reactant2IsState] = deal(reactant2IsState, reactant1IsState);
            end
        else
            if m.d(m.vuInd(reactant1)) >= m.d(m.vxInd(reactant2))
                [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
                [reactant1IsState reactant2IsState] = deal(reactant2IsState, reactant1IsState);
            end
        end
        
        if reactant2IsState
            % A3 reaction
            % Increment entries
            nAdd = 1 + product1IsState + product2IsState;
            nA3Entries = nA3Entries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(dA3dkEntries,1);
            if nA3Entries > currentLength
                addlength = max(currentLength, nAdd);
                dA3dkEntries = [dA3dkEntries; zeros(addlength,2)];
                dA3dkValues  = [dA3dkValues;  zeros(addlength,1)];
            end
            
            % Add entries
            % Substract reactant 2
            dA3dkEntries(nA3Entries-nAdd+1,1) = sub2ind([nx,nx,nu], reactant2, reactant2, reactant1);
            dA3dkEntries(nA3Entries-nAdd+1,2) = parameter;
            dA3dkValues(nA3Entries-nAdd+1)    = -1;
            
            % Add product 1
            if product1IsState
                dA3dkEntries(nA3Entries-nAdd+2,1) = sub2ind([nx,nx,nu], product1, reactant2, reactant1);
                dA3dkEntries(nA3Entries-nAdd+2,2) = parameter;
                dA3dkValues(nA3Entries-nAdd+2)    = 1;
            end
            
            % Add product 2
            if product2IsState
                dA3dkEntries(nA3Entries-nAdd+product1IsState+2,1) = sub2ind([nx,nx,nu], product2, reactant2, reactant1);
                dA3dkEntries(nA3Entries-nAdd+product1IsState+2,2) = parameter;
                dA3dkValues(nA3Entries-nAdd+product1IsState+2)    = 1;
            end
        else%reactant1IsState
            % A4 reaction
            % Increment entries
            nAdd = 1 + product1IsState + product2IsState;
            nA4Entries = nA4Entries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(dA4dkEntries,1);
            if nA4Entries > currentLength
                addlength = max(currentLength, nAdd);
                dA4dkEntries = [dA4dkEntries; zeros(addlength,2)];
                dA4dkValues  = [dA4dkValues;  zeros(addlength,1)];
            end
            
            % Add entries
            % Substract reactant 1
            dA4dkEntries(nA4Entries-nAdd+1,1) = sub2ind([nx,nu,nx], reactant1, reactant2, reactant1);
            dA4dkEntries(nA4Entries-nAdd+1,2) = parameter;
            dA4dkValues(nA4Entries-nAdd+1)    = -1;
            
            % Add product 1
            if product1IsState
                dA4dkEntries(nA4Entries-nAdd+2,1) = sub2ind([nx,nu,nx], product1, reactant2, reactant1);
                dA4dkEntries(nA4Entries-nAdd+2,2) = parameter;
                dA4dkValues(nA4Entries-nAdd+2)    = 1;
            end
            
            % Add product 2
            if product2IsState
                dA4dkEntries(nA4Entries-nAdd+product1IsState+2,1) = sub2ind([nx,nu,nx], product2, reactant2, reactant1);
                dA4dkEntries(nA4Entries-nAdd+product1IsState+2,2) = parameter;
                dA4dkValues(nA4Entries-nAdd+product1IsState+2)    = 1;
            end
        end
    elseif reactant1Exists && reactant2Exists && ~reactant1IsState && ~reactant2IsState
        % A5 reaction
        % Order reactants so that freest species is second
        if m.d(m.vuInd(reactant1)) > m.d(m.vuInd(reactant2)) || (m.d(m.vuInd(reactant1)) == m.d(m.vuInd(reactant2)) && reactant1 > reactant2)
            [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
        end
        
        % Increment entries
        nAdd = product1IsState + product2IsState;
        nA5Entries = nA5Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(dA5dkEntries,1);
        if nA5Entries > currentLength
            addlength = max(currentLength, nAdd);
            dA5dkEntries = [dA5dkEntries; zeros(addlength,2)];
            dA5dkValues  = [dA5dkValues;  zeros(addlength,1)];
        end
        
        % Add entries
        % Add product 1
        if product1IsState
            dA5dkEntries(nA5Entries-nAdd+1,1) = sub2ind([nx,nu,nu], product1, reactant2, reactant1);
            dA5dkEntries(nA5Entries-nAdd+1,2) = parameter;
            dA5dkValues(nA5Entries-nAdd+1)    = 1;
        end
        
        % Add product 2
        if product2IsState
            dA5dkEntries(nA5Entries-nAdd+product1IsState+1,1) = sub2ind([nx,nu,nu], product2, reactant2, reactant1);
            dA5dkEntries(nA5Entries-nAdd+product1IsState+1,2) = parameter;
            dA5dkValues(nA5Entries-nAdd+product1IsState+1)    = 1;
        end
    elseif reactant1Exists && ~reactant2Exists && ~reactant1IsState
        % A6 reaction
        nAdd = product1IsState + product2IsState;
        nA6Entries = nA6Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(dA6dkEntries,1);
        if nA6Entries > currentLength
            addlength = max(currentLength, nAdd);
            dA6dkEntries = [dA6dkEntries; zeros(addlength,2)];
            dA6dkValues  = [dA6dkValues;  zeros(addlength,1)];
        end
        
        % Add entries
        % Add product 1
        if product1IsState
            dA6dkEntries(nA6Entries-nAdd+1,1) = sub2ind([nx,nu], product1, reactant1);
            dA6dkEntries(nA6Entries-nAdd+1,2) = parameter;
            dA6dkValues(nA6Entries-nAdd+1)  = 1;
        end
        
        % Add product 2
        if product2IsState
            dA6dkEntries(nA6Entries-nAdd+product1IsState+1,1) = sub2ind([nx,nu], product2, reactant1);
            dA6dkEntries(nA6Entries-nAdd+product1IsState+1,2) = parameter;
            dA6dkValues(nA6Entries-nAdd+product1IsState+1)    = 1;
        end
    elseif ~reactant1Exists && ~reactant2Exists
        % a reaction
        nAdd = product1IsState + product2IsState;
        naEntries = naEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(dadkEntries,1);
        if naEntries > currentLength
            addlength = max(currentLength, nAdd);
            dadkEntries = [dadkEntries; zeros(addlength,2)];
            dadkValues  = [dadkValues;  zeros(addlength,1)];
        end
        
        % Add entries
        % Add product 1
        if product1IsState
            dadkEntries(naEntries-nAdd+1,1) = product1;
            dadkEntries(naEntries-nAdd+1,2) = parameter;
            dadkValues(naEntries-nAdd+1)  = 1;
        end
        
        % Add product 2
        if product2IsState
            dadkEntries(naEntries-nAdd+product1IsState+1,1) = product2;
            dadkEntries(naEntries-nAdd+product1IsState+1,2) = parameter;
            dadkValues(naEntries-nAdd+product1IsState+1)    = 1;
        end
    else
        error('Unknown error occured while processing reactions')
    end
end

% Construct derivative matrices
m.dA1dk = sparse(dA1dkEntries(1:nA1Entries,1), dA1dkEntries(1:nA1Entries,2), dA1dkValues(1:nA1Entries), nx*nx,    nk);
m.dA2dk = sparse(dA2dkEntries(1:nA2Entries,1), dA2dkEntries(1:nA2Entries,2), dA2dkValues(1:nA2Entries), nx*nx*nx, nk);
m.dA3dk = sparse(dA3dkEntries(1:nA3Entries,1), dA3dkEntries(1:nA3Entries,2), dA3dkValues(1:nA3Entries), nx*nu*nx, nk);
m.dA4dk = sparse(dA4dkEntries(1:nA4Entries,1), dA4dkEntries(1:nA4Entries,2), dA4dkValues(1:nA4Entries), nx*nx*nu, nk);
m.dA5dk = sparse(dA5dkEntries(1:nA5Entries,1), dA5dkEntries(1:nA5Entries,2), dA5dkValues(1:nA5Entries), nx*nu*nu, nk);
m.dA6dk = sparse(dA6dkEntries(1:nA6Entries,1), dA6dkEntries(1:nA6Entries,2), dA6dkValues(1:nA6Entries), nx*nu,    nk);
m.dadk  = sparse(dadkEntries(1:naEntries,1),   dadkEntries(1:naEntries,2),   dadkValues(1:naEntries),   nx,       nk);

% Alternative shape
m.dA1dk_fk_x  = spermute132(m.dA1dk, [nx,nx,nk],    [nx*nk,nx]);
m.dA2dk_fk_xx = spermute132(m.dA2dk, [nx,nx*nx,nk], [nx*nk,nx*nx]);
m.dA3dk_fk_ux = spermute132(m.dA3dk, [nx,nu*nx,nk], [nx*nk,nu*nx]);
m.dA4dk_fk_xu = spermute132(m.dA4dk, [nx,nx*nu,nk], [nx*nk,nx*nu]);
m.dA5dk_fk_uu = spermute132(m.dA5dk, [nx,nu*nu,nk], [nx*nk,nu*nu]);
m.dA6dk_fk_u  = spermute132(m.dA6dk, [nx,nu,nk],    [nx*nk,nu]);

% Construct fullest possible matrices
m.A1 = reshape(m.dA1dk * rand(nk,1), nx,nx);
m.A2 = reshape(m.dA2dk * rand(nk,1), nx,nx*nx);
m.A3 = reshape(m.dA3dk * rand(nk,1), nx,nu*nx);
m.A4 = reshape(m.dA4dk * rand(nk,1), nx,nx*nu);
m.A5 = reshape(m.dA5dk * rand(nk,1), nx,nu*nu);
m.A6 = reshape(m.dA6dk * rand(nk,1), nx,nu);
m.a  = m.dadk  * rand(nk,1);

% Determine used columns of bimolecular matrices
[unused A2UsedColumns] = find(m.A2);
[unused A3UsedColumns] = find(m.A3);
[unused A4UsedColumns] = find(m.A4);
[unused A5UsedColumns] = find(m.A5);

A2UsedColumns = unique(A2UsedColumns);
A3UsedColumns = unique(A3UsedColumns);
A4UsedColumns = unique(A4UsedColumns);
A5UsedColumns = unique(A5UsedColumns);

[A2UsedSpecies2 A2UsedSpecies1] = ind2sub([nx,nx], A2UsedColumns);
[A3UsedSpecies2 A3UsedSpecies1] = ind2sub([nx,nu], A3UsedColumns);
[A4UsedSpecies2 A4UsedSpecies1] = ind2sub([nu,nx], A4UsedColumns);
[A5UsedSpecies2 A5UsedSpecies1] = ind2sub([nu,nu], A5UsedColumns);

%% Empty out added items
m.add.nv = 0;
m.add.nxu = 0;
m.add.ny = 0;
m.add.nk = 0;
m.add.nr = 0;
m.add.Compartments = growCompartments([],0);
m.add.Species      = growSpecies([],0);
m.add.Outputs      = growOutputs([],0);
m.add.Parameters   = growParameters([],0);
m.add.Reactions    = growReactions([],0);

%% Final build of model
m = final(m, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2);

end

function m = final(m, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2)

% Constants
nx = m.nx;
nu = m.nu;
nk = m.nk;
isu = cat(1, m.Species.IsInput, false(0,1));

% Build kronecker matrices
sparsek = sparse(m.k);
m.A1 = reshape(m.dA1dk * sparsek, nx,nx);
m.A2 = reshape(m.dA2dk * sparsek, nx,nx*nx);
m.A3 = reshape(m.dA3dk * sparsek, nx,nu*nx);
m.A4 = reshape(m.dA4dk * sparsek, nx,nx*nu);
m.A5 = reshape(m.dA5dk * sparsek, nx,nu*nu);
m.A6 = reshape(m.dA6dk * sparsek, nx,nu);
m.a  = m.dadk * m.k;

% Handles
m.drdx = @drdx;
m.drdu = @drdu;
m.drdk = @drdk;

m.d2rdx2  = @d2rdx2;
m.d2rdk2  = @d2rdk2;
m.d2rdxdk = @d2rdxdk;
m.d2rdkdx = @d2rdkdx;

m.f = fHidden(m.A1, m.A2, m.A3, m.A4, m.A5, m.A6, m.a, m.B1, m.B2, m.b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2, m.vxInd, m.vuInd);

m.dfdx = dfdxHidden(m.A1, m.A2, m.A3, m.A4, m.B1, m.B2, m.b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, m.vxInd, m.vuInd);
m.dfdu = dfduHidden(m.A3, m.A4, m.A5, m.A6, m.B1, m.B2, m.b, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2, m.vxInd, m.vuInd);
m.dfdk = dfdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.dadk, m.B1, m.B2, m.b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2, m.vxInd, m.vuInd);

m.d2fdx2  = d2fdx2Hidden(m.A2, m.B1, m.B2, m.b, A2UsedColumns, A2UsedSpecies2, m.vxInd);
m.d2fdk2  = d2fdk2Hidden(m.nx, m.nk);
m.d2fdxdk = d2fdxdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.B1, m.B2, m.b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdx = d2fdkdxHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.B1, m.B2, m.b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, m.vxInd, m.vuInd);

m.Ready = true;
m.Update = @Update;

    function mout = Update(k, x0, q)
        % Copy existing model
        mout = m;
        
        % Apply changes
        mout.k = k;
        mout.x0 = x0;
        mout.q = q;
        
        % Distribute values
        k = num2cell(k);
        [mout.Parameters.Value] = k{:};
        
        x0 = num2cell(x0);
        [mout.Species(~isu).InitialValue] = x0{:};
        
        qIndex = 0;
        uInd = find(isu);
        for iu = 1:nu
            mout.Species(uInd(iu)).Value.Parameters = q(qIndex+1:qIndex+m.nqu(iu));
            qIndex = qIndex + m.nqu(iu);
        end
        
        % Rebuild model
        mout = final(mout, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2);
    end
end

function handle = fHidden(A1, A2, A3, A4, A5, A6, a, B1, B2, b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @f;

    function val = f(t, x, u)
%   f = A1*x + A2*(x kron x/vx) + A3*(u kron x/vx) + A4*(x kron u/vu) + A5*(u kron u/vu) + A6*u + a
        % Compartment column
        v = B1 * x + B2 * u + b;
        xvx = x ./ v(vxInd);
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(A2UsedColumns, ones(numel(A2UsedColumns),1), x(A2UsedSpecies1) .* xvx(A2UsedSpecies2), nx*nx,1);
        ukronx = sparse(A3UsedColumns, ones(numel(A3UsedColumns),1), u(A3UsedSpecies1) .* xvx(A3UsedSpecies2), nx*nu,1);
        xkronu = sparse(A4UsedColumns, ones(numel(A4UsedColumns),1), x(A4UsedSpecies1) .* uvu(A4UsedSpecies2), nu*nx,1);
        ukronu = sparse(A5UsedColumns, ones(numel(A5UsedColumns),1), u(A5UsedSpecies1) .* uvu(A5UsedSpecies2), nu*nu,1);
        
        val = A1 * x + A2 * xkronx + A3 * ukronx + A4 * xkronu + A5 * ukronu + A6 * u + a; % f_
    end
end

function handle = dfdxHidden(A1, A2, A3, A4, B1, B2, b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdx;

    function val = dfdx(t,x,u)
%   dfdx = A1 + A2*(Ix kron x/vx) + A2*(x kron diag(1/vx)) + A3*(u kron diag(1/vx))
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(A2UsedColumns, A2UsedSpecies1, xvx(A2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(A2UsedColumns, A2UsedSpecies2, x(A2UsedSpecies1) .* vx(A2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(A3UsedColumns, A3UsedSpecies2, u(A3UsedSpecies1) .* vx(A3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(A4UsedColumns, A4UsedSpecies1, uvu(A4UsedSpecies2), nx*nu,nx);
        
        val = A1 + A2 * (Ixkronxvx + xkron1vx) + A3 * ukron1vx + A4 * Ixkronuvu; % f_x
    end
end

function handle = dfduHidden(A3, A4, A5, A6, B1, B2, b, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdu;

    function val = dfdu(t, x, u)
%   dfdu = A3*(Iu kron x/vx) + A4*(x kron diag(1/vu)) + A5*(Iu kron u/vu) + A5*(u kron diag(1/vu)) + A6
        % Compartment column
        v = B1 * x + B2 * u + b;
        xvx = x ./ v(vxInd);
        vu = 1 ./ v(vuInd);
        uvu = u .* vu;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(A3UsedColumns, A3UsedSpecies1, xvx(A3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(A4UsedColumns, A4UsedSpecies2, x(A4UsedSpecies1) .* vu(A4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(A5UsedColumns, A5UsedSpecies1, uvu(A5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(A5UsedColumns, A5UsedSpecies2, u(A5UsedSpecies1) .* vu(A5UsedSpecies2), nu*nu,nu);
        
        val = A3 * Iukronxvx + A4 * xkron1vu + A5 * (Iukronuvu + ukron1vu) + A6; % f_u
    end
end

function handle = dfdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, dadk, B1, B2, b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, A5UsedColumns, A5UsedSpecies1, A5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @dfdk;

    function val = dfdk(t, x, u)
%   dfdk = dA1dk*x + dA2dk*(x kron x/vx) + dA3dk*(u kron x/vx) + dA4dk*(x kron u/vu) + dA5dk*(u kron u/vu) + dA6dk*u + dadk
        % Compartment column
        v = B1 * x + B2 * u + b;
        xvx = x ./ v(vxInd);
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(A2UsedColumns, ones(numel(A2UsedColumns),1), x(A2UsedSpecies1) .* xvx(A2UsedSpecies2), nx*nx,1);
        ukronx = sparse(A3UsedColumns, ones(numel(A3UsedColumns),1), u(A3UsedSpecies1) .* xvx(A3UsedSpecies2), nx*nu,1);
        xkronu = sparse(A4UsedColumns, ones(numel(A4UsedColumns),1), x(A4UsedSpecies1) .* uvu(A4UsedSpecies2), nu*nx,1);
        ukronu = sparse(A5UsedColumns, ones(numel(A5UsedColumns),1), u(A5UsedSpecies1) .* uvu(A5UsedSpecies2), nu*nu,1);
        
        val = dA1dk_fk_x * sparse(x) + dA2dk_fk_xx * xkronx + dA3dk_fk_ux * ukronx + dA4dk_fk_xu * xkronu + dA5dk_fk_uu * ukronu + dA6dk_fk_u * sparse(u); % fk_
        val = reshape(val, nx,nk) + dadk; % f_k
    end
end

function handle = d2fdx2Hidden(A2, B1, B2, b, A2UsedColumns, A2UsedSpecies2, vxInd)
% Constants
nx = size(A2,1);

% Reverse the internal subscripts in the UsedColumns linear index
[x1ind,x2ind] = ind2sub([nx,nx], A2UsedColumns);
A2UsedColumnsReverse = sub2ind([nx, nx], x2ind, x1ind);

% Return handle
handle = @d2fdx2;

    function val = d2fdx2(t, x, u)
%   d2fdx2 = 2*A2*(Ix kron diag(1/vx))
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        
        % Second-order derivative of (x kron x/vx)
        Ixkron1vx = sparse([A2UsedColumns; A2UsedColumns], [A2UsedColumns; A2UsedColumnsReverse], [vx(A2UsedSpecies2); vx(A2UsedSpecies2)], nx*nx,nx*nx);
        
        val = A2 * Ixkron1vx; % f_xx
        val = reshape(val, nx*nx, nx); % fx_x
    end
end

function handle = d2fdk2Hidden(nx, nk)
% Return handle
handle = @d2fdk2;

    function val = d2fdk2(t, x, u)
%   d2fdk2 = 0
        val = sparse(nx*nk, nk); % fk_k
    end
end

function handle = d2fdxdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, B1, B2, b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdxdk;

    function val = d2fdxdk(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(A2UsedColumns, A2UsedSpecies1, xvx(A2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(A2UsedColumns, A2UsedSpecies2, x(A2UsedSpecies1) .* vx(A2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(A3UsedColumns, A3UsedSpecies2, u(A3UsedSpecies1) .* vx(A3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(A4UsedColumns, A4UsedSpecies1, uvu(A4UsedSpecies2), nu*nx,nx);
        
        val = dA1dk_fk_x + dA2dk_fk_xx * (Ixkronxvx + xkron1vx) + dA3dk_fk_ux * ukron1vx + dA4dk_fk_xu * Ixkronuvu; % fx_k
    end
end

function handle = d2fdkdxHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, B1, B2, b, A2UsedColumns, A2UsedSpecies1, A2UsedSpecies2, A3UsedColumns, A3UsedSpecies1, A3UsedSpecies2, A4UsedColumns, A4UsedSpecies1, A4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @d2fdkdx;

    function val = d2fdkdx(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(A2UsedColumns, A2UsedSpecies1, xvx(A2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(A2UsedColumns, A2UsedSpecies2, x(A2UsedSpecies1) .* vx(A2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(A3UsedColumns, A3UsedSpecies2, u(A3UsedSpecies1) .* vx(A3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(A4UsedColumns, A4UsedSpecies1, uvu(A4UsedSpecies2), nu*nx,nx);
        
        val = dA1dk_fk_x + dA2dk_fk_xx * (Ixkronxvx + xkron1vx) + dA3dk_fk_ux * ukron1vx + dA4dk_fk_xu * Ixkronuvu; % fx_k
        val = spermute132(val, [nx,nk,nx], [nx*nx,nk]); % fk_x 
    end
end
