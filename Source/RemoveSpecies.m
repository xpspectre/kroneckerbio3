function m = RemoveSpecies(m, name)

% Find added instances of this name
indAdd = strcmp(name, {m.add.Species.Name});

% Remove all mention of this output
m.Species(strcmp(name, {m.Species.Name})) = [];
m.add.Species(indAdd) = [];
m.add.ny = m.add.ny - nnz(indAdd);

m.Ready = false;