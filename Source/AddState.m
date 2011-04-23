function m = AddState(m, name, compartment, initialValue, units)

% Clean-up inputs
if nargin < 5
    units = [];
    if nargin < 4
        initialValue = [];
    end
end

if isempty(initialValue)
    initialValue = 0;
end
if isempty(units)
    units = '';
end

% Increment counter
nxu = m.add.nxu + 1;
m.add.nxu = nxu;
m.add.Species = growSpecies(m.add.Species, nxu);

% Add item
m.add.Species(nxu).Name = fixName(name);
m.add.Species(nxu).Compartment = compartment;
m.add.Species(nxu).IsInput = false;
m.add.Species(nxu).Value = initialValue;
m.add.Species(nxu).Units = fixSpeciesUnits(units);

m.Ready = false;