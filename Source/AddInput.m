function m = AddInput(m, name, compartment, value, parameters, units)

% Clean-up inputs
if nargin < 6
    units = [];
    if nargin <5
        parameters = [];
        if nargin < 4
            value = [];
        end
    end
end

% Increment counter
nxu = m.add.nxu + 1;
m.add.nxu = nxu;
m.add.Species = growSpecies(m.add.Species, nxu);

% Add item
m.add.Species(nxu).Name = fixName(name);
m.add.Species(nxu).Compartment = compartment;
m.add.Species(nxu).IsInput = true;
m.add.Species(nxu).Value.Function = fixInputValueFunction(value);
m.add.Species(nxu).Value.Parameters = fixInputParameters(parameters);
m.add.Species(nxu).Units = fixSpeciesUnits(units);

m.Ready = false;