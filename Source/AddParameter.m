function m = AddParameter(m, name, value, units)

% Clean-up inputs
if nargin < 4
    units = [];
end

% Increment counter
nk = m.add.nk + 1;
m.add.nk = nk;
m.add.Parameters = growParameters(m.add.Parameters, nk);

% Add item
m.add.Parameters(nk).Name = fixName(name);
m.add.Parameters(nk).Value = fixParameterValue(value, fixSpeciesUnits(units));

m.Ready = false;