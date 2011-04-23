function m = AddOutput(m, name, expressions, values, units)

% Clean-up inputs
if nargin < 5
    units = [];
end

% Increment counter
ny = m.add.ny + 1;
m.add.ny = ny;
m.add.Outputs = growOutputs(m.add.Outputs, ny);

% Add item
m.add.Outputs(ny).Name = fixName(name);
m.add.Outputs(ny).Expressions = fixExpression(expressions);
m.add.Outputs(ny).Values = fixOutputValues(values, fixUnits(units));

m.Ready = false;