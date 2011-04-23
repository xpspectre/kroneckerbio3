function m = AddCompartment(m, name, dimension, expressions, values, units)

% Clean-up inputs
if nargin < 6
    units = [];
end

% Validate compartment name
assert(isempty(regexp(name, '\.|,', 'once')), 'KroneckerBio:AddCompartment:CompartmentName', 'Compartment names cannot contain "." or ","')

% Increment counter
nv = m.add.nv + 1;
m.add.nv = nv;
m.add.Compartments = growCompartments(m.add.Compartments, m.add.nv);

% Add item
m.add.Compartments(nv).Name = fixName(name);
m.add.Compartments(nv).Expressions = fixExpression(expressions);
m.add.Compartments(nv).Values = fixCompartmentValues(values, fixUnits(units));
m.add.Compartments(nv).Dimension = dimension;

m.Ready = false;