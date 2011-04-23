function m = AddReaction(m, name, compartment, reactant1, reactant2, product1, product2, kForward, kReverse)

% Clean-up
if nargin < 9
    kReverse = [];
    if nargin < 8
        kForward = [];
    end
end

% Standardize names
[name1 name2] = fixReactionName(name, kForward, kReverse);
reactant1 = fixReactionSpecies(reactant1);
reactant2 = fixReactionSpecies(reactant2);
product1  = fixReactionSpecies(product1);
product2  = fixReactionSpecies(product2);
kForward  = fixReactionParameter(kForward);
kReverse  = fixReactionParameter(kReverse);
compartment = fixReactionCompartment(compartment);

% Increment counter
nr = m.add.nr + 1;
m.add.nr = nr;
m.add.Reactions = growReactions(m.add.Reactions, nr);

% Add item
m.add.Reactions(nr).Names = {name1, name2};
m.add.Reactions(nr).Compartment  = compartment;

if ~isempty(reactant1)
    m.add.Reactions(nr).Reactants = {reactant1, reactant2};
else
    % Prevent an order of {0, reactant}
    m.add.Reactions(nr).Reactants = {reactant2, reactant1};
end

if ~isempty(product1)
    m.add.Reactions(nr).Products = {product1, product2};
else
    % Prevent an order of {0, product}
    m.add.Reactions(nr).Products = {product2, product1};
end

m.add.Reactions(nr).Parameters = {kForward, kReverse};

m.Ready = false;