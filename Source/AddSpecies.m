function m = AddSpecies(m, name, compartment, value, units, isInput, parameters)

% Clean-up inputs
if nargin < 7
    parameters = [];
    if nargin < 6
        isInput = [];
        if nargin < 5
            units = [];
            if nargin < 4
                value = [];
            end
        end
    end
end

% Default inputs
if isempty(isInput)
    isInput = false;
end

% Simply defer calculation to the approriate functions
if isInput
    m = AddInput(m, name, compartment, value, parameters, units);
else
    m = AddState(m, name, compartment, value, units);
end
