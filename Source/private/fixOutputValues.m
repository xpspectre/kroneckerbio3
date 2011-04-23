function values = fixOutputValues(values, units)

for i = 1:numel(units)
    if ~isempty(units{i}); error; end
end