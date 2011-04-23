function u = completeInputFunction(inputs)

nu = numel(inputs);

% Complete function
if nu == 0
    u = @(t,q)zeros(0,numel(t));
else
    % Extract input values
    inputs = cat(1, inputs.Value);
    
    % Stack each input function
    uEach = {inputs.Function}.';
    
    % Stack the parameters in a cell vector
    q = {inputs.Parameters}.';
    u = @(t)cellfun(@feval, uEach, repmat({t}, nu,1), q);
end