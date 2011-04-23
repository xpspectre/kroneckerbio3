function [func handled] = fixInputValueFunction(value, units, v0)

if nargin < 2
    units = [];
end

% TODO: units
if ~isempty(units); error('Units not yet implemented, tell David if you would like to volunteer'); end

% TODO: deal with the derivatives
if iscell(value); value = value{1}; end

handled = true;

if isnumeric(value)
    func = str2func(sprintf('@(t,q)repmat(%g, 1,numel(t))', value));
elseif isa(value, 'function_handle')
    func = value;
elseif ischar(value)
    if value(1) =='@'
        % Already a function handle
        func = str2func(value);
    else
        % Needs the handle
        value = ['@(t,q)' value];
        func = str2func(value);
    end
else
    error('KroneckerBio:InputValueFunction', 'The class of the input value function was invalid')
end
