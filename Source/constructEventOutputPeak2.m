function eventFcn = constructEventOutputPeak2(m, outputIndex)
% 
% 
% function eventFcn = constructEventOutputPeak(m, outputIndex)
% 
% This function constructs an events function that looks for local maxima
% of the model output.  It stops integration when it finds one.
% 
% Inputs
%       m               -   The model from which output is collected
%       outputIndex     -   The index of the output for which the ISE will
%                           be computed
%       maxOrMin        -   Specifies whether to collect maxima or minima.
%                           This argument takes 3 values:   1 = maxima
%                                                           0 = all extrema
%                                                          -1 = minima
% 
% See MATLAB's ODE solver documentation (for instance, type 'doc ode15s')
% for more information on events functions.
% 
% See also:
%       constructEventSpeciesPeak
%       ode15s

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

% argument processing and error checking
if nargin<2 || isempty(outputIndex)
    outputIndex = 1;
end

assert(isnumeric(outputIndex), '', 'The second argument to constructEventOutputPeak must be an output index.');
assert(outputIndex <= m.nY, '', 'The output index %d is greater than the total number of outputs %d.', outputIndex, m.nY);


% return event function
eventFcn.events = @events;
eventFcn.update = @update;

M = length(outputIndex);

    function eventFcn = update(mNew)
        eventFcn = constructEventOutputPeak2(mNew, outputIndex);
    end

%% The events function that finds peaks in the output trajectory
    function [val isTerm direction] = events(t, x, u)
        
        val       = [m.c(outputIndex, :)*m.f(t, x, u);
                     m.c(outputIndex, :)*m.f(t, x, u)];
        
        isTerm    = ones(2*M, 1);
        direction = [-1*ones(M,1)
                      1*ones(M,1)];
        
%         val         = [val(1:2:end);       val(2:2:end)];
%         isTerm      = [isTerm(1:2:end);    isTerm(2:2:end)];
%         direction   = [direction(1:2:end); direction(2:2:end)];
    end
end