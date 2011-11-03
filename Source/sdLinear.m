function handle = sdLinear(a, b)

handle = @sdHandle;


    function [sigma dsigmady d2sigmady2] = sdHandle(t, yInd, yVal)
        sigma = a + b .* yVal;
        var = sigma.^2;
        
        if nargout >= 2
            % First derivative
            dsigmady = 0.5 .* var.^-0.5 .* propVar .* 2 .* yVal;
            if nargout >=3
                % Second derivative
                d2sigmady2 = 0.5 .* var.^-0.5 .* propVar .* 2 + -0.25 .* var.^-1.5 .* propVar.^2 .* 4 .* yVal.^2;
            end
        end
    end
end