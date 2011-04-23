function absTol = fixAbsTol(absTol, order, integrateObj, nX, nCon, useAdjoint, useParams, useICs, useModelICs, selfSensitivities)
%FIXABSTOL Standardize the presentation of AbsTol
%
%   There are many ways to present the absolute integration tolerance to
%   Kronecker Bio. This function is the processing center for these
%   different presentations. The many inputs define the type of problem
%   that the AbsTol is needed for; essentially, the length. The standard
%   presentation is a cell array of vectors of the correct length. Each
%   cell corresponds to an experimental condition.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 10
    selfSensitivities = false;
    if nargin < 9
        useModelICs = [];
        if nargin < 8
            useICs = [];
            if nargin < 7
                useParams = [];
                if nargin < 6
                    useAdjoint = [];
                end
            end
        end
    end
end

% The default AbsTol
if isempty(absTol) || (~iscell(absTol) && all(vec(isnan(absTol))))
    absTol = 1e-9;
end

% Constants
nVP = sum(useParams);
nVX = sum(sum(useICs));
nV = nVP + nVX;

abstol = cell(nCon,1);

switch order
    case 1
        if ~selfSensitivities
            % System integration
            if ~iscell(absTol)
                % AbsTol is the same for all conditions
                if isscalar(absTol)
                    % Copy the value to all species and conditions
                    for i = 1:nCon
                        abstol{i} = zeros(nX+integrateObj(i),1) + absTol;
                    end
                    %absTol = repmat({zeros(nX+integrateObj,1) + absTol}, nCon,1);
                else %isvector
                    % Copy the vector. Only if integrateObj is zero allow the
                    % wrong length vector to be supplied.
                    assert(all(integrateObj == 0) || all(length(absTol) == nX+integrateObj), 'KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
                    abstol = repmat({absTol(1:nX+integrateObj)}, nCon,1);
                end
            end
        else %selfSensitivities
        end
    case 2
        if ~selfSensitivities
            % Sensitivity integration
            if ~iscell(absTol)
                % AbsTol is the same for all conditions
                if useAdjoint
                    if isscalar(absTol)
                        % Copy the value to all species and conditions
                        abstol = repmat({zeros(nX+1+nX+nV,1) + absTol}, nCon,1);
                    else %isvector
                        % Copy the vector. Only if integrateObj is zero allow the
                        % wrong length vector to be supplied.
                        assert(length(absTol) == nX+1+nX+nV, 'KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
                        abstol = repmat({absTol}, nCon,1);
                    end
                else %~useAdjoint
                    if isscalar(absTol)
                        % AbsTol is starting scalar
                        for i = 1:nCon
                            if useModelICs
                                % Number of parameters is constant across conditions
                                inV = nV;
                            else
                                % Number of IC parameters may vary across conditions
                                inV = sum(useICs(:,i)) + nVP;
                            end
                            abstol{i} = zeros(nX+integrateObj(i)+nX*inV+integrateObj(i)*inV,1) + absTol;
                        end
                    elseif numel(absTol) == nX*(nV+1)
                        if any(integrateObj)
                            error('KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective and gradient')
                        else
                            if nCon == 1 || useModelICs || nVX == 0
                                % Number of parameters is the same across conditions.
                                abstol = repmat({vec(absTol)}, nCon,1);
                            else
                                error('KroneckerBio:AbsTol:VectorWithMultipleVariableExperimentICs', 'A vector cannot be provided for AbsTol if UseModelICs is false and there are variable ICs')
                            end
                        end
                    elseif numel(absTol) == nX+1+nX*nV+nV
                        for i = 1:nCon
                            if nCon == 1 || useModelICs || nVX == 0
                                % Number of parameters is the same across conditions.
                                abstol{i} = [absTol(1:nX);
                                    repmat(absTol(nX+1), integrateObj(i),1);
                                    absTol(nX+1+1:nX+1+nX*nV);
                                    repmat(absTol(nX+1+nX*nV+1:nX+1+nX*nV+nV), integrateObj(i),1)];
                            else
                                error('KroneckerBio:AbsTol:VectorWithMultipleVariableExperimentICs', 'A vector cannot be provided for AbsTol if UseModelICs is false and there are variable ICs')
                            end
                        end
                    elseif all(numel(absTol) == nX+integrateObj+nX*nV+nV*integrateObj)
                        % AbsTol is the correct length for all experiments
                        abstol = repmat({absTol}, nCon,1);
                    else
                        error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
                    end
                end
            end
        else %selfSensitivities
            if ~iscell(absTol)
                if useAdjoint
                    error('Adjoint method with optimal AbsTol not yet supported.')
                else
                    if isscalar(absTol)
                        % AbsTol is starting scalar
                        for i = 1:nCon
                            if useModelICs
                                % Number of parameters is constant across conditions
                                inV = nV;
                            else
                                % Number of IC parameters may vary across conditions
                                inV = sum(useICs(:,i)) + nVP;
                            end
                            abstol{i} = zeros(nX+integrateObj(i)+nX*inV+integrateObj(i)*inV+nX*nX+nX*inV*nX,1) + absTol;
                        end
                    elseif all(numel(absTol) == nX+integrateObj+nX*nV+nV*integrateObj+nX*nX+nX*nV*nX)
                        % AbsTol is the correct length for all experiments
                        abstol = repmat({absTol}, nCon,1);
                    else
                        error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
                    end
                end
            end
        end
    case 3
        % Hessian integration
        error('An unsupported order was passed to fixAbsTol.')
    otherwise
        error('An unsupported order was passed to fixAbsTol.')
end

if ~iscell(absTol) % As of now, no changes are make to a cellular absTol
    absTol = abstol;
end