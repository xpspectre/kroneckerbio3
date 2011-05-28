function absTol = fixAbsTol(absTol, order, integrateObj, nx, nCon, useAdjoint, useParams, useICs, useModelICs, selfSensitivities)
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
nTk = sum(useParams);
nTx = sum(sum(useICs));
nT = nTk + nTx;

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
                        abstol{i} = zeros(nx+integrateObj(i),1) + absTol;
                    end
                    %absTol = repmat({zeros(nx+integrateObj,1) + absTol}, nCon,1);
                else %isvector
                    % Copy the vector. Only if integrateObj is zero allow the
                    % wrong length vector to be supplied.
                    assert(all(integrateObj == 0) || all(numel(absTol) == nx+integrateObj), 'KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
                    abstol = repmat({absTol(1:nx+integrateObj)}, nCon,1);
                end
            end
        else %selfSensitivities
            error('Self sensitivities not yet supported for first order.')
        end
    case 2
        if ~selfSensitivities
            % Sensitivity integration
            if ~iscell(absTol)
                % AbsTol is the same for all conditions
                if useAdjoint
                    if isscalar(absTol)
                        % Copy the value to all species and conditions
                        for i = 1:nCon
                            abstol{i} = zeros(nx+integrateObj(i)+nx+nT,1) + absTol;
                        end
                    else %isvector
                        % Copy the vector. Only if integrateObj is zero allow the
                        % wrong length vector to be supplied.
                        assert(all(numel(absTol) == nx+integrateObj+nx+nT), 'KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
                        abstol = repmat({absTol}, nCon,1);
                    end
                else %~useAdjoint
                    if isscalar(absTol)
                        % AbsTol is starting scalar
                        for i = 1:nCon
                            if useModelICs
                                % Number of parameters is constant across conditions
                                inT = nT;
                            else
                                % Number of IC parameters may vary across conditions
                                inT = sum(useICs(:,i)) + nTk;
                            end
                            abstol{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT,1) + absTol;
                        end
                    elseif numel(absTol) == nx*(nT+1)
                        if any(integrateObj)
                            error('KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective and gradient')
                        else
                            if nCon == 1 || useModelICs || nTx == 0
                                % Number of parameters is the same across conditions.
                                abstol = repmat({vec(absTol)}, nCon,1);
                            else
                                error('KroneckerBio:AbsTol:VectorWithMultipleVariableExperimentICs', 'A vector cannot be provided for AbsTol if UseModelICs is false and there are variable ICs')
                            end
                        end
                    elseif numel(absTol) == nx+1+nx*nT+nT
                        for i = 1:nCon
                            if nCon == 1 || useModelICs || nTx == 0
                                % Number of parameters is the same across conditions.
                                abstol{i} = [absTol(1:nx);
                                    repmat(absTol(nx+1), integrateObj(i),1);
                                    absTol(nx+1+1:nx+1+nx*nT);
                                    repmat(absTol(nx+1+nx*nT+1:nx+1+nx*nT+nT), integrateObj(i),1)];
                            else
                                error('KroneckerBio:AbsTol:VectorWithMultipleVariableExperimentICs', 'A vector cannot be provided for AbsTol if UseModelICs is false and there are variable ICs')
                            end
                        end
                    elseif all(numel(absTol) == nx+integrateObj+nx*nT+nT*integrateObj)
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
                                inT = nT;
                            else
                                % Number of IC parameters may vary across conditions
                                inT = sum(useICs(:,i)) + nTk;
                            end
                            abstol{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT+nx*nx+nx*inT*nx,1) + absTol;
                        end
                    elseif all(numel(absTol) == nx+integrateObj+nx*nT+nT*integrateObj+nx*nx+nx*nT*nx)
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