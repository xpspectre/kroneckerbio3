function [useControls, nTq] = fixUseControls(useControls, useModelInputs, nCon, nqModel, nqCon)
%TODO: process linear indexes
if useModelInputs
    if iscell(useControls)
        assert(numel(useControls)==1)
        useControls = useControls{1};
    end
    nTq = nnz(useControls);
else%~useModelInputs
    if ~iscell(useControls)
        % Repeat across experiments
        if islogical(useControls)
            assert(all(numel(useControls) == nqCon), 'KroneckerBio:UseControls:LogicalLength', 'If a logical vector is supplied for UseControls when UseModelInputs is false, then the number of input control parameters must be the same for all conditions.')
            useControls = repmat({useControls}, nCon,1);
        else%isindexes
            temp = cell(nCon,1);
            for iCon = 1:nCon
                temp{iCon} = false(nqCon(iCon));
                temp{iCon}(useControls) = true;
            end
            useControls = temp;
        end
    end
    
    % Repeat a singular array
    if numel(useControls) ~= nCon
        assert(numel(useControls) == 1 , 'KroneckerBio:UseControls:UseControlsLength', 'If UseControls is provided as a cell array, its length must be equal to the length of con or 1')
        useControls = repmat(useControls, nCon,1);
    end
    
    nTq = 0;
    for iCon = 1:nCon
        nTq = nTq + nnz(useControls{iCon});
    end
end

