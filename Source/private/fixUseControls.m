function [useControls, nTq] = fixUseControls(useControls, useModelInputs, nCon, nqModel, nqCon)
%TODO: process linear indexes
if useModelInputs
    if iscell(useControls)
        assert(numel(useControls)==1)
        useControls = useControls{1};
    end
    nTq = nnz(useControls);
else
    if iscell()
        %TODO:repeat singular cell array
        nTq = 0;
        for iCon = 1:nCon
            nTq = nTq + nnz(useControls{iCon});
        end
    else
        %TODO:noncell processes
    end
end

