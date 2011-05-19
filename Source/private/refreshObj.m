function newObj = refreshObj(m, con, obj, UseParams, UseICs, UseControls)
%REFRESHOBJ Update Kronecker Bio objective functions according to a new
%   model
%
%   Whenever either the topolgy or parameters of a model change, the
%   objective function structures should be updated because they may
%   depend on the changes. This function loops over the objective functions
%   in order to accomplish this common task.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Make options consistent as cell arrays
if ~iscell(UseParams)
    UseParams = {UseParams};
end
if ~iscell(UseICs)
    UseParams = {UseICs};
end
if ~iscell(UseControls)
    UseParams = {UseControls};
end


nTop = numel(m);
nCon = size(con, 1);
nObj = size(obj, 1);
newObj = Gzero([nObj,nCon,nTop]);

for iTop = 1:nTop
    for iCon = 1:nCon
        for iObj = 1:nObj
            if size(UseICs{iTop},2) == 1 && size(UseControls,1) == 1 && nCon > 1
                % UseModelICs == true && UseModelInputs == true
                useICsi = UseICs{iTop};
                useControlsi = UseControls{1,iTop};
            elseif size(UseICs{iTop},2) == 1 && nCon > 1
                % UseModelsICs == true && UseModelInputs == false
                useICsi = UseICs{iTop};
                useControlsi = UseControls{iCon,iTop};
            elseif size(UseControls,1) == 1 && nCon > 1
                % UseModelsICs == false && UseModelInputs == true
                useICsi = UseICs{iTop}(:,iCon);
                useControlsi = UseControls{iCon,iTop};
            else% (size(UseICs{iTop},2) ~= 1 && size(UseControls,1) ~= 1) || nCon == 1
                % UseModelsICs == false && UseModelInputs == false
                useICsi = UseICs{iTop}(:,iCon);
                useControlsi = UseControls{iCon,iTop};
            end
            
            % Refresh obj
            newObj(iObj,iCon,iTop) = pastestruct(Gzero(m(iTop)), obj(iObj,iCon,iTop).Update(m(iTop), con(iCon,iTop), UseParams{iTop}, useICsi, useControlsi));
        end
    end
end
