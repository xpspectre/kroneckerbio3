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

nTop = numel(m);
nCon = size(obj, 1);
nObj = size(obj, 1);
newObj = Gzero([nObj,nCon,nTop]);
for iTop = 1:nTop
    for iCon = 1:nCon
        for iObj = 1:nObj
            newObj(iObj,iCon,iTop) = pastestruct(Gzero(m(iTop)), obj(iObj,iCon,iTop).Update(m(iTop), con(iCon,iTop), UseParams{iTop}, UseICs{iTop}(:,iCon), UseControls{iCon,iTop}));
        end
    end
end
