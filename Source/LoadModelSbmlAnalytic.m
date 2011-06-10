function kronModel = LoadModelSbmlAnalytic(simbioModel, uNames, yNames, yMembers, yValues, opts)
%LoadModelSbmlAnalytic converts a Simbiology model into a pseudo-kronecker
%   model, which interacts with the Kronecker Bio toolbox exactly like a
%   Kronecker model, but with substantial performance reductions.
% 
%   m = LoadModelSbmlAnalytic(SimbioModel, uNames, yNames, yMembers, opts)
%
%   Inputs
%       simbioModel - A Simbiology model object or SBML model file name
%       uNames   - The names of the species you want to designate as the
%                  input. Can be a string, a vector of symbolics, a cell
%                  vector of strings, or a scalar vector corresponding to
%                  the indexes of the input variables. Default = {}
%       yNames   - The names you want to give to the outputs. If yMembers
%                  is not provided, the names are also interpreted as the
%                  species that will be represented by the outputs.
%                  Naturally, if yMembers is not provided, then each
%                  output will only reflect the concentration of a single
%                  species. Can be a string, a vector of symbolics, a cell
%                  vector of strings, or a scalar vector corresponding to
%                  the index of the species if yMembers is not provided.
%                  Default = 1:nX
%       yMembers - The names of the species you want to be included in each
%                  output. Must be a cell array of cell arrays of strings.
%       opts     - Scalar structure of options
%           Verbose [ true | {false} ]
%               Print progress to command window
%           Order [ 0 | 1 | {2} | 3 ]
%               Determines how deep the derivatives should be taken with
%               respect to x and p. Each level increases the cost
%               exponentially, but increases the number of Kronecker Bio
%               functions that can be run on the model.
%
%   Outputs
%       m - A psuedo-kronecker model
%
%   Limitations:
%   Not all Simbiology features are compatible with this converter. This
%   function ignores any events, some rules, and functions of the model.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 6
    opts = [];
    if nargin < 5
        yValues = [];
        if nargin < 4
            yMembers = [];
            if nargin < 3
                yNames = [];
                if nargin < 2
                    uNames = [];
                end
            end
        end
    end
end

% Load sbml if file is provided
if ischar(simbioModel)
    simbioModel = sbmlimport(simbioModel);
end

% Use sbml2Symbolic to convert an SBML model to a symbolic model
symModel = simbio2Symbolic(simbioModel, opts);

% Use symbolic2Kronecker to convert a symbolic model to a psuedo kronecker
% model
kronModel = symbolic2PseudoKronecker(symModel, uNames, yNames, yMembers, yValues, opts);