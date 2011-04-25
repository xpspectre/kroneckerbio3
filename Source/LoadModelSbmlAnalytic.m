function kronModel = LoadModelSbmlAnalytic(SimbioModel, Unames, Ynames, Ymembers, opts)
%LoadModelSbmlAnalytic converts a Simbiology model into a pseudo-kronecker
%   model, which interacts with the Kronecker Bio toolbox exactly like a
%   Kronecker model, but with substantial performance reductions.
% 
%   m = loadModelSbml(SimbioModel, Unames, Ynames, Ymembers, opts)
%
%   Inputs
%       SimbioModel - A Simbiology model object or SBML model file name
%       Unames   - The names of the species you want to designate as the
%                  input. Can be a string, a vector of symbolics, a cell
%                  vector of strings, or a scalar vector corresponding to
%                  the indexes of the input variables. Default = {}
%       Ynames   - The names you want to give to the outputs. If Ymembers
%                  is not provided, the names are also interpreted as the
%                  species that will be represented by the outputs.
%                  Naturally, if Ymembers is not provided, then each
%                  output will only reflect the concentration of a single
%                  species. Can be a string, a vector of symbolics, a cell
%                  vector of strings, or a scalar vector corresponding to
%                  the index of the species if Ymembers is not provided.
%                  Default = 1:nX
%       Ymembers - The names of the species you want to be included in each
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
%   function ignores any events, rules, and functions of the model.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 5
    opts = [];
    if nargin < 4
        Ymembers = [];
        if nargin < 3
            Ynames = [];
            if nargin < 2
                Unames = [];
            end
        end
    end
end

% Load sbml if file is provided
if ischar(SimbioModel)
    SimbioModel = sbmlimport(SimbioModel);
end

% Use sbml2Symbolic to convert an SBML model to a symbolic model
symModel = simbio2Symbolic(SimbioModel, opts);

% Use symbolic2Kronecker to convert a symbolic model to a psuedo kronecker
% model
kronModel = symbolic2PseudoKronecker(symModel, Unames, Ynames, Ymembers, opts);