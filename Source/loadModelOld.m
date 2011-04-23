function m = loadModel(varargin)
% myModel = loadModel(...)
% 
% 
% This function loads a KroneckerBio mass action model that is located in
% the tab-delimited text files RxnFile, ParamFile, and SpeciesFile.  The
% fourth argument specifies the name for this loaded model.
% 
% RxnFile     -   A tab-delimited text file containing the species and
%                 parameters participating in all reactions.
%
% ParamFile   -   A tab-delimited text file containing the names and values
%                 for all parameters.
% 
% SpeciesFile -   A tab-delimited text file containing a list of the names
%                 and initial conditions of all species, inputs and outputs
%                 for the model.  The Columns of the file are:
%
% Usage:
% This function can be called with either of two syntaxes:
% 
% myModel = loadModel('model1_reactions.txt', 'model1_parameters.txt', 'model1_species.txt', 'model1'); 
% or
% myModel = loadModel('model1');
% 
% In the latter case, loadModel assumes that the reactions, parameters and
% species files are those with the names <name>_reactions.txt,
% <name>_parameters.txt and <name>_species.txt, respectively.
%
% To generate text file stubs with the proper format, use generateModelFile
%
% 
% See also: 
%       generateModelFile

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

if nargin == 1
    Name        = varargin{1};
    RxnFile     = [Name '_reactions.txt'];
    ParamFile   = [Name '_parameters.txt'];
    SpeciesFile = [Name '_species.txt'];
elseif nargin == 3
    RxnFile     = varargin{1};
    ParamFile   = varargin{2};
    SpeciesFile = varargin{3};
    Name        = '';
elseif nargin == 4
    RxnFile     = varargin{1};
    ParamFile   = varargin{2};
    SpeciesFile = varargin{3};
    Name        = varargin{4};
else
    error('KroneckerBio:loadModel', 'Invalid arguments passed to loadModel');
end

assert(exist(RxnFile, 'file'), 'KroneckerBio:loadModel', 'Reaction file %s not found', RxnFile);
assert(exist(ParamFile, 'file'), 'KroneckerBio:loadModel', 'Parameter file %s not found', ParamFile);
assert(exist(SpeciesFile, 'file'), 'KroneckerBio:loadModel', 'Species file %s not found', SpeciesFile);

% Read in the reactions
[Rxn.R1Name, Rxn.R2Name, Rxn.PName, Rxn.FwdName, Rxn.RevName] ...
    = textread(RxnFile, '%q%q%q%q%q%*[^\n]',  'delimiter', ' \t',  'headerlines', 1);

% Read in the parameters
[Rate.Name,	 Rate.Value] = textread(ParamFile, '%q%f%*[^\n]', 'delimiter',' \t',  'headerlines', 1);

% Read in the species (excluding the outputs)
[Init.Name, Init.Conc, Init.Input] ...
    = textread(SpeciesFile, '%q%f%d%*[^\n]',  'delimiter', ' \t',  'headerlines', 1);

% Read in the output variables
OutputVarStart = 4; %Column where output listing begins
headerline     = textread(SpeciesFile, '%s', 1, 'delimiter', '\n');
matches        = regexp(headerline, '\S*', 'match');
matches        = matches{:}; %Decompose (1x1 cell of 1x(3+nY) cell) into just (1x(3+nY) cell)
Output.Names   = matches(OutputVarStart:end);
% Need the -1 because dlmread is 0 based
% TODO: Need to allow for space seperator in outputs values
Output.Values  = dlmread(SpeciesFile, '\t', 1, (OutputVarStart-1));

% check if 'Comments' or 'Comment' is the last field.  If so, kill it!
ind = strcmp(Output.Names, 'Comment') | strcmp(Output.Names, 'Comments');
if ind
    Output.Names = Output.Names(1:find(ind,1)-1);
    Output.Values = Output.Values(:, 1:find(ind,1)-1);
end

% Get the sizes and setup the model data structure
nX     = length(find(Init.Input==0));
nU     = length(find(Init.Input==1));
nY     = size(Output.Values, 2);
nP     = length(Rate.Value);
m      = initModel(nX, nU, nY, nP);
m.name = Name;

% Copy the parameter names and values into the model
m.pNames = Rate.Name;
m.p      = Rate.Value;

% Create a table that connects the species index to the index in the 
% state vector. Each value in stateVars indicates the index of the species
% that it corresponds to.
stateVars = find(Init.Input==0);

% Create a table that connects the species index to the index in the 
% input vector
inputVars = find(Init.Input==1);

% Copy the initial condition names and values into the model.  Input
% variables go into u and state variables go into x.
for i = 1:length(Init.Name)
    if Init.Input(i)
        ind               = toInputIndex(i);
        m.uNames{ind}     = Init.Name{i};
    else
        ind               = toStateIndex(i);
        m.ic(ind)         = Init.Conc(i);
        m.xNames{ind}     = Init.Name{i};
    end
    
end

% Loop through the lines of the reaction file
for i=1:length(Rxn.R1Name);
    % Match the entries in the reaction list to entries in the species file
    % and the parameter file.  This fails if there is a missing species or
    % parameter.  the indXXX are the index into x for regular
    % concentrations and u for inputs as denoted byt the inputXX flag.
    R1   = getReactantByName(Rxn.R1Name{i});
    R2   = getReactantByName(Rxn.R2Name{i});
    P    = getReactantByName(Rxn.PName{i});
    kFwd = getParamNumByName(Rxn.FwdName{i});
    kRev = getParamNumByName(Rxn.RevName{i});
    

    % Determine if this is a forward, revers, bidirectional reaction
    % or if there is no reaction
    if exists(kFwd)
        if exists(kRev)
            direction = 'b';
        else
            direction = 'f';
        end
    else
        if(exists(kRev))
            direction = 'r';
        else
            direction = '0';
        end
    end
    
    
    % To cut down on the cases we make sure that if there is only one
    % reactant it is R1 and if there is one one input it is also in R1
    if(R1.type == '0' || R2.type == 'u')
        tmp = R1;
        R1  = R2;
        R2  = tmp;
    end
    
    
    % Assemble a string to indicate the type of reaction
    type = [R1.type, R2.type, P.type, direction]; 
    
    r1   = R1.ind;
    r2   = R2.ind;
    p    = P.ind;
    f    = kFwd;
    r    = kRev;
    
    switch(type)
        % 0 <--> X
        case '00xf', m = addZero2One(m, p, f);            
        case '00xr', m = addOne2Zero(m, p, r);
        case '00xb', m = addOne2ZeroRev(m, p, r, f);
        
        % X <-> 0
        case 'x00f', m = addOne2Zero(m, r1, f);            
        case 'x00r', m = addZero2One(m, r1, r);
        case 'x00b', m = addOne2ZeroRev(m, r1, f, r);
        
        % X1 <-> X2
        case 'x0xf', m = addOne2One(m, r1,  p, f);         
        case 'x0xr', m = addOne2One(m,  p, r1, r);
        case 'x0xb', m = addOne2OneRev(m, r1,  p, f, r);
            
        % X1 + X2 <-> X3
        case 'xxxf', m = addTwo2One(m, r1, r2,  p, f);
        case 'xxxr', m = addOne2Two(m,  p, r1, r2, r);
        case 'xxxb', m = addTwo2OneRev(m, r1, r2,  p, f, r);

        % U <-> X
        case 'u0xf', m = addInputOne2One(m, r1, p, f);
        case 'u0xr', m = addOne2Zero(m, p, r);
        case 'u0xb', m = addInputOne2OneRev(m, r1, p, f, r);

        % X <-> U
        case 'x0uf', m = addOne2Zero(m, r1, f);
        case 'x0ur', m = addInputOne2One(m, p, r1, r);
        case 'x0ub', m = addInputOne2OneRev(m, p, r1, r, f);
            
        % U  + X1 <-> X2
        case 'uxxf', m = addInputTwo2One(m, r1, r2, p, f);
        case 'uxxr', m = addOne2One(m, p, r2, r);
        case 'uxxb', m = addInputTwo2OneRev(m, r1, r2, p, f, r);  
            
        otherwise
            switch(type(4))
                case '0', aStr = ' -- ';
                case 'f', aStr = ' -->';
                case 'r', aStr = '<-- ';
                case 'b', aStr = '<-->';
                otherwise
                    error('Invalid Reaction Flag');
            end

            rxnAssert(0, i, 'reactions of the form %s1 + %s2 %s %s3',...
                type(1), type(2),aStr, type(3));            
    end
end
 



% Copy the output matrix into the model along with the output variable
% names
assert(all(all(Output.Values(inputVars, :)==0)), 'kronecker:error', 'Input variables cannot be outputs');
m.c(:, 1:nX)  = Output.Values(stateVars, :)';
m.yNames      = Output.Names;



% Computes some quantities that will help us later, and add function
% handles
m = postProcessModel(m);


%--------------------------------------------------------------------------
% Nested functions
%--------------------------------------------------------------------------
    function ind = getParamNumByName(name)
         % Special case for 0 which means no parameter
        if(strcmp(name, '0'))
            ind = [];
            return
        end

        % find the parameter with the given name
        ind = find(strcmp(Rate.Name, name));

        assert(~isempty(ind), 'kronecker:modelerror', ...
        'Rate constant %s used in reaction %d could not be found in the parameter file.', name, i);

        assert(length(ind)==1, 'kronecker:modelerror', ...
            'Rate constant %s used in reaction %d appears multiple times in the parameter file', name, i);
    end


    function r = getReactantByName(name)
        % Special case for 0 which means no species
        if(strcmp(name, '0'))
            r.ind     = 0;
            r.isInput = 0;
            r.type    = '0';
            return
        end 
        
        speciesIndex = find(strcmp(Init.Name, name));

        assert(~isempty(speciesIndex), 'kronecker:modelerror', ...
        'Species constant %s used in reaction %d could not be found in the species file.', name, i);

        assert(length(speciesIndex)==1, 'kronecker:modelerror', ...
            'Species constant %s used in reaction %d appears multiple times in the species file', name, i);

        r.isInput = Init.Input(speciesIndex);
        
        if(r.isInput)
            r.ind = toInputIndex(speciesIndex);
            r.type = 'u';
        else
            r.ind = toStateIndex(speciesIndex);
            r.type = 'x';
        end
    end


    function ind=toInputIndex(speciesIndex)
        % Finds the position in the input vector corresponding to the 
        % position in the species list
        ind = find(inputVars == speciesIndex);
        assert(exists(ind), 'kronecker:modelerror', ...
            'Concentration index %d is not an input variable', speciesIndex);
    end


    function ind=toStateIndex(speciesIndex)
        % Finds the position in the state vector corresponding to the 
        % position in the species list
        if(speciesIndex == 0)
            ind = [];
        else
            ind = find(stateVars == speciesIndex);
            assert(exists(ind), 'kronecker:modelerror', ...
                'Concentration index %d is not a state variable', speciesIndex);
        end
    end
        
    function val=exists(obj)
        val = ~isempty(obj);
    end

    %%---------------------------------------------------------------------
    %% Error Reporting
    %%---------------------------------------------------------------------

    function rxnAssert(test, num, cause, varargin)
        causestr = sprintf(cause, varargin{:});
        assert(test, 'kronecker:modelerror', ...
            'Error loading Rxn %04d : %s', num, causestr);
    end
end

