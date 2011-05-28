function m = symbolic2PseudoKronecker(SymModel, Unames, Ynames, Ymembers, opts)
%SYMBOLIC2KRONECKER converts a symbolic model into a pseudo-kronecker
%   model, which interacts with the Kronecker Bio toolbox exactly like a
%   Kronecker model, but with substantial performance reductions.
% 
%   m = symbolic2Kronecker(SymModel, Unames, Ynames, Ymembers, opts)
% 
%   Inputs
%       SymModel - A symbolic model
%       Unames   - The names of the species you want to designate as the
%                  input. Can be a string, a vector of symbolics, a cell
%                  vector of strings, or a scalar vector corresponding to
%                  the indexes of the input variables.
%       Ynames   - The names you want to give to the outputs. If Ymembers
%                  is not provided, the names are also interpreted as the
%                  species that will be represented by the outputs.
%                  Naturally, if Ymembers is not provided, then each
%                  output can only reflect the concentration of a single
%                  species. Can be a string, a vector of symbolics, a cell
%                  vector of strings, or a scalar vector corresponding to
%                  the index of the species if Ymembers is not provided.
%                  Default = 1:nX
%       Ymembers - The names of the species you want to be included in each
%                  output. Must be a cell array of cell arrays of strings.
%       opts     - Optional function options
%           .Verbose [ true | {false} ]
%               Print progress to command window
%           .Order [ 0 | 1 | {2} | 3 ]
%               Determines how deep the derivatives should be taken with
%               respect to x and p. Each level increases the cost
%               exponentially, but increases the number of Kronecker Bio
%               functions that can be run on the model.
%
%   Outputs
%       m - A psuedo-kronecker model

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.


%TODO: Error checking that Unames and Ynames are in Xnames

%% Work-up
% Clean up inputs
if nargin < 5
    opts = [];
    if nargin < 4
        Ymembers = [];
        if nargin < 3
            Ynames = [];
        end
    end
end

%Options for InputType and displaying progress
defaultOpts.Verbose = 0;
defaultOpts.Order   = 2;

opts = mergestruct(defaultOpts, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 1: Generate Differential Equations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xnames      = SymModel.Xnames;
Knames      = SymModel.Knames;
x0          = SymModel.x0;
k           = SymModel.k;
Xsyms       = SymModel.Xsyms;
Ksyms       = SymModel.Ksyms;
f           = SymModel.f;

% Add compartment volumes as parameters
for i = 1:nv
    nk = nk + 1;
    kNames{nk,1} = vNames{i};
    k(nk,1) = vValues(i);
    kNicestrs{nk,1} = vNicestrs{i};
    kSyms(nk,1) = vSyms(i);
end

if opts.Order >= 1
%Gradient of r with respect to x
drdx = jacobian(Rsyms, xSyms);

%Gradient of r with respect to k
drdk = jacobian(Rsyms, kSyms);

%Gradient of f with respect to x
dfdx = S*drdx;

%Gradient of f with respect to k
dfdk = S*drdk;
else
dfdx        = '';
dfdk        = '';
end

if opts.Order >= 2
%Gradient of drdx with respect to x
d2rdx2 = jacobian(vec(drdx), xSyms);

%Gradient of drdk with respect to k
d2rdk2 = jacobian(vec(drdk), kSyms);

%Gradient of drdx with respect to k
d2rdkdx = jacobian(vec(drdx), kSyms);

%Gradient of drdk with respect to x
d2rdxdk = jacobian(vec(drdk), xSyms);

%Gradient of dfdx with respect to x
d2fdx2 = S*reshape(d2rdx2, nr,nx*nx);
d2fdx2 = reshape(d2fdx2, nx*nx,nx);

%Gradient of dfdk with respect to k
d2fdk2 = S*reshape(d2rdk2, nr,nk*nk);
d2fdk2 = reshape(d2fdk2, nx*nk,nk);

%Gradient of dfdx with respect to k
d2fdkdx = S*reshape(d2rdkdx, nr,nx*nk);
d2fdkdx = reshape(d2fdkdx, nx*nx,nk);

%Gradient of dfdk with respect to x
d2fdxdk = S*reshape(d2rdxdk, nr,nk*nx);
d2fdxdk = reshape(d2fdxdk, nx*nk,nx);
else
d2fdx2      = '';
d2fdk2      = '';
d2fdkdx     = '';
d2fdxdk     = '';
end

if opts.Order >= 3
%Gradient of d2rdx2 with respect to x
d3rdx3 = jacobian(vec(d2rdx2), xSyms);

%Gradient of d2rdx2 with respect to k
d3rdkdx2 = jacobian(vec(d2rdx2), kSyms);

%Gradient of d2fdx2 with respect to x
d3fdx3 = S*reshape(d3rdx3, nr,nx*nx*nx);
d3fdx3 = reshape(d3fdx3, nx*nx*nx,nx);

%Gradient of d2fdx2 with respect to k
d3fdkdx2 = S*reshape(d3rdkdx2, nr,nx*nx*nk);
d3fdkdx2 = reshape(d3fdkdx2, nx*nx*nx,nk);
else
d3fdx3      = '';
d3fdkdx2    = '';
end

nX          = length(Xsyms);
nK          = length(Ksyms);

Xnicestrs = cell(nX,1);
for i = 1:nX
    Xnicestrs{i} = char(Xsyms(i));
end

Knicestrs = cell(nK,1);
for i = 1:nK
    Knicestrs{i} = char(Ksyms(i));
end

% Clear symbolic variables to save space
clear SymModel Xsyms Ksyms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 2: Include the inputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deconvolute the input
% Standardize Unames as a cell array of strings
if isa(Unames, 'char')
    Unames = {Unames};
    Uind = find(strcmp(Xnames, Unames{1}));
elseif isa(Unames, 'sym')
    temp = Unames;
    Unames = cell(length(Unames),1);
    for i = 1:length(Unames)
        Unames{i} = char(temp(i));
        Uind = find(strcmp(Xnames, Unames{i}));
    end
elseif isa(Unames, 'numeric')
    Uind = Unames;
    Unames = Xnames(Uind);
elseif isa(Unames, 'cell')
    Uind = zeros(length(Unames),1);
    for i = 1:length(Unames)
        Uind(i) = find(strcmp(Xnames, Unames{i}));
    end
end
Unicestr = Xnicestrs(Uind);
nU = length(Uind);

%% Delete the rows and columns corresponding to the input

Xnames(Uind)    = [];
x0(Uind)        = [];

Xnicestrs(Uind) = [];

% f
f(Uind,:)       = [];

% dfdx, keeping dfdu
if ~isempty(dfdx)
    dfdx(Uind,:)    = [];
    dfdu = dfdx(:,Uind);
    dfdx(:,Uind) = [];
end

% dfdk
if ~isempty(dfdk)
    dfdk(Uind,:)    = [];
end

% d2fdx2
if ~isempty(d2fdx2)
    d2fdx2 = reshape(d2fdx2, nX,nX,nX);
    d2fdx2(Uind,:,:) = [];
    d2fdx2(:,Uind,:) = [];
    d2fdx2(:,:,Uind) = [];
    d2fdx2 = reshape(d2fdx2, (nX-nU)*(nX-nU),(nX-nU));
end

% d2fdk2
if ~isempty(d2fdk2)
    d2fdk2 = reshape(d2fdk2, nX,nK,nK);
    d2fdk2(Uind,:,:) = [];
    d2fdk2 = reshape(d2fdk2, (nX-nU)*nK,nK);
end

% d2fdkdx
if ~isempty(d2fdkdx)
    d2fdkdx = reshape(d2fdkdx, nX,nX,nK);
    d2fdkdx(Uind,:,:) = [];
    d2fdkdx(:,Uind,:) = [];
    d2fdkdx = reshape(d2fdkdx, (nX-nU)*(nX-nU),nK);
end

% d2fdxdk
if ~isempty(d2fdxdk)
    d2fdxdk = reshape(d2fdxdk, nX,nK,nX);
    d2fdxdk(Uind,:,:) = [];
    d2fdxdk(:,:,Uind) = [];
    d2fdxdk = reshape(d2fdxdk, (nX-nU)*nK,(nX-nU));
end

% d3fdx3
if ~isempty(d3fdx3)
    d3fdx3 = reshape(d3fdx3, nX,nX,nX,nX);
    d3fdx3(Uind,:,:,:) = [];
    d3fdx3(:,Uind,:,:) = [];
    d3fdx3(:,:,Uind,:) = [];
    d3fdx3(:,:,:,Uind) = [];
    d3fdx3 = reshape(d3fdx3, (nX-nU)*(nX-nU)*(nX-nU),(nX-nU));
end

% d3fdkdx2
if ~isempty(d3fdkdx2)
    d3fdkdx2 = reshape(d3fdkdx2, nX,nX,nX,nK);
    d3fdkdx2(Uind,:,:,:) = [];
    d3fdkdx2(:,Uind,:,:) = [];
    d3fdkdx2(:,:,Uind,:) = [];
    d3fdkdx2 = reshape(d3fdkdx2, (nX-nU)*(nX-nU)*(nX-nU),nK);
end

% Update size of x
nX = nX - nU;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 3: Include the outputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deconvolute the output
% Standardize  Ynames as a cell array of strings
% and Ymembers as a cell array of cell arrays of strings
if isempty(Ymembers)
    % Ynames are the Ymembers
    if isempty(Ynames)
        % All species are outputs
        Ynames = 1:nX;
    end
    
    if isa(Ynames, 'char')
        nY = 1;
        Ynames = {Ynames};
        C = zeros(nY, nX);
        C(strcmp(Xnames, Ynames{1})) = 1;
    elseif isa(Ynames, 'sym')
        nY = length(Ynames);
        temp = Ynames;
        Ynames = cell(nY,1);
        C = zeros(nY, nX);
        for i = 1:nY
            Ynames{i} = char(temp(i));
            C(i, strcmp(Xnames, Ynames{i})) = C(i, strcmp(Xnames, Ynames{i})) + 1;
        end
    elseif isa(Ynames, 'numeric')
        nY = length(Ynames);
        Ynames = Xnames(Ynames);
        C = zeros(nY, nX);
        for i = 1:nY
            C(i, strcmp(Xnames, Ynames{i})) = C(i, strcmp(Xnames, Ynames{i})) + 1;
        end
    else%if isa(Ynames, 'cell')
        nY = length(Ynames);
        Ynames = vec(Ynames);
        C = zeros(nY, nX);
        for i = 1:nY
            C(i, strcmp(Xnames, Ynames{i})) = C(i, strcmp(Xnames, Ynames{i})) + 1;
        end
    end
else
    % Ymembers is provided
    if isa(Ynames, 'char')
        nY = 1;
        Ynames = {Ynames};
    elseif isa(Ynames, 'sym')
        nY = length(Ynames);
        temp = Ynames;
        Ynames = cell(nY,1);
        for i = 1:nY
            Ynames{i} = char(temp(i));
        end
    elseif isa(Ynames, 'numeric')
        nY = length(Ynames);
        Ynames = Xnames(Ynames);
    else%if isa(Ynames, 'cell')
        nY = length(Ynames);
        Ynames = vec(Ynames);
    end
    
    % Ymembers must be a cell array of cells
    C = zeros(nY, nX);
    for i = 1:nY
        nMem = length(Ymembers{i});
        for j = 1:nMem
            C(i, strcmp(Xnames, Ymembers{i}{j})) = C(i, strcmp(Xnames, Ymembers{i}{j})) + 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 4: Format variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replace Verbose Names with Systematic names
if opts.Verbose; fprintf('Converting symbolics to strings...\n'); end
% Convert the symbolics into strings
f           = strtrim(evalc('disp(f)'));
if ~isempty(dfdx)
    dfdx    = strtrim(evalc('disp(dfdx)'));
    if nU > 0
        % Symbolic toolbox fails if Uind is empty
        dfdu    = strtrim(evalc('disp(dfdu)'));
    else
        % Use "zeros" if it will be empty
        dfdu    = sprintf('zeros(%d,0)',nX);
    end
else
    dfdu = '';
end
if ~isempty(dfdk)
    dfdk     = strtrim(evalc('disp(dfdk)'));
end
if ~isempty(d2fdx2)
    d2fdx2   = strtrim(evalc('disp(d2fdx2)'));
end
if ~isempty(d2fdk2)
    d2fdk2   = strtrim(evalc('disp(d2fdk2)'));
end
if ~isempty(d2fdkdx)
    d2fdkdx  = strtrim(evalc('disp(d2fdkdx)'));
end
if ~isempty(d2fdxdk)
    d2fdxdk  = strtrim(evalc('disp(d2fdxdk)'));
end
if ~isempty(d3fdx3)
    d3fdx3   = strtrim(evalc('disp(d3fdx3)'));
end
if ~isempty(d3fdkdx2)
    d3fdkdx2 = strtrim(evalc('disp(d3fdkdx2)'));
end

% Remove spaces *too slow and not worth it
% f       = regexprep(f, ' ', '');
% dfdx    = regexprep(dfdx, ' ', '');
% dfdu    = regexprep(dfdu, ' ', '');
% dfdk    = regexprep(dfdk, ' ', '');
% d2fdx2  = regexprep(d2fdx2, ' ', '');
% d2fdk2  = regexprep(d2fdk2, ' ', '');
% d2fdkdx = regexprep(d2fdkdx, ' ', '');
% d2fdxdk = regexprep(d2fdxdk, ' ', '');

% Replace species names with vector index names
if opts.Verbose; fprintf('Replacing names with vectorized variables...\n'); end
if opts.Verbose; fprintf('   species...\n');end
for i = 1:nX
    name = sprintf('x(%d)', i);
    f        = regexprep(f, Xnicestrs{i}, name, 0);
    dfdx     = regexprep(dfdx, Xnicestrs{i}, name, 0);
    dfdu     = regexprep(dfdu, Xnicestrs{i}, name, 0);
    dfdk     = regexprep(dfdk, Xnicestrs{i}, name, 0);
    d2fdx2   = regexprep(d2fdx2, Xnicestrs{i}, name, 0);
    d2fdk2   = regexprep(d2fdk2, Xnicestrs{i}, name, 0);
    d2fdkdx  = regexprep(d2fdkdx, Xnicestrs{i}, name, 0);
    d2fdxdk  = regexprep(d2fdxdk, Xnicestrs{i}, name, 0);
    d3fdx3   = regexprep(d3fdx3, Xnicestrs{i}, name, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, Xnicestrs{i}, name, 0);
end

% Replace input namus with vector index names
if opts.Verbose; fprintf('   input...\n'); end
for i = 1:nU
    name = sprintf('u(%d)', i);
    f        = regexprep(f, Unicestr{i}, name, 0);
    dfdx     = regexprep(dfdx, Unicestr{i}, name, 0);
    dfdu     = regexprep(dfdu, Unicestr{i}, name, 0);
    dfdk     = regexprep(dfdk, Unicestr{i}, name, 0);
    d2fdx2   = regexprep(d2fdx2, Unicestr{i}, name, 0);
    d2fdk2   = regexprep(d2fdk2, Unicestr{i}, name, 0);
    d2fdkdx  = regexprep(d2fdkdx, Unicestr{i}, name, 0);
    d2fdxdk  = regexprep(d2fdxdk, Unicestr{i}, name, 0);
    d3fdx3   = regexprep(d3fdx3, Unicestr{i}, name, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, Unicestr{i}, name, 0);
end

% Replace parameters with vector index names
if opts.Verbose; fprintf('   parameters...\n');end
for i = 1:nK
    name = sprintf('k(%d)', i);
    f        = regexprep(f, Knicestrs{i}, name, 0);
    dfdx     = regexprep(dfdx, Knicestrs{i}, name, 0);
    dfdu     = regexprep(dfdu, Knicestrs{i}, name, 0);
    dfdk     = regexprep(dfdk, Knicestrs{i}, name, 0);
    d2fdx2   = regexprep(d2fdx2, Knicestrs{i}, name, 0);
    d2fdk2   = regexprep(d2fdk2, Knicestrs{i}, name, 0);
    d2fdkdx  = regexprep(d2fdkdx, Knicestrs{i}, name, 0);
    d2fdxdk  = regexprep(d2fdxdk, Knicestrs{i}, name, 0);
    d3fdx3   = regexprep(d3fdx3, Knicestrs{i}, name, 0);
    d3fdkdx2 = regexprep(d3fdkdx2, Knicestrs{i}, name, 0);
end

if opts.Verbose; fprintf('   ...done.\n'); end

%% Convert strings into function handles
if opts.Verbose; fprintf('Converting expressions into handles...'); end
f = eval(['@(t,x,u,k) [' f ']']);

if ~isempty(dfdx)
    dfdx = regexprep(dfdx, '[', ''); %remove extra "[" from front of lines
    dfdx = regexprep(dfdx, ']', ''); %and "]"
    dfdx = strtrim(dfdx);            %trim excess new lines from the end
    dfdx = eval(['@(t,x,u,k) sparse([' dfdx '])']);
end

if ~isempty(dfdu)
    dfdu = regexprep(dfdu, '[', '');
    dfdu = regexprep(dfdu, ']', '');
    dfdu = strtrim(dfdu);
    dfdu = eval(['@(t,x,u,k) sparse([' dfdu '])']);
end

if ~isempty(dfdk)
    dfdk = regexprep(dfdk, '[', '');
    dfdk = regexprep(dfdk, ']', '');
    dfdk = strtrim(dfdk);
    dfdk = eval(['@(t,x,u,k) sparse([' dfdk '])']);
end

if ~isempty(d2fdx2)
    d2fdx2 = regexprep(d2fdx2, '[', '');
    d2fdx2 = regexprep(d2fdx2, ']', '');
    d2fdx2 = strtrim(d2fdx2);
    d2fdx2 = eval(['@(t,x,u,k) sparse([' d2fdx2 '])']);
end

if ~isempty(d2fdk2)
    d2fdk2 = regexprep(d2fdk2, '[', '');
    d2fdk2 = regexprep(d2fdk2, ']', '');
    d2fdk2 = strtrim(d2fdk2);
    d2fdk2 = eval(['@(t,x,u,k) sparse([' d2fdk2 '])']);
end

if ~isempty(d2fdkdx)
    d2fdkdx = regexprep(d2fdkdx, '[', '');
    d2fdkdx = regexprep(d2fdkdx, ']', '');
    d2fdkdx = strtrim(d2fdkdx);
    d2fdkdx = eval(['@(t,x,u,k) sparse([' d2fdkdx '])']);
end

if ~isempty(d2fdxdk)
    d2fdxdk = regexprep(d2fdxdk, '[', '');
    d2fdxdk = regexprep(d2fdxdk, ']', '');
    d2fdxdk = strtrim(d2fdxdk);
    d2fdxdk = eval(['@(t,x,u,k) sparse([' d2fdxdk '])']);
end

if ~isempty(d3fdx3)
    d3fdx3 = regexprep(d3fdx3, '[', '');
    d3fdx3 = regexprep(d3fdx3, ']', '');
    d3fdx3 = strtrim(d3fdx3);
    d3fdx3 = eval(['@(t,x,u,k) sparse([' d3fdx3 '])']);
end

if ~isempty(d3fdkdx2)
    d3fdkdx2 = regexprep(d3fdkdx2, '[', '');
    d3fdkdx2 = regexprep(d3fdkdx2, ']', '');
    d3fdkdx2 = strtrim(d3fdkdx2);
    d3fdkdx2 = eval(['@(t,x,u,k) sparse([' d3fdkdx2 '])']);
end

if opts.Verbose; fprintf('done.\n'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 5: Build Psuedo Kronecker Model %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.Verbose; fprintf('Evaluating handles...'); end
%% Replace rate constants with their values
% Clear everything but variables necessary for updating
temp = opts.Verbose;
clear opts
opts.Verbose = temp;
clear Knicestrs Uind Unicestr Xnicestrs Ymembers defaultOpts i j nMem name temp

m.nX            = nX;
m.nP            = nK;
m.nU            = nU;
m.nY            = nY;

m.xNames        = Xnames.';
m.yNames        = Ynames.';
m.pNames        = Knames.';
m.uNames        = Unames.';

m.c             = C;
m.ic            = x0;
m.p             = k;

clear nX nP nU nY xNames yNames pNames uNames C

m.f         = @(t,x,u)f(t,x,u,k);
if ~isempty(dfdx)
    m.dfdx      = @(t,x,u)dfdx(t,x,u,k);
end
if ~isempty(dfdk)
    m.dfdp      = @(t,x,u)dfdk(t,x,u,k);
end
if ~isempty(dfdu)
    m.dfdu      = @(t,x,u)dfdu(t,x,u,k);
end
if ~isempty(d2fdx2)
    m.d2fdx2    = @(t,x,u)d2fdx2(t,x,u,k);
end
if ~isempty(d2fdk2)
    m.d2fdp2    = @(t,x,u)d2fdk2(t,x,u,k);
end
if ~isempty(d2fdkdx)
    m.d2fdpdx   = @(t,x,u)d2fdkdx(t,x,u,k);
end
if ~isempty(d2fdxdk)
    m.d2fdxdp   = @(t,x,u)d2fdxdk(t,x,u,k);
end
if ~isempty(d3fdx3)
    m.d3fdx3    = @(t,x,u)d3fdx3(t,x,u,k);
end
if ~isempty(d3fdkdx2)
    m.d3fdpdx2  = @(t,x,u)d3fdkdx2(t,x,u,k);
end
m.update    = @update;

if opts.Verbose; fprintf('done.\n'); end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Update function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function varargout = update(newk, newx0)
        
        % Apply inputs
        k = newk;
        x0 = newx0;
        
        m.ic        = x0;
        m.p         = k;
        
        m.f         = @(t,x,u)f(t,x,u,k);
        if ~isempty(dfdx)
            m.dfdx      = @(t,x,u)dfdx(t,x,u,k);
        end
        if ~isempty(dfdk)
            m.dfdp      = @(t,x,u)dfdk(t,x,u,k);
        end
        if ~isempty(dfdu)
            m.dfdu      = @(t,x,u)dfdu(t,x,u,k);
        end
        if ~isempty(d2fdx2)
            m.d2fdx2    = @(t,x,u)d2fdx2(t,x,u,k);
        end
        if ~isempty(d2fdk2)
            m.d2fdp2    = @(t,x,u)d2fdk2(t,x,u,k);
        end
        if ~isempty(d2fdkdx)
            m.d2fdpdx   = @(t,x,u)d2fdkdx(t,x,u,k);
        end
        if ~isempty(d2fdxdk)
            m.d2fdxdp   = @(t,x,u)d2fdxdk(t,x,u,k);
        end
        if ~isempty(d3fdx3)
            m.d3fdx3    = @(t,x,u)d3fdx3(t,x,u,k);
        end
        if ~isempty(d3fdkdx2)
            m.d3fdpdx2  = @(t,x,u)d3fdkdx2(t,x,u,k);
        end
        m.update    = @update;
        
        varargout{1} = m;
    end

end
