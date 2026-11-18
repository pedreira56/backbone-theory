function [matFile, outDir] = bbFindMat(stem, hint)
%BBFINDMAT  Locate the .mat file of an experiment, wherever this file sits.
%
%   [MATFILE,OUTDIR] = BBFINDMAT('E1') returns the full path of E1.mat and
%   the folder holding it, searching, in order:
%
%       <this folder>\..\Result\E1.mat     (when this folder is named Cod)
%       <this folder>\Result\E1.mat
%       <this folder>\E1.mat
%       <parent>\*\Result\E1.mat           (sibling experiment folders)
%       <grandparent>\*\Result\E1.mat
%
%   The last two are what make a post-processing script work from a folder
%   of its own: the results stay in  ...\Exp 1\Result\ ,  ...\Exp 3\Result\
%   and so on, and the script finds them without being moved.
%
%   [...] = BBFINDMAT(STEM,HINT) searches HINT first, which may be the .mat
%   file itself or the folder holding it.
%
%   The search never recurses: it looks one folder deep under each root, so
%   a large tree costs nothing.

here = fileparts(mfilename('fullpath'));
cand = {};

if nargin > 1 && ~isempty(hint)
    if isfile(hint)
        cand{end+1} = hint;
    else
        cand{end+1} = fullfile(hint,[stem '.mat']);
        cand{end+1} = fullfile(hint,'Result',[stem '.mat']);
    end
end

[parent, leaf] = fileparts(here);
if strcmpi(leaf,'Cod')
    cand{end+1} = fullfile(parent,'Result',[stem '.mat']);
end
cand{end+1} = fullfile(here,'Result',[stem '.mat']);
cand{end+1} = fullfile(here,[stem '.mat']);

roots = {parent, fileparts(parent)};
for r = 1:numel(roots)
    if isempty(roots{r}) || ~isfolder(roots{r}), continue, end
    d = dir(roots{r});
    for i = 1:numel(d)
        if d(i).isdir && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..')
            cand{end+1} = fullfile(roots{r},d(i).name,'Result',[stem '.mat']); %#ok<AGROW>
            cand{end+1} = fullfile(roots{r},d(i).name,[stem '.mat']);          %#ok<AGROW>
        end
    end
end

for i = 1:numel(cand)
    if isfile(cand{i})
        matFile = cand{i};
        outDir  = fileparts(matFile);
        fprintf('  using %s\n',matFile);
        return
    end
end

error('bbFindMat:notFound', ...
    ['Could not find %s.mat.\n' ...
     'Looked under %s and its sibling folders.\n' ...
     'Pass the folder explicitly, for example:  %s(''G:\\...\\Exp 1\\Result'')'], ...
     stem, here, stem);
end
