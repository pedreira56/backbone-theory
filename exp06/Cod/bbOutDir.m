function outDir = bbOutDir(callerFile)
%BBOUTDIR  Standard output folder for an experiment, created if absent.
%
%   OUTDIR = BBOUTDIR(mfilename('fullpath')) returns the folder where an
%   experiment driver should write its .mat files and figures, creating it
%   if it does not exist yet.
%
%   Layout convention:
%
%       ...\Exp1\Cod\      code (this file, bb*.m, E*_*.m, runE*.m)
%       ...\Exp1\Result\   output          <- returned by this function
%
%   The rule is: if the folder holding the calling script is named 'Cod',
%   the output goes to a sibling folder named 'Result'; otherwise it goes to
%   a 'Result' subfolder of the code folder itself.  Either way the output
%   stays inside the experiment tree and never lands one level too high,
%   which is the failure mode this function exists to prevent.
%
%   Example, inside runE1.m:
%       outDir = bbOutDir(mfilename('fullpath'));
%       save(fullfile(outDir,'E1.mat'),'-struct','out');
%
%   See also RUNE1, FULLFILE, MKDIR.

codeDir = fileparts(callerFile);
[parentDir, leaf] = fileparts(codeDir);

if strcmpi(leaf,'Cod')
    outDir = fullfile(parentDir,'Result');
else
    outDir = fullfile(codeDir,'Result');
end

if ~exist(outDir,'dir')
    [ok,msg] = mkdir(outDir);
    if ~ok
        error('bbOutDir:mkdirFailed', ...
            'Could not create the output folder "%s": %s', outDir, msg);
    end
end
end
