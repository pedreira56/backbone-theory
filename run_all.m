function run_all(quick)
%RUN_ALL  Reproduce every numerical experiment of the paper, in order.
%
%   RUN_ALL()      full runs, as reported in the paper.  Hours: see README.
%   RUN_ALL(true)  quick pass over everything, for a smoke test.
%
%   Each experiment lives in its own folder with its own copy of the shared
%   helpers, because those helpers evolved while the work was done and each
%   experiment is kept with the version it ran with.  This driver therefore
%   adds one code folder at a time and removes it afterwards, so that no
%   experiment can pick up another one's version of bbParams or bbSample.
%
%   The order is a dependency order, not a numbering: runE8 must precede
%   baselinesE8, which reads E8.mat; runE6 must precede E9, E10, E10b and
%   E11, which read E6.mat; E12_nonGaussian must precede postE12_dispersao.
%
%   Failures do not stop the sequence.  Each step is reported as it ends and
%   a summary is printed at the end, so that one missing toolbox or one long
%   experiment interrupted by hand does not discard the rest of the run.
%
%   TWO THINGS TO KNOW ABOUT THE QUICK PASS
%   It OVERWRITES the stored Result files with smoke-test output, because each
%   driver writes to its own Result folder.  Under version control the stored
%   files are recovered with a checkout; outside it, keep a copy first.
%   It also SKIPS E10b and E11.  A quick run of runE6 produces a universe of
%   ten attributes instead of fourteen, and those two scripts are tied to the
%   full universe: E11 by the realized backbone it examines, E10b by needing
%   the pairing term of every subset of the complete enumeration.  Both run
%   normally in the full pass, and both now stop with an explicit message if
%   given an E6.mat from a quick run.

if nargin < 1 || isempty(quick), quick = false; end
root = fileparts(mfilename('fullpath'));

%   folder(s) to add        function              full args  quick args  skip if quick
S = {
  {'exp01'},               'runE1',               {},        {true},     false
  {'exp02'},               'runE2',               {},        {true},     false
  {'exp03'},               'runE3',               {},        {true},     false
  {'exp04'},               'runE4',               {},        {true},     false
  {'exp05'},               'runE5',               {},        {true},     false
  {'exp06'},               'runE6',               {},        {true},     false
  {'exp07'},               'runE7',               {},        {true},     false
  {'exp08'},               'runE8',               {},        {true},     false
  {'exp10_12'},            'baselinesE8',         {},        {},         false
  {'exp09'},               'E9_pilotSelection',   {},        {5},        false
  {'exp10_12'},            'E10_fullyEmpirical',  {},        {5},        false
  {'exp10_12'},            'E10b_greedyOnly',     {},        {},         true
  {'exp10_12'},            'E11_withinDesign',    {},        {},         true
  {'exp10_12'},            'E12_nonGaussian',     {},        {5},        false
  {'exp10_12'},            'postE4',              {},        {},         false
  {'exp10_12'},            'postE6',              {},        {},         false
  {'exp10_12'},            'postE12_dispersao',   {},        {},         false
  {'exp10_12','exp01'},    'figE1',               {},        {},         false
  {'exp10_12','exp03'},    'figE3',               {},        {},         false
};

nS = size(S,1);
status = cell(nS,1);
tAll = tic;
if quick
    fprintf(['\nPassagem rapida: sobrescreve os arquivos guardados em Result com saida ' ...
             'de teste,\ne salta o E10b e o E11, que exigem a enumeracao completa do ' ...
             'runE6.  Ver a ajuda.\n']);
end

for i = 1:nS
    folders = S{i,1};  fname = S{i,2};
    args = S{i, 3 + double(quick)};
    if quick && S{i,5}
        status{i} = sprintf('salta            %s  (exige a enumeracao completa do runE6)', fname);
        fprintf('\n======== %d/%d  %s: saltado na passagem rapida ========\n', i, nS, fname);
        continue
    end
    paths = cellfun(@(f) fullfile(root, f, 'Cod'), folders, 'UniformOutput', false);
    fprintf('\n======== %d/%d  %s ========\n', i, nS, fname);
    for p = 1:numel(paths), addpath(paths{p}); end
    t0 = tic;
    try
        f = str2func(fname);
        f(args{:});
        status{i} = sprintf('ok    %6.1f min  %s', toc(t0)/60, fname);
    catch err
        status{i} = sprintf('FALHA %6.1f min  %s: %s', toc(t0)/60, fname, err.message);
        fprintf(2, 'falhou: %s\n', err.message);
    end
    for p = 1:numel(paths), rmpath(paths{p}); end
end

fprintf('\n======== resumo (%.1f min no total) ========\n', toc(tAll)/60);
fprintf('%s\n', status{:});
end
