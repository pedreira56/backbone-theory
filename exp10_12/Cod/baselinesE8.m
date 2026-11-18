function T = baselinesE8(capHours, hint)
%BASELINESE8  What is lost by omitting one component of the selection criterion.
%
%   T = BASELINESE8() reads E8.mat (see BBFINDMAT) and compares the backbone
%   recommended by the full criterion against three rules, each omitting one
%   component, across the price grid of E8 at a reference bench capacity.
%   Writes E8_baselines.txt next to the .mat.  Runs in seconds: no simulation
%   is repeated, every number comes from what E8 already stored.
%   T = BASELINESE8(CAPHOURS) picks the capacity level (default: the 60-hour
%   level used in the paper, or the first available).
%   T = BASELINESE8(CAPHOURS,FOLDER) searches FOLDER first.
%
%   THE DECOMPOSITION THIS EXPLOITS
%   For each candidate, E8 stored the risk of value transfer at the donor
%   counts the capacity allows, and the floor of that risk with effectively
%   unlimited donors.  Hence
%       floor            = information: the irreducible conditional variance
%       risk - floor     = geometry: the price of finite donors
%       costB            = measurement cost
%   and the full criterion is risk + lambda*cost.  Dropping one term gives a
%   rule that ignores exactly one component:
%       omit geometry     minimize floor + lambda*cost
%       omit information  minimize (risk - floor) + lambda*cost
%       omit cost         minimize risk
%   Each rule's recommendation is then scored BY THE FULL CRITERION, so the
%   loss is reported in one currency: risk in units of outcome variance plus
%   price times cost.  A rule that ignores a component is not being asked to
%   optimize something else; it is being asked what its choice costs under the
%   objective the design actually faces.
%
%   Mutual information is deliberately absent.  It would need its own
%   nonparametric estimator in a heavy-tailed population, and in the
%   linear-Gaussian experiment of Section 9.8 it selected the same backbone as
%   residual-variance selection, which is the omit-geometry rule here.

if nargin < 1, capHours = []; end
if nargin < 2, hint = ''; end
addpath(fileparts(mfilename('fullpath')));
[matFile, outDir] = bbFindMat('E8', hint);
E8 = load(matFile);

caps = E8.caps(:).';
if isempty(capHours)
    k = find(abs(caps - 60) < 1e-9, 1);
    if isempty(k), k = 1; end
else
    [~, k] = min(abs(caps - capHours));
end
capUsed = caps(k);
risk = E8.risk(:,k);  flo = E8.floor(:);  cost = E8.costB(:);
card = E8.card(:);    subsets = E8.subsets;
names = E8.meas.name;
lam = E8.lam(:).';
fprintf('baselinesE8: capacity %g h per test; %d candidates; %d price levels\n', ...
        capUsed, numel(risk), numel(lam));

exc = risk - flo;                       % geometric component
nameOf = @(b) ternary(isempty(b), '—', strjoin(cellstr(names(b)), ', '));

rows = {};
for i = 1:numel(lam)
    L = lam(i);
    full = risk + L*cost;
    [vFull, kFull] = min(full);
    [~, kGeo]  = min(flo + L*cost);         % omits geometry
    [~, kInfo] = min(exc + L*cost);         % omits information
    [~, kCost] = min(risk);                 % omits cost
    for rule = 1:4
        switch rule
            case 1, kk = kFull;  lab = 'critério completo';
            case 2, kk = kGeo;   lab = 'omite a geometria';
            case 3, kk = kInfo;  lab = 'omite a informação';
            case 4, kk = kCost;  lab = 'omite o custo';
        end
        rows(end+1,:) = {L, lab, card(kk), nameOf(subsets{kk}), ...
                         full(kk), full(kk) - vFull, ...
                         100*(full(kk) - vFull)/vFull};               %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'preco','regra','cardinalidade','backbone','criterio_completo','perda','perda_pct'});
disp(' ');
disp('=== E8: what each omission costs, scored by the full criterion ===');
disp(T);

fid = fopen(fullfile(outDir,'E8_baselines.txt'),'w');
fprintf(fid,'E8 baselines - cost of omitting one component of the criterion\n\n');
fprintf(fid,'bench capacity: %g hours per destructive test\n', capUsed);
fprintf(fid,'full criterion: risk of value transfer + price * measurement cost\n');
fprintf(fid,'omit geometry:     minimize floor + price*cost   (ignores finite donors)\n');
fprintf(fid,'omit information:  minimize excess + price*cost   (ignores predictive power)\n');
fprintf(fid,'omit cost:         minimize risk                  (ignores cost and capacity)\n');
fprintf(fid,'every recommendation is scored by the FULL criterion, so losses are comparable\n\n');
fprintf(fid,'%s\n', evalc('disp(T)'));
fclose(fid);
fprintf('\nWrote E8_baselines.txt to:\n  %s\n', outDir);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
