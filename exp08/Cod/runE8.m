function out = runE8(quick)
%RUNE8  Driver for experiment E8 (a worked design problem, Section 8).
%
%   OUT = RUNE8()      full run   (use this for the paper)
%   OUT = RUNE8(true)  quick run  (smoke test, under a minute)
%
%   Layout: put this file and E8_design.m in  ...\Exp 8\Cod\  and results
%   appear in  ...\Exp 8\Result\ , created automatically.  No other file is
%   needed: the population of this experiment is its own.
%
%   Outputs: E8.mat, E8_fig.pdf/.png, E8_table.txt
%
%   Panel (a) separates the two effects that a larger backbone has under a
%   fixed bench-time budget: it lowers the irreducible part of the risk, and
%   it raises the part due to matching, because fewer specimens are tested.
%   Panel (b) shows every candidate at one capacity, against its monetary
%   cost, with the frontier of candidates that no cheaper one beats.
%   Panel (c) shows the recommended backbone as the price of a measurement
%   rises, for the three capacities.

if nargin < 1 || isempty(quick), quick = false; end

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = localOutDir(codeDir);
fprintf('[E8] code   : %s\n[E8] results: %s\n\n', codeDir, outDir);

cfg = struct();
if quick
    cfg.nPilot = 5e4; cfg.nRecv = 1500; cfg.nFloor = 6000;
    cfg.repRisk = 2;  cfg.repFloor = 2;
else
    cfg.nPilot = 3e5; cfg.nRecv = 4000; cfg.nFloor = 20000;
    cfg.repRisk = 5;  cfg.repFloor = 4;
end
cfg.caps = [15 60 300];

tStart = tic;
out = E8_design(cfg);
out.outDir  = outDir;
out.seconds = toc(tStart);
fprintf('\n[E8] elapsed: %.1f min\n', out.seconds/60);

nS = numel(out.subsets);  caps = out.caps;  cMid = 2;
nameOf = @(b) strjoin(out.meas.name(b),' + ');

% best candidate at each cardinality, middle capacity
pMax = max(out.card);
bestAtCard = zeros(pMax+1,1);
for b = 0:pMax
    idx = find(out.card == b);
    [~,k] = min(out.risk(idx,cMid));
    bestAtCard(b+1) = idx(k);
end
[~,bestOverall] = min(out.risk(:,cMid));

f = figure('Position',[80 80 1350 420]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

% ---- (a) the two effects of a larger backbone --------------------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
bb = (0:pMax).';
fl = out.floor(bestAtCard);
rk = out.risk(bestAtCard,cMid);
bar(ax, bb, [fl, rk-fl], 'stacked', 'EdgeColor',[.3 .3 .3],'LineWidth',0.4);
ch = get(ax,'Children');
set(ch(2),'FaceColor',[0.72 0.78 0.86],'DisplayName','irreducible part');
set(ch(1),'FaceColor',[0.93 0.69 0.53],'DisplayName','matching, at the capacity allowed');
plot(ax, bb, rk, 'o-','Color',[0.2 0.2 0.2],'MarkerSize',4.5,'LineWidth',1.3, ...
     'MarkerFaceColor','w','DisplayName','criterion');
set(ax,'FontSize',9);
xlabel(ax,'number of measurements in the backbone');
ylabel(ax,'summed risk, in units of outcome variance');
text(ax,-0.02,1.05,'(a)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','northeast','FontSize',8,'Box','off');

% ---- (b) every candidate against its monetary cost ---------------------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
scatter(ax, out.costB, out.risk(:,cMid), 16, [0.55 0.6 0.68],'filled', ...
        'MarkerFaceAlpha',0.55,'DisplayName','all candidates');
% the frontier: candidates that no cheaper candidate beats
[~,ord] = sort(out.costB);
front = []; best = inf;
for i = ord.'
    if out.risk(i,cMid) < best - 1e-12, best = out.risk(i,cMid); front(end+1) = i; end %#ok<AGROW>
end
plot(ax, out.costB(front), out.risk(front,cMid), 'o-','Color',[0.85 0.33 0.10], ...
     'MarkerSize',5,'LineWidth',1.3,'MarkerFaceColor','w','DisplayName','frontier');
plot(ax, out.costB(bestOverall), out.risk(bestOverall,cMid),'p', ...
     'Color',[0.00 0.45 0.74],'MarkerSize',13,'MarkerFaceColor',[0.00 0.45 0.74], ...
     'DisplayName','minimiser, no monetary cost');
set(ax,'FontSize',9);
xlabel(ax,'monetary cost of the backbone');
ylabel(ax,'summed risk, in units of outcome variance');
text(ax,-0.02,1.05,'(b)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','northeast','FontSize',8,'Box','off');

% ---- (c) the recommendation as the price of a measurement rises --------
ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
% The three curves coincide over wide stretches of the price axis, so they
% are drawn with decreasing width and distinct dashes: where they overlap,
% all three remain visible instead of the last one hiding the others.
lw = [2.6 1.7 1.1];  ls = {'-','--',':'};
for c = 1:numel(caps)
    stairs(ax, out.lam, out.bestCard(:,c), ls{min(c,3)}, 'Color',cols(c,:), ...
           'LineWidth',lw(min(c,3)), ...
           'DisplayName',sprintf('%d bench hours per test',caps(c)));
end
set(ax,'FontSize',9);
xlabel(ax,'price of one measurement, in units of outcome variance');
ylabel(ax,'measurements in the recommended backbone');
ylim(ax,[-0.4, pMax+0.4]);
text(ax,-0.02,1.05,'(c)','Units','normalized','FontSize',11,'FontWeight','bold');
legend(ax,'Location','northeast','FontSize',8,'Box','off');

exportgraphics(f, fullfile(outDir,'E8_fig.pdf'),'ContentType','vector');
exportgraphics(f, fullfile(outDir,'E8_fig.png'),'Resolution',200);

% ---- console and table --------------------------------------------------
L = {};
L{end+1} = 'E8 - a worked design problem with costs and capacity';
L{end+1} = '';
L{end+1} = 'Candidate measurements';
L{end+1} = sprintf('  %-14s %6s %8s','measurement','cost','minutes');
for k = 1:numel(out.meas.cost)
    L{end+1} = sprintf('  %-14s %6d %8d',out.meas.name{k},out.meas.cost(k),out.meas.time(k));
end
L{end+1} = '';
L{end+1} = sprintf('Destructive tests: %s', strjoin(out.tests.name,', '));
L{end+1} = sprintf('  bench minutes per specimen: %s', num2str(out.tests.time));
L{end+1} = '';
L{end+1} = sprintf('Best backbone at each cardinality, %d bench hours per test',caps(cMid));
L{end+1} = sprintf('  %2s %9s %9s %9s %6s %8s   %s','b','risk','floor','excess','cost','donors','backbone');
for b = 0:pMax
    s = bestAtCard(b+1);
    L{end+1} = sprintf('  %2d %9.4f %9.4f %9.4f %6d %8d   %s', ...
        b, out.risk(s,cMid), out.floor(s), out.risk(s,cMid)-out.floor(s), ...
        out.costB(s), out.nDon(s,1,cMid), nameOf(out.subsets{s}));
end
L{end+1} = '';
L{end+1} = sprintf('Five best candidates at %d bench hours per test',caps(cMid));
[~,ordR] = sort(out.risk(:,cMid));
for i = 1:5
    s = ordR(i);
    L{end+1} = sprintf('  %d. risk %7.4f  cost %2d  donors %4d/%4d/%4d  %s', ...
        i, out.risk(s,cMid), out.costB(s), out.nDon(s,1,cMid), ...
        out.nDon(s,2,cMid), out.nDon(s,3,cMid), nameOf(out.subsets{s}));
end
L{end+1} = '';
L{end+1} = 'Recommendation at each capacity, with no monetary cost';
for c = 1:numel(caps)
    [~,s] = min(out.risk(:,c));
    L{end+1} = sprintf('  %3d bench hours: risk %7.4f  donors %4d/%4d/%4d  %s', ...
        caps(c), out.risk(s,c), out.nDon(s,1,c), out.nDon(s,2,c), out.nDon(s,3,c), ...
        nameOf(out.subsets{s}));
end
L{end+1} = '';
L{end+1} = 'Recommendation as the price of a measurement rises';
for c = 1:numel(caps)
    L{end+1} = sprintf('  %3d bench hours:',caps(c));
    prev = -1;
    for i = 1:numel(out.lam)
        s = out.bestIdx(i,c);
        if s ~= prev
            L{end+1} = sprintf('     price >= %5.3f : %s', out.lam(i), nameOf(out.subsets{s}));
            prev = s;
        end
    end
end
L{end+1} = '';
L{end+1} = sprintf('Elapsed: %.1f min',out.seconds/60);

fprintf('\n'); fprintf('%s\n',L{:});
fid = fopen(fullfile(outDir,'E8_table.txt'),'w');
fprintf(fid,'%s\n',L{:});
fclose(fid);

save(fullfile(outDir,'E8.mat'),'-struct','out');
fprintf('\nWrote E8.mat, E8_fig.pdf, E8_fig.png, E8_table.txt to:\n  %s\n', outDir);
end

% =======================================================================
function outDir = localOutDir(codeDir)
%LOCALOUTDIR  Sibling 'Result' folder when the code sits in 'Cod'.
[parentDir, leaf] = fileparts(codeDir);
if strcmpi(leaf,'Cod'), outDir = fullfile(parentDir,'Result');
else,                   outDir = fullfile(codeDir,'Result');
end
if ~exist(outDir,'dir')
    [ok,msg] = mkdir(outDir);
    if ~ok, error('runE8:mkdirFailed','Could not create "%s": %s',outDir,msg); end
end
end
